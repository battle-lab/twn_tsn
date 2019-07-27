
#This recovers tissue-specific networks for tissues from Bicmix output
##### Chuan Gao, edited and updated by A. Gewirtz #####

library("GeneNet")
library(Matrix)
library(reshape2)
library(ggplot2)


# argument parsing
args <- arg_parser("program");
args <- add_argument(args, "-out", help="output directory")
args <- add_argument(args, "-gn", help="file with gene names in order used in
                                     expression matrix, one per line")
args <- add_argument(args, "-cov", help="covariate matrix with samples as columns
                                    and covariates as rows")
args <- add_argument(args, "-cn", help="covariate names")
args <- add_argument(args, "-rd", help="directory with bicmix run outputs(each
                                     as a directory)")
args <- add_argument(args, "-it", help="iteration at which to use bicmix 
                                      output")
args <- add_argument(args, "-nr", help="number of runs of bicmix that you ran")
args <- add_argument(args, "-gn", help="GeneNet probability that an edge is 
                                      non-zero")
args <- add_argument(args, "-thresh", help="number of times you require an edge
                                       be duplicated to be included in network 
                                       = ceiling ( number of runs a covariate-
                                       specific factor was found in / THRESH)",
                                        default=4)

argv     <- parse_args(args)
cov_fn   <- argv$cov
cn_fn    <- argv$cn
dir      <- argv$rd
itr      <- argv$it
gn_fn    <- argv$gn
num_runs <- argv$nr
thresh   <- argv$thresh
outdir   <- argv$out
gn_prob  <- argv$gn


# read in the true data
geneNames <- read.table(gn_fn, as.is = T, stringsAsFactors = F)$V1
cov <- read.table(cov_fn, as.is = T)
covNames <- read.table(cn_fn, as.is = T)[,1]


# some basic parameters of the data
nG <- length(geneNames)      ## gene number
nS <- ncol(cov)              ## sample size
nCovTotal <- nrow(cov)       ## number of covariates


#################All the functions############################
##############################################################

## read in the results into a big list
read.data <- function(dir, sd, itr, nS, nG, nCovTotal){
  if(!file.exists(paste(dir, sd, "/LAM_", itr, sep = ""))){
    return("noFile")
  }else if(file.info(paste(dir, sd, "/LAM_", itr, sep = ""))$size == 0){
    return("noContent")
  }else{
    lamCov <- t(matrix(scan(paste(dir, sd, "/LAM_", itr, sep = "")), 
                ncol = nG))
    lam <- lamCov[1:nG,]
    z <- t(matrix(scan(paste(dir, sd, "/Z_", itr, sep = "")), ncol = 2))
    o <- t(matrix(scan(paste(dir, sd, "/O_", itr, sep = "")), ncol = 2))
    ex <- t(matrix(scan(paste(dir, sd, "/EX_", itr, sep = "")), nrow = nS))
    EXX <- scan(paste(dir, sd, "/EXX_", itr, sep = ""))
    EXX <- t(matrix(EXX, nrow = sqrt(length(EXX))))
    PSI <- scan(paste(dir, sd,"/PSI_",itr,sep = ""))
    return(list(lam = lam, z = z, ex = ex, o = o, EXX = EXX, PSI = PSI))
  }
}


# function to identify covariate specific factors: returns two list corresponding to 
#the two levels of the covariate; each list contains the index of the factors that 
#are specific to that level, and a matrix containing the number of samples that are 
#non-zeros for that factor level in the first row, and the number of samples that are 
#non-zeros for the other factor level in the second row 

find.fac.spec <- function(data, covSpec){
  lam <- data$lam; ex <- data$ex; z <- data$z; o <- data$o
  lamCount <- apply(lam, 2, function(x) {return(sum(x != 0))})
  sumCount <- apply(ex, 1, function(x) {
    count1 <- sum(x[covSpec == 1] != 0)
    count2 <- sum(x[covSpec == 0] != 0)
    return(c(count1, count2))
  })
  whichFac1 <- which(sumCount[1,] != 0 & sumCount[2,] == 0 & lamCount < 1000)
  Count1 <- sumCount[,whichFac1]
  whichFac2 <- which(sumCount[1,] == 0 & sumCount[2,] != 0 & lamCount < 1000)
  Count2 <- sumCount[,whichFac2]
  return(list(count1 <- list(whichFac = whichFac1, Count = Count1),
    count2 <- list(whichFac = whichFac2, Count = Count2))
  )
}


# Build the gene network that is specific to the covariate
get.network.all.component.spec <- function(data, facSpec, gn_prob){
  lam <- data$lam; ex <- data$ex; z <- data$z; o <- data$o; EXX <- data$EXX
  PSI <- data$PSI
  VXX <- EXX[facSpec, facSpec, drop = F] - ex[facSpec, , drop = F] %*%
       t(ex[facSpec, , drop = F])

  if(nrow(VXX) == 1){
    indexGene <- which((lam[, facSpec] != 0), arr.ind = T)
  }else{
    indexGene <- which((lam[,facSpec] != 0), arr.ind = T)[,1]
  }
  if(length(facSpec) == 1){
    indexGene <- matrix(indexGene, ncol = 1)
  }
  indexGene <- as.vector(indexGene)
  indexGene <- indexGene[!duplicated(indexGene)]
  indexGene <- sort(indexGene)
  if(length(facSpec) == 1 & length(indexGene) == 1){
    return(0)
  }

  LXXL <- lam[indexGene,facSpec] %*% VXX %*% t(lam[indexGene, facSpec]) + 
            diag(PSI[indexGene])
  LXXLi <- solve(LXXL)
  posrows <- which(diag(LXXL) > 0)
  precision <- cov2cor(LXXLi[posrows, posrows])
  nanrows <- which(is.na(precision))
  if(length(nanrows) > 1){
	 precision <- precision[-(nanrows),]
  }

  exitStatus <- tryCatch({
    arth.edges <- network.test.edges(precision, direct = FALSE, plot = FALSE)     
    if(nrow(precision) == 2){
      arth.edges <- t(arth.edges)
      arth.net <- arth.edges[arth.edges[,6] > gn_prob,]
    }else{
      arth.net <- arth.edges[arth.edges[,6] > gn_prob,]
    }
  },error = function(e){return("Error")})
  
  if(length(arth.net) == 0){
    return(all.net="no edges")
  }
  # If you would like to see the number of edges recovered, uncomment below
  # cat("number of edges ", nrow(arth.net), "\n")
  nname <- (1:nrow(lam))[indexGene]
  arth.net[,2] <- nname[arth.net[,2]]
  arth.net[,3] <- nname[arth.net[,3]]
  return(arth.net)
}

# vote on the duplicated edges using the ensemble method. 
# nDup is the number of times each edge is replicated. 
vote.network<- function(all.net, geneNames, nDup){
  all.net.out <- list()
  all.net.bak <- all.net
  
  edges <- data.frame(seed = all.net$sd, edges = paste(
            all.net$node1, "_", all.net$node2, sep = ""), stringsAsFactors = F)
  edgesND <- edges[!duplicated(edges$edges),]
  rownames(edgesND) <- edgesND$edges
  Dtable <- table(edges[,1:2])
  
  countDupEdges <- apply(Dtable, 2, function(x){return(sum(x > 0))})
  dupSeed <- apply(Dtable, 2, function(x){return(paste(names(which(x != 0)),
              collapse = ","))})
  
  Dindex <- (countDupEdges > nDup)
  if(sum(Dindex) == 0){return("noEdges")}
  goodEdges <- data.frame(edgesND[names(Dindex)[Dindex],], countEdges =
                countDupEdges[Dindex], DupSeed = dupSeed[Dindex], 
                stringsAsFactors = F)
  
  n <- nrow(goodEdges)
  node1 <- sapply(goodEdges[,2], function(x){return(strsplit(x, "_")[[1]][1])})
  node2 <- sapply(goodEdges[,2], function(x){return(strsplit(x, "_")[[1]][2])})
  all.net.out$edges <- data.frame(Source = node1, Target = node2, 
                        Weight = goodEdges[,3])
  
  nodeAll <- sort(as.numeric(as.character(c(node1, node2))))
  nodeAll <- nodeAll[!duplicated(nodeAll)]
  
  AllSd <- sapply(nodeAll, function(x){
    index <- grep(x, goodEdges[,2])
    AllSdi <- lapply(goodEdges$DupSeed[index], function(x){
                                                return(strsplit(x, ",")[[1]])})
    AllSdi <- unlist(AllSdi)
    AllSdi <- sort(as.numeric(AllSdi[!duplicated(AllSdi)]))
    AllSdi <- paste(AllSdi, collapse = ",")
    return(AllSdi)
  })
  
  nodeFull <- data.frame(Id = nodeAll, Label = geneNames[nodeAll], DupSeed = 
                AllSd, stringsAsFactors = F)
  all.net.out$nodes <- nodeFull
  return(all.net.out)
}


###########################Functions over###############################
#########################################################################

all.net <- rep(list(c()), nCovTotal)
sdFac <- rep(list(c()), nCovTotal)
names(all.net) <- covNames

seeds2use <- c(1:num_runs)

# find the tissue-specific factors from each run
for(i in 1:length(seeds2use)){
  isd <- seeds2use[i]   
  data <- read.data(dir, isd, itr, nS, nG, nCovTotal);
  lam <- data$lam; ex <- data$ex; z <- data$z; o <- data$o; EXX <- data$EXX
  PSI <- data$PSI
  lamCount <- apply(lam, 2, function(x){return(sum(x != 0))})
  xCount <- apply(ex, 1, function(x){return(sum(x != 0))})

  # iTis below controls which covariates you find tissue-specific networks for
  for(iTis in c(1:nCovTotal)){
    facSpec <- find.fac.spec(data, cov[,iTis]);
    if(length(facSpec[[1]]whichFac) == 0){
      next
    }
    arth.net <- get.network.all.component.spec(data, facSpec[[1]]whichFac,
                 gn_prob)
    if(length(arth.net) == 1){
      next
    }
    all.net[[iTis]] <- rbind(all.net[[iTis]], data.frame(arth.net, 
                        sd = rep(isd, nrow(arth.net))))

    sdFac[[iTis]] <- rbind(sdFac[[iTis]], data.frame(
                      sd = rep(isd, length(facSpec[[1]]whichFac)), 
                      fac = facSpec[[1]]whichFac, 
                      lamCount = lamCount[facSpec[[1]]whichFac],
                      xCount = xCount[facSpec[[1]]whichFac]))
  }
}

# select final edges to be in each covariate specific network
for(iTis in c(1:nCovTotal)){
  # get the number of runs of bicmix that produced a covariate-specific factor
  sdTisi <- all.net[[iTis]][,7] 
  sdTisi <- sdTisi[!duplicated(sdTisi)]
  nSd <- length(sdTisi)
  if(nSd == 0){
    next
  }
  # number of duplications required for an edge to be included
  nDup <- ceiling(nSd / thresh) 

  voteNetwork <- vote.network(all.net[[iTis]], geneNames, nDup)
  if(length(voteNetwork) == 1){
    next
  }

  write.table(voteNetwork$nodes[,1:2], paste(outdir, covNames[iTis], 
    "_spec_ensemble_nDup", nDup,"_nodes.csv", sep = ""), sep=";", quote = F,
    row.names = F, col.names = names(voteNetwork$nodes[,1:2]))

  write.table(voteNetwork$edges, paste(outdir, covNames[iTis], 
    "_spec_ensemble_nDup", nDup, "_edges.csv", sep = ""), sep = ";", quote = F,
    row.names = F, col.names = names(voteNetwork$edges))
}
