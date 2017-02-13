all.pkg <- rownames(installed.packages())
if (! 'argparser' %in% all.pkg){
  print('argparser package was not installed in R.')
  stop('twn:argparsernotexist - argparser package was not installed in R.')
}

if (! 'data.table' %in% all.pkg){
  print('data.table package was not installed in R.')
  stop('twn:argparsernotexist - data.table package was not installed in R.')
}

library(argparser)
library(data.table)

args <- arg_parser('program')
args <- add_argument(args, '-gene_annot',
                     help='gene annotation file',
                     default='/scratch1/battle-fs1/ashis/progdata/gtex_hidden_factor/rsem/annot/gene_annot.txt')
args <- add_argument(args, '-trans_annot',
                     help='gene annotation file',
                     default='/scratch1/battle-fs1/ashis/progdata/gtex_hidden_factor/rsem/annot/transcript_annot.txt')
args <- add_argument(args, '-conflict',
                     help='conflicting genes file',
                     default='/scratch0/battle-fs1/annotation/mappability/pairwise_conflict.txt')
args <- add_argument(args, '-overlap',
                     help='positional overlap file',
                     default='/scratch0/battle-fs1/annotation/mappability/positional_overlap.txt')

argv <- parse_args(args)
gene_annot_fn <- argv$gene_annot
trans_annot_fn <- argv$trans_annot
conflict_fn <- argv$conflict
positional_overlap_fn <- argv$overlap

# check if gene annotation file is ok
gene_annot <- tryCatch(fread(input = gene_annot_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F), 
                             error = function(e){stop(paste0('Error: could not read the gene annotation file - ', e$message))})
if(length( intersect(colnames(gene_annot), c('gene_id', 'ensembl_gene_id'))) != 2)
    stop('gene annotation file is not properly formatted')
if(length(unique(gene_annot$gene_id)) != nrow(gene_annot))
  stop('gene ids must be unique in the gene annotation file')

# check if gene annotation file is ok
trans_annot <- tryCatch(fread(input = trans_annot_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F), 
                       error = function(e){stop(paste0('Error: could not read the isoform annotation file - ', e$message))})
if(length( intersect(colnames(trans_annot), c('transcript_id', 'gene_id', 'ensembl_gene_id'))) != 3)
  stop('isoform annotation file is not properly formatted')
if(length(unique(trans_annot$transcript_id)) != nrow(trans_annot))
  stop('transcript ids must be unique in the isoform annotation file')

# check if pairwise conflict file is ok
pairwise_conflicts <- tryCatch(fread(input = conflict_fn, sep = '\t', header = T, stringsAsFactors = F, colClasses = 'character', check.names = F, data.table = F), 
                        error = function(e){stop(paste0('Error: could not read the cross mappability file - ', e$message))})
if(length( intersect(colnames(pairwise_conflicts), c('gene1', 'gene2'))) != 2)
  stop('isoform annotation file is not properly formatted')

# check if pairwise conflict file is ok
positional_overlap <- tryCatch(fread(input = positional_overlap_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F), 
                               error = function(e){stop(paste0('Error: could not read the positional overlap file - ', e$message))})
if(length( intersect(colnames(positional_overlap), c('gene1', 'gene2'))) != 2)
  stop('isoform annotation file is not properly formatted')

cat('R installations are OK.\n')