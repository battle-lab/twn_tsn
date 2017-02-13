library(argparser)
library(data.table)

args <- arg_parser('program')
args <- add_argument(args, '-gene_annot',
                     help='gene annotation file',
                     default='data/gene_annot.txt')
args <- add_argument(args, '-trans_annot',
                     help='gene annotation file',
                     default='data/transcript_annot.txt')
args <- add_argument(args, '-conflict',
                     help='conflicting genes file',
                     default='data/cross_mappable_genes.txt')
args <- add_argument(args, '-overlap',
                     help='positional overlap file',
                     default='data/positional_overlap.txt')
args <- add_argument(args, '-te',
                     help='total expression file',
                     default='test/TE_sample1.txt')
args <- add_argument(args, '-ir',
                     help='isoform ratio file',
                     default='test/IR_sample1.txt')
args <- add_argument(args, '-o',
                     help='output satus file',
                     default='test/data_check_status.txt')

argv <- parse_args(args)
gene_annot_fn <- argv$gene_annot
trans_annot_fn <- argv$trans_annot
conflict_fn <- argv$conflict
positional_overlap_fn <- argv$overlap
te_fn <- argv$te
ir_fn <- argv$ir
out_fn <- argv$o

# checking fails by default
cat(0, file=out_fn)
  
# check if gene annotation file is ok
gene_annot <- tryCatch(fread(input = gene_annot_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F),
                             error = function(e){stop(paste0('Error: could not read the gene annotation file - ', e$message))})
if(length( intersect(colnames(gene_annot), c('gene_id', 'ensembl_gene_id'))) != 2)
    stop('gene annotation file is not properly formatted')
if(length(unique(gene_annot$gene_id)) != nrow(gene_annot))
  stop('gene ids must be unique in the gene annotation file')
rownames(gene_annot) <- gene_annot$gene_id

# check if gene annotation file is ok
trans_annot <- tryCatch(fread(input = trans_annot_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F),
                       error = function(e){stop(paste0('Error: could not read the isoform annotation file - ', e$message))})
if(length( intersect(colnames(trans_annot), c('transcript_id', 'gene_id', 'ensembl_gene_id'))) != 3)
  stop('isoform annotation file is not properly formatted')
if(length(unique(trans_annot$transcript_id)) != nrow(trans_annot))
  stop('transcript ids must be unique in the isoform annotation file')
rownames(trans_annot) <- trans_annot$transcript_id

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

# check if every gene has an entry in the gene annotation file
te_data <- tryCatch(fread(input = te_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F),
                    error = function(e){stop(paste0('Error: could not read the total expression data file - ', e$message))})
rownames(te_data) <- te_data[,1]
te_data <- te_data[,2:ncol(te_data)]
te_genes <- colnames(te_data)
if(length(intersect(te_genes, gene_annot$gene_id)) != length(te_genes))
  stop('not all genes in the total expression data file are present in the gene annotation file.')

# check if every transcript has an entry in the transcript annotation file
ir_data <- tryCatch(fread(input = ir_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F),
                    error = function(e){stop(paste0('Error: could not read the isoform ratio data file - ', e$message))})
rownames(ir_data) <- ir_data[,1]
ir_data <- ir_data[,2:ncol(ir_data)]
ir_transcripts <- colnames(ir_data)
if(length(intersect(ir_transcripts, trans_annot$transcript_id)) != length(ir_transcripts))
  stop('not all isoforms in the isoform ratio data file are present in the transcript annotation file.')

# check if every transcript belong to a gene with an entry in the gene annotatio file
ir_genes <- unique(trans_annot[ir_transcripts, 'gene_id'])
if(length(intersect(ir_genes, gene_annot$gene_id)) != length(ir_genes))
  stop('the gene of a transcript must have an entry in the gene annotation file.')

# check if there is any missing entry
te_missing = sum(!is.finite(as.matrix(te_data)))
if (te_missing > 0)
  stop('missing value in total expression data.')

ir_missing = sum(!is.finite(as.matrix(ir_data)))
if (ir_missing > 0)
  stop('missing value in isoform ratio data.')


# check if the data are normally distributed, otherwise give warning.


# data checking: passed
cat(1, file=out_fn)
