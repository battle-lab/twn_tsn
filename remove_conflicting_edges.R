### this script removes possibly conflicting edges due to alignment/mappability and positional overlap

library(argparser)
library(data.table)
library(ggplot2)
library(cowplot)

args <- arg_parser('program')
args <- add_argument(args, '-net',
                     help='network file',
                     default = '/scratch1/battle-fs1/ashis/results/gtex_twn/rsem/results_15000/quic_scale_free_combined/selected/diff_gene_net/WholeBlood.out.txt')
args <- add_argument(args, '-mnew',
                     help='new mappability file',
                     default = '/scratch0/battle-fs1/annotation/mappability/avg_mappability_Exon_UTR.txt')
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
args <- add_argument(args, '-o',
                     help='output network file',
                     default='results/artifact_removed_net.out.txt')

argv <- parse_args(args)
net_fn <- argv$net
mappability_new_fn <- argv$mnew
gene_annot_fn <- argv$gene_annot
trans_annot_fn <- argv$trans_annot
conflict_fn <- argv$conflict
positional_overlap_fn <- argv$overlap
out_fn <- argv$o

plt_prefix <- out_fn

#### read inputs
net <- fread(input = net_fn, sep = '\t', header = T, stringsAsFactors = F, colClasses = c('character', 'character', 'integer', 'numeric'), col.names = c('Name1', 'Name2', 'Edge.type', 'Edge.weight'), check.names = T, data.table = F)
#dim(net)
#head(net)

mappability_new <- fread(input = mappability_new_fn, sep = '\t', header = F, stringsAsFactors = F, colClasses = c('character', 'numeric'), col.names = c('gene', 'mappability'), data.table = F)
#dim(mappability_new)
#head(mappability_new)
mappability_new <- mappability_new[!is.na(mappability_new$mappability), ]  # removed NAs

gene_annot <- fread(input = gene_annot_fn, sep = '\t', header = T, stringsAsFactors = F, colClasses = c(chr='character'), check.names = F, data.table = F)
rownames(gene_annot) <- gene_annot$id
#dim(gene_annot)
#head(gene_annot)

trans_annot <- fread(input = trans_annot_fn, sep = '\t', header = T, stringsAsFactors = F, colClasses = c(chr='character'), check.names = F, data.table = F)
rownames(trans_annot) <- trans_annot$id
#dim(trans_annot)
#head(trans_annot)

pairwise_conflicts <- fread(input = conflict_fn, sep = '\t', header = F, stringsAsFactors = F, colClasses = 'character', col.names=c('gene1','gene2'), check.names = F, data.table = F)
#dim(pairwise_conflicts)
#head(pairwise_conflicts)

positional_overlap <- fread(input = positional_overlap_fn, sep = '\t', header = T, stringsAsFactors = F, check.names = F, data.table = F)
#dim(positional_overlap)
#head(positional_overlap)



##### get enesemble ids of genes in net
te_te_net <- net[net$Edge.type==1, ]
te_ir_net <- net[net$Edge.type==2, ]
ir_ir_net <- net[net$Edge.type==3, ]

all_genes <- unique(c(te_te_net$Name1, te_te_net$Name2, te_ir_net$Name1))
gene_ensgid_df <- merge(data.frame(sym=all_genes, stringsAsFactors = F), gene_annot, by.x='sym', by.y='sym', all.x=T,  all.y=F)

if(nrow(gene_ensgid_df) > length(all_genes)){
  # all geneids could not be identified unambiguously
  gene_counts <- tapply(gene_ensgid_df$sym, gene_ensgid_df$sym, length)
  ambiguous_genes <- names(which(gene_counts>1))
  warning(paste('these geneids could not be identified unambiguously: ', paste(ambiguous_genes, sep=',', collapse = ',')))
  # breaking ambiguity
  selected_idx <- tapply(1:nrow(gene_ensgid_df), gene_ensgid_df$sym, function(x) x[1])
  gene_ensgid_df <- gene_ensgid_df[selected_idx,]
}
rownames(gene_ensgid_df) <- gene_ensgid_df$sym

te_te_net['ensgid1'] <- gene_ensgid_df[te_te_net$Name1, 'id']
te_te_net['ensgid2'] <- gene_ensgid_df[te_te_net$Name2, 'id']
te_ir_net['ensgid1'] <- gene_ensgid_df[te_ir_net$Name1, 'id']
te_ir_net['ensgid2'] <- trans_annot[te_ir_net$Name2, 'geneid']
ir_ir_net['ensgid1'] <- trans_annot[ir_ir_net$Name1, 'geneid']
ir_ir_net['ensgid2'] <- trans_annot[ir_ir_net$Name2, 'geneid']

##### filter conflict data based on edges in the net
all_ensgids <- unique(c(te_te_net$ensgid1, te_te_net$ensgid2, 
                        te_ir_net$ensgid1, te_ir_net$ensgid2, 
                        ir_ir_net$ensgid1, ir_ir_net$ensgid2))

pairwise_conflicts <- pairwise_conflicts[(pairwise_conflicts$gene1 %in% all_ensgids) & (pairwise_conflicts$gene2 %in% all_ensgids), ]

##### filter out all conflicting edges (edge="gene1|gene2") - both pairwise mappability conflicts and positional overlap conflicts
conflict_edges <- paste(c(pairwise_conflicts$gene1, pairwise_conflicts$gene2, positional_overlap$gene1, positional_overlap$gene2) , c(pairwise_conflicts$gene2, pairwise_conflicts$gene1, positional_overlap$gene2, positional_overlap$gene1), sep='|')

te_te_edges <- paste(te_te_net$ensgid1, te_te_net$ensgid2, sep='|')
te_te_net_filtered <- te_te_net[!(te_te_edges %in% conflict_edges), ]

te_ir_edges <- paste(te_ir_net$ensgid1, te_ir_net$ensgid2, sep='|')
te_ir_net_filtered <- te_ir_net[!(te_ir_edges %in% conflict_edges), ]

ir_ir_edges <- paste(ir_ir_net$ensgid1, ir_ir_net$ensgid2, sep='|')
ir_ir_net_filtered <- ir_ir_net[!(ir_ir_edges %in% conflict_edges), ]

net_filtered <- rbind(te_te_net_filtered, te_ir_net_filtered, ir_ir_net_filtered)

##### save filtered net
write.table(net_filtered[,c('Name1', 'Name2', 'Edge.type', 'Edge.weight')], file=out_fn, col.names = c('Name1', 'Name2', 'Edge type', 'Edge weight'), row.names = F, quote = F, sep='\t')
