'''
some analyses require graphical lasso network to be in a readable format like below:

Name1   Name2   Edge type   Edge weight
MNF1    PER1    1           0.0015

this script converts matrix market file to such format.
it also removes edges with zero weight.
'''

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', '--input',
                    help='input matrix market file.',
                    default='/scratch1/battle-fs1/ashis/results/gtex_twn/glasso_4000/Muscle-Skeletal_l0.2.out.mm')
parser.add_argument('-expr',
                    help='expression file used for creating the mm file. needed for gene id.',
                    default='/scratch0/battle-fs1/GTEx_Analysis_2015-01-12/processedData/data_used_for_twn/gene_rpkm_4000/Muscle-Skeletal.txt')
parser.add_argument('-iso',
                    help='isoform file used for creating the mm file. needed for isoform id.',
                    default='/scratch0/battle-fs1/GTEx_Analysis_2015-01-12/processedData/data_used_for_twn/isoform_ratio_4000/Muscle-Skeletal.txt')
parser.add_argument('-o', '--output',
                    help='output path to save the network file.',
                    default='results/network.txt')

args = parser.parse_args()

''' setting variables '''
in_fn = args.input
expr_fn = args.expr
iso_fn = args.iso
out_fn = args.output

# read inputs
mm_data = pd.read_table(in_fn, sep='\\s+', header=None, index_col=None, skiprows=3, engine='python')
expr_data = pd.read_table(expr_fn, sep='\t', header=0, index_col=0, nrows=1)
iso_data = pd.read_table(iso_fn, sep='\t', header=0, index_col=0, nrows=1)

# check if input files have right dimensions
genes = expr_data.columns.values.tolist()
isoforms = iso_data.columns.values.tolist()

# remove unnamed genes
genes = [g for g in genes if not g.startswith('Unnamed')]
isoforms = [i for i in isoforms if not i.startswith('Unnamed')]

# combine features
features = genes + isoforms

n_genes = len(genes)
n_isoforms = len(isoforms)

with open(in_fn) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith('%'):
            continue
        node_counts = line.split(' ')
        if int(node_counts[0]) != n_genes + n_isoforms:
            raise Exception('matrix market file is supposed to have exactly all the genes and isoforms.')
        break

###### convert mm file
# remove zero weights
non_zeros = mm_data.iloc[:,2].values != 0
mm_data = mm_data.loc[non_zeros,:]

# smaller index node comes first
node1s_idx = mm_data.iloc[:,:2].min(axis=1) - 1;
node2s_idx = mm_data.iloc[:,:2].max(axis=1) -1;

# get node names
node1s = [features[idx] for idx in node1s_idx]
node2s = [features[idx] for idx in node2s_idx]

# make a new data frame
df = pd.DataFrame()
df['Name1'] = node1s
df['Name2'] = node2s
df['Edge type'] = 3 # default iso-iso, set later
df['Edge weight'] = mm_data.iloc[:,2]


# Set edge type
is_node1_gene = node1s_idx < n_genes
is_node2_gene = node2s_idx < n_genes
gene_gene = is_node1_gene & is_node2_gene
gene_iso = is_node1_gene & np.invert(is_node2_gene)
df.loc[gene_gene, 'Edge type'] = 1
df.loc[gene_iso, 'Edge type'] = 2

# save file
df.to_csv(out_fn, sep='\t', header=True, index=False)
