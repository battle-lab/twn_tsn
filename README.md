# twn
Transcriptome-wide network

## Installation

All the scripts here were written and run in Linux environment with matlab 2013a, python 2.7.7, and R 3.1.1.

1. QUIC() function must be running in matlab. Please see instructions here: http://www.cs.utexas.edu/~sustik/QUIC/. For convenience, we included QUIC package in our repository. You may need to compile it following the instructions in README file inside QUIC package.
2. 3 packages must be installed in python: numpy, pandas, and argparse.
3. argparser package must be installed in R.


##Files and Formats

###Total Expression File
Tab delimited file containing corrected total expression data. Each row represent a sample and each column represents a gene. First row and first column contain gene ids and sample ids, respectively.

###Isoform Ratio File
Tab delimited file containing corrected isoform ratio data. Each row represent a sample and each column represents an isoform. First row and first column contain isoform ids and sample ids, respectively. Samples in both Total Expression and Isoform Ratio files have to be in the same order.

###Gene Annotation File
Tab delimited file with two columns: gene_id, and ensembl_gene_id. 

###Isoform Annotation File
Tab delimited file with three columns: transcript_id (isoform id),  gene_id, and ensembl_gene_id. 

###Positional Overlap File
Tab delimited file with two columns containing pair of genes (ensembl gene ids) with positional overlap in the genome.

###Cross Mappability File
Tab delimited file with two columns containing pair of genes (ensembl gene ids) with cross mappability (see the paper for details).


##Settings
settings.sh file contains necessary information to run twn.sh. You need to edit this file to customize your settings.


##Are the installations and the settings are OK?
To check if all pre-requisites have been successfully installed, and the setting file contains valid configuration, run the following command.

```
sh ./check_prerequisites.sh 
```

##How to construct a TWN?
You have to run the script twn.sh with the following arguments.  
  * Total expression data file path
  * Isoform ratio data file path
  * Output file prefix
  * 5 penalty parameters (lambda_tt, lambda_ti, lambda_ii, lambda_d, lambda_s). see the paper for details.

###Sample shell script code
```
te_fn="/Users/ashissaha/github/twn/test/TE_sample1.txt"  #total expression file
ir_fn="/Users/ashissaha/github/twn/test/IR_sample1.txt"  #isoform ratio file
out_fn="/Users/ashissaha/github/twn/results/sample1"     #output file prefix
l_tt=0.5    # penalty parameter
l_ti=0.4    # penalty parameter
l_ii=0.4    # penalty parameter
l_d=0       # penalty parameter
l_s=0.05    # penalty parameter

cd /Users/ashissaha/github/twn    # move to the twn source directory 
sh ./twn.sh $te_fn $ir_fn $out_fn $l_tt $l_ti $l_ii $l_d $l_s     # run twn
```

You will find a number of files with starting with the given output file prefix and the following suffixes:
.final.txt
.txt
.***

*.final.txt is a tab delimited file with four columns representing the constructed transcriptome-wide network.

The first two columns, containing either a total expression id or an isoform id, together represent an edge. The third column represents the type of the edge: 1 for an edge strictly between two total expressions, 2 for an edge between a total expression and an isoform, 3 for an edge strictly between two isoforms. The fourth column represents the edge weight.


## How to cite (it will be updated once the paper is published)
Saha A, Kim Y, Gewirtz ADH, Jo B, Gao C, McDowell IC, GTEx Consortium, Engelhardt BE, Battle A. 2016. Co-expression networks reveal the tissue-specific regulation of transcription and splicing. bioRxiv. http://biorxiv.org/content/early/2016/10/02/078741.abstract.
