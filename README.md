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
Tab delimited file with two columns: _gene_id_, and _ensembl_gene_id_. 

###Isoform Annotation File
Tab delimited file with three columns: _transcript_id_ (isoform id),  _gene_id_, and _ensembl_gene_id_. 

###Positional Overlap File
Tab delimited file with two columns (_gene1_, _gene2_) containing pair of genes (ensembl gene ids) with positional overlap in the genome.

###Cross Mappability File
Tab delimited file with two columns (_gene1_, _gene2_) containing pair of genes (ensembl gene ids) with cross mappability (see the paper for details).

### Example Data
You may download example data from [here](https://drive.google.com/file/d/0B4XmrKDM9Pe3bmFaMWhqYXdmbjQ/view?usp=sharing). Please unzip the file and put inside the repository directory (in the directory where twn.sh file is). If you do not keep data in this directory, please update the settings file accordingly.

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
# parameters
twn_dir='/home/asaha6/github/twn'            # twn repository directory
te_fn="$twn_dir/data/demo/TE_demo.txt"       # total expression file
ir_fn="$twn_dir/data/demo/IR_demo.txt"       # isoform ratio file
out_fn_prefix="$twn_dir/demo/output_demo"    # output file prefix
l_tt=0.5      # penalty parameter
l_ti=0.4      # penalty parameter
l_ii=0.4      # penalty parameter
l_d=0         # penalty parameter
l_s=0.05      # penalty parameter

# move to the twn source directory 
cd $twn_dir   

# run twn
sh ./twn.sh $te_fn $ir_fn $out_fn_prefix $l_tt $l_ti $l_ii $l_d $l_s  
```

For convenience, a sample script has been provided in the demo folder to construct a TWN. After a successful run, you will find a number of files with starting with the given output file prefix (e.g., _*.twn.txt_, _*.iter_ etc.). 

_*.twn.txt_ is a tab delimited file with four columns representing the constructed transcriptome-wide network. Here, the first two columns, containing either a total expression id or an isoform id, together represent an edge. The third column represents the type of the edge: 1 for an edge strictly between two total expressions, 2 for an edge between a total expression and an isoform, 3 for an edge strictly between two isoforms. The fourth column represents the edge weight.  

Among other files, _*.quic.txt_ file represents a network similar to _*.twn.txt_, but it contains edges between overlapped or cross-mappable genes. _*.obj_ file contains the optimum objective value obtained from QUIC.  _*.iter_ file contains the number of iteration needed for the optimization. _*.time_ file cotains the time needed to finish the optimization.


## How to cite (it will be updated once the paper is published)
Saha A, Kim Y, Gewirtz ADH, Jo B, Gao C, McDowell IC, GTEx Consortium, Engelhardt BE, Battle A. 2016. Co-expression networks reveal the tissue-specific regulation of transcription and splicing. bioRxiv. http://biorxiv.org/content/early/2016/10/02/078741.abstract.

## Contact
Ashis Saha (ashis@jhu.edu)
