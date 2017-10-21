this folder should contain TSN codes and data.

# Tissue-Specific Network (TSN)

A Tissue Specific Network (TSN) is an undirected network constructed over total gene expressions using [Bicmix](https://www.cs.princeton.edu/~bee/software.html). It contains an edges between two nodes if the nodes are co-expressed only in the tissue of interest. Codes to reconstruct a TSN reside in the _tsn_ folder, here referred as _tsn directory_. 

### Installation

All the scripts were written and run in Linux environment with R 3.3.1.

1. This assumes that you have already run Bicmix. Bicmix code is available here:
	https://www.cs.princeton.edu/~bee/software.html
	A sample batch perl script for creating scripts for multiple bicmix jobs is below:

```
##### START SAMPLE PERL SCRIPT

$heredoc2 = <<END;

bicluster_mixture_simul_up  --y \$data_i --nf \$fac --sep tab --out \${dir}/\$result --interval \${interval} 
END

my $fac=1000;
my $interval=150;
my $dir="name of output directory"; #TODO: EDIT 
`mkdir $dir`;
my $data_i =("name of expression data file"); #TODO: EDIT

  for(my $iseed=1;$iseed<40;$iseed++){ #TODO: Change iseed limit based on how many bicmix runs are desired
        my $result="${data_i}_fac${fac}_interval${interval}_seed${iseed}"; #Recommended dir names to be compatible later
        `mkdir $dir/$result/`;
        my $file="${result}.sh";
        print "$file\n";
        open FILE,">$dir/$file" or die;
        print FILE "#!/bin/bash\n\n";
        print FILE "data_i=$data_i\n";
        print FILE "fac=$fac\n";
        print FILE "interval=$interval\n";
        print FILE "dir=$dir\n";
        print FILE "result=$result\n";
        print FILE $heredoc2;
        close FILE;
        #`./$dir/$file`; TODO: These shell scripts all need to be run
  }

####### END SAMPLE PERL SCRIPT -remember to run all the shell scripts you created in $dir
```
2. argparser package must be installed in R.


### Files and Formats


* Output Directory (-out): Path to directory in which results will be written
* File with Gene Names (Row labels) for the expression matrix (-gn): Path to file with one ensembl gene label per line, in the same order as in the expression matrix 
* File with Sample Size for each tissue used (-ss): Path to file with each line containing the number of samples from that tissue. The order of tissues is the order in which they appear in the expression matrix
* File with Covariate Matrix (-cov): Path to tab delimited file with n rows and x columns, where n is the number of samples in the expression matrix and x is the number of covariates. The matrix is filled in with the corresponding value of the covariate for each sample. There are no column or row names in this file. NOTE: recover_TSN.R will attempt to recover networks for each of these covariates. If you have more covariates than you are interested in getting specific networks for, you should delete those columns of the matrix.
* File with covariate names (-cn): Path to file with x rows, where each row lists the label of the corresponding column in the covariate matrix file. In other words, this should have the column names of the -cov file, with one label per row. NOTE: If you have edited the covariate matrix, remember to also edit the covariate names file.
* Directory with Bicmix outputs (-rd): Path to directory of bicmix results. This directory should contain one directory for each run of bicmix. The path to this directory should also include the name of the subdirectories (results for each run of bicmix) up to the unique number at the end. For example, in the sample perl file above we named these subdirectories ${data_i}_fac${fac}_interval${interval}_seed${iseed}, so this path would be "/path to larger results directory/${data_i}_fac${fac}_interval${interval}_seed". (obviously those are perl variables so the values would be there instead)
* Iteration (-it): Iteration number (ie 300) at which to use the results from Bicmix, since each run will run at its own pace and most likely finish at a different iteration.
* Number of runs of bicmix (-nr): This assumes that the directories partially named by the -rd path are numbered in order from 1 to whatever number you input here.
* Duplication Threshold (-thresh): Default is 4 -- meaning that edges in the TSNs are required to have shown up in at least (100/4)=25% of the bicmix runs which identified significant tissue-specific edges.
* GeneNet probability that an edge is nonzero (-gn): Default is 0.8


### How to reconstruct a TSN?
You have to run the script recover_TSN.R with the arguments listed above
  
After running that script, you will find several files in the given output directory. For each covariate where a specific network was recovered, there will be two files: both will start with the name of the covariate and one will end with _nodes.csv, the other with _edges.csv. Those two files will contain the information for the covariate specific networks.


## How to cite (it will be updated once the paper is published)
Saha A, Kim Y, Gewirtz ADH, Jo B, Gao C, McDowell IC, GTEx Consortium, Engelhardt BE, Battle A. 2016. Co-expression networks reveal the tissue-specific regulation of transcription and splicing. bioRxiv. http://biorxiv.org/content/early/2016/10/02/078741.abstract.

## Contact
Ashis Saha (ashis@jhu.edu)  
Ariel Gewirtz (gewirtz@princeton.edu)  


