#!/bin/sh

### read inputs
if [ ! $# -eq 8 ]; then
  # there must have 8 parameters
  echo "invalid number of parameters: $#\nrequired parameters:\n 1) total expression file (string)\n 2) isoform ratio file (string)\n 3) output file prefix (string)\n 4) lambda_tt (float)\n 5) lambda_ti (float)\n 6) lambda_ii (float)\n 7) lambda_d (float)\n 8) lambda_s (float)";
  exit 1;
fi

te_fn=$1;
ir_fn=$2;
out_fn_prefix=$3;
l_tt=$4;
l_ti=$5;
l_ii=$6;
l_d=$7;
l_s=$8;

### load settings
twn_directory=$(cd $(dirname "$0") && pwd -P);
settings_fn="$twn_directory/settings.sh";
source $settings_fn;

### check if input files are ok

### run quic
cd "$twn_directory";
$matlab_path -nodisplay -nosplash -singleCompThread -r "try, twn $te_fn $ir_fn $out_fn_prefix  $l_tt $l_ti $l_ii $l_d $l_s $n_iteration $threshold 1 $standardize_data $isoform_annotation $quic_directory, catch err, disp(err), end, quit";

### convert matrix market file to readable network format
$python_path matrix_market2txt.py -in "$out_fn_prefix.mm" -expr $te_fn -iso $ir_fn -o "$out_fn_prefix.txt";

### remove conflicting and overlapped edges
$rscript_path remove_conflicting_edges.R -net "$out_fn_prefix.txt"  -gene_annot "$gene_annotation" -trans_annot "$isoform_annotation" -conflict "$mappability_conflict" -overlap "$positional_overlap" -o "$out_fn_prefix.final.txt";
