#!/bin/sh

### read inputs
if [ $# -lt 8 ] || [ $# -gt 9 ] ; then
  # there must have at 8 or 9 parameters
  echo "invalid number of parameters: $#\nrequired parameters:\n 1) total expression file (string)\n 2) isoform ratio file (string)\n 3) output file prefix (string)\n 4) lambda_tt (float)\n 5) lambda_ti (float)\n 6) lambda_ii (float)\n 7) lambda_d (float)\n 8) lambda_s (float)\n 9) settings_fn (string) [optional]";
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

twn_directory=$(cd $(dirname "$0") && pwd -P);

### load settings
if [ $# -ge 9 ]; then
  # settings parameter (9th) is optional
  settings_fn=$9
  # convert to absolute path
  case $settings_fn in
    /*) settings_fn=$settings_fn;;
    *) settings_fn=$PWD/$settings_fn;;
  esac
else
  # default settings file
  settings_fn="$twn_directory/settings.sh";
fi
source $settings_fn;

cd "$twn_directory";

### check if input files are ok
echo "checking data ...";
data_status_fn="${out_fn_prefix}_data_status.txt"
$rscript_path check_data.R -gene_annot "$gene_annotation" -trans_annot "$isoform_annotation" -conflict "$cross_mappable_genes" -overlap "$positional_overlap" -te "$te_fn" -ir "$ir_fn" -o $data_status_fn;
data_status=$(head -n 1 $data_status_fn)

if(( $data_status == '1' ))
then
    ### run quic
    echo "running QUIC ...";
    $matlab_path -nodisplay -nosplash -singleCompThread -r "try, twn $te_fn $ir_fn $out_fn_prefix  $l_tt $l_ti $l_ii $l_d $l_s $n_iteration $threshold 1 $standardize_data $isoform_annotation $quic_directory, catch err, disp(err), end, quit";

    ### filter unwanted edges
    echo "filtering unwanted edges ...";
    $rscript_path filter_edges.R -net "$out_fn_prefix.quic.txt"  -gene_annot "$gene_annotation" -trans_annot "$isoform_annotation" -conflict "$cross_mappable_genes" -overlap "$positional_overlap" -o "$out_fn_prefix.twn.txt";

    ### inform final status
    echo "Please find TWN here: $out_fn_prefix.twn.txt. Thanks!";
fi

