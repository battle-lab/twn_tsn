#!/bin/sh
twn_directory=$(cd $(dirname "$0") && pwd -P);
settings_fn="$twn_directory/settings.sh";
source $settings_fn;


### check matlab
echo "========== checking matlab installations... ==========";
cd "$twn_directory";
$matlab_path -nodisplay -nosplash -singleCompThread -r "try, check_prerequisites $isoform_annotation  $quic_directory, catch err, disp(err), end, quit";

### check R
echo "========== checking R installations... ==========";
$rscript_path check_prerequisites.R -gene_annot "$gene_annotation" -trans_annot "$isoform_annotation" -conflict "$cross_mappable_genes" -overlap "$positional_overlap";

echo "========== prerquisite checking done: please see messages above. ==========";
