#!/bin/sh
twn_directory=$(cd $(dirname "$0") && pwd -P);
settings_fn="$twn_directory/settings.sh";
source $settings_fn;


### check matlab
echo "checking matlab installations..."
cd "$twn_directory";
$matlab_path -nodisplay -nosplash -singleCompThread -r "try, check_prerequisites $isoform_annotation  $quic_directory, catch err, disp(err), end, quit";

### check python
#$python_path matrix_market2txt.py -in "$out_fn_prefix.mm" -expr $te_fn -iso $ir_fn -o "$out_fn_prefix.txt";

### check R
#$rscript_path remove_conflicting_edges.R -net "$out_fn_prefix.txt"  -gene_annot "$gene_annotation" -trans_annot "$isoform_annotation" -conflict "$cross_mappable_genes" -overlap "$positional_overlap" -o "$out_fn_prefix.final.txt";
