settings_dir=$(cd $(dirname "$0") && pwd -P);
quic_directory="$settings_dir/QUIC";
gene_annotation="$settings_dir/data/gene_annot.txt"; 
isoform_annotation="$settings_dir/data/transcript_annot.txt";
positional_overlap="$settings_dir/data/positional_overlap.txt";
cross_mappable_genes="$settings_dir/data/cross_mappable_genes.txt";

n_iteration=50;
threshold=1e-4;
standardize_data=1;     # 1: standardize data, 0: do not standardize data
verbose=4;              # verbosity: 0 (lowest) to 4 (highest)

matlab_path="matlab";   # "matlab" to use default matlab
rscript_path="Rscript"; # "Rscript" to use default R
