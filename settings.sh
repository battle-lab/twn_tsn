#settings_dir=$(cd $(dirname "$0") && pwd -P)/$(basename "$0");
settings_dir=$(cd $(dirname "$0") && pwd -P);
#echo $settings_dir;
directory_separator="/";  # for linux and mac: "/", for windows: "\"
quic_directory="$settings_dir/QUIC_1.2"; #"/Users/ashissaha/github/twn/QUIC_1.2";
gene_annotation="$settings_dir/data/gene_annot_unambiguous.txt"; 
isoform_annotation="$settings_dir/data/transcript_annot.txt";
positional_overlap="$settings_dir/data/positional_overlap.txt";
mappability_conflict="$settings_dir/data/avg_mappability_Exon_UTR.txt";

n_iteration=50;
threshold=1e-4;
standardize_data=1;     # 1: standardize data, 0: do not standardize data

matlab_path="matlab";   # "matlab" to use default matlab
python_path="python";   # "python" to use default python
rscript_path="Rscript"; # "Rscript" to use default R
