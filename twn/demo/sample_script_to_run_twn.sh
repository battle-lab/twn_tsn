# parameters
twn_dir='/home/asaha6/github/twn_tsn/twn'    # twn directory
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
