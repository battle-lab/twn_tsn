twn_dir='/home/asaha6/github/twn'
te_fn="$twn_dir/data/demo/TE_demo.txt"
ir_fn="$twn_dir/data/demo/IR_demo.txt"
out_fn="$twn_dir/demo/output_demo"
l_tt=0.5
l_ti=0.4
l_ii=0.4
l_d=0
l_s=0.05

cd $twn_dir
sh ./twn.sh $te_fn $ir_fn $out_fn $l_tt $l_ti $l_ii $l_d $l_s
