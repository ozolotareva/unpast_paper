#!/usr/bin/bash

readonly a=0.05 # 0.5, 0.1, 0.05
readonly q=0.1 # 0.5, 0.25, 0.1
readonly b=1.0
readonly p_val=0.01
readonly network=$HOME/Downloads/unpast_trans/data/bicon_network.tsv
readonly base_dir=$HOME/Downloads/desmod_run/D1
readonly data_dir=$HOME/Downloads/unpast_trans/data
readonly script_path=$HOME/Dropbox/Doutorado/Coop/hamburgo/DESMOND

declare -A datasets=([GDC]="TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv" [Mbr]="METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv")

for i in "${!datasets[@]}"; do
    # mkdir $base_dir/$i
    python3 $script_path/DESMOND.py --exprs $data_dir/${datasets[$i]} \
        --network $network --basename $i --out_dir $base_dir/$i \
        --alpha $a --p_val $p_val -q $q --direction UP \
        --verbose >$base_dir/$i/$i.UP.LOG 2>$base_dir/$i/$i.UP.ERR &

    python3 $script_path/DESMOND.py --exprs $data_dir/${datasets[$i]} \
        --network $network --basename $i --out_dir $base_dir/$i --alpha $a \
        --p_val $p_val -q $q --direction DOWN \
        --verbose >$base_dir/$i/$i.DOWN.LOG 2>$base_dir/$i/$i.DOWN.ERR

    python3 $script_path/post-processing.py \
        --up $base_dir/$i/$i.'alpha='$a',beta_K='$b',direction=UP,p_val='$p_val',q='$q.biclusters.tsv \
        --down $base_dir/$i/$i.'alpha='$a',beta_K='$b',direction=DOWN,p_val='$p_val',q='$q.biclusters.tsv \
        --exprs $data_dir/${datasets[$i]} \
        --network $network --SNR_file $base_dir/$i/$i',q='$q'.SNR_threshold.txt' \
        --out $base_dir/$i/$i.'alpha='$a',beta_K='$b',p_val='$p_val',q='$q.biclusters.permutations.tsv \
        --verbose >$base_dir/$i/$i.permutations.LOG 2>$base_dir/$i/$i.permutations.ERR

done

echo "done"
