#!/usr/bin/bash

readonly a=0.5 # 0.5, 0.1, 2.5
readonly q=0.5 # fixed
readonly b=1.0
readonly p_val=0.01 # fixed
readonly network=$HOME/Downloads/unpast_trans/data/bicon_network.tsv
readonly base_dir=$HOME/Downloads/desmod_run/D1
readonly data_dir=$HOME/Downloads/unpast_trans/data
readonly script_path=$HOME/Dropbox/Doutorado/Coop/hamburgo/DESMOND

declare -A datasets=([GDC]="TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv" [Mbr]="METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv")

for i in "${!datasets[@]}"; do
    for ((n=0;n<4;n++)); do
        if [ -d  $base_dir/$i/$n ]; then
            echo "folder exists"
        else
            mkdir -p $base_dir/$i/$n
        fi

        python3 $script_path/DESMOND.py --exprs $data_dir/${datasets[$i]} \
            --network $network --basename $i --out_dir $base_dir/$i/$n \
            --alpha $a --p_val $p_val -q $q --direction UP \
            --verbose >$base_dir/$i/$n/$i.UP.LOG 2>$base_dir/$i/$n/$i.UP.ERR &

        python3 $script_path/DESMOND.py --exprs $data_dir/${datasets[$i]} \
            --network $network --basename $i --out_dir $base_dir/$i/$n \
            --alpha $a --p_val $p_val -q $q --direction DOWN \
            --verbose >$base_dir/$i/$n/$i.DOWN.LOG 2>$base_dir/$i/$n/$i.DOWN.ERR

        python3 $script_path/post-processing.py \
            --up $base_dir/$i/$n/$i.'alpha='$a',beta_K='$b',direction=UP,p_val='$p_val',q='$q.biclusters.tsv \
            --down $base_dir/$i/$n/$i.'alpha='$a',beta_K='$b',direction=DOWN,p_val='$p_val',q='$q.biclusters.tsv \
            --exprs $data_dir/${datasets[$i]} \
            --network $network --SNR_file $base_dir/$i/$n/$i',q='$q'.SNR_threshold.txt' \
            --out $base_dir/$i/$n/$i.'alpha='$a',beta_K='$b',p_val='$p_val',q='$q.biclusters.permutations.tsv \
            --verbose >$base_dir/$i/$n/$i.permutations.LOG 2>$base_dir/$i/$n/$i.permutations.ERR

    done
done

echo "done"
