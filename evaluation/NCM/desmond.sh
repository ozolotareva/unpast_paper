a=0.5
b=1.0
q=0.5
p_val=0.01
network=$HOME/Downloads/unpast_trans/data/bicon_network.tsv
baseDir=$HOME/Downloads/desmod_run/D1
dataDir=$HOME/Downloads/unpast_trans/data

declare -A datasets=([GDC]="TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv" [Mbr]="METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv")

for i in "${!datasets[@]}"; do
    mkdir $baseDir/$i
    python3 $HOME/projects/DESMOND/DESMOND.py --exprs $dataDir/${datasets[$i]} \
        --network $network --basename $i --out_dir $baseDir/$i \
        --alpha $a --p_val $p_val -q $q --direction UP \
        --verbose >$baseDir/$i/$i.UP.LOG 2>$baseDir/$i/$i.UP.ERR

    python3 $HOME/projects/DESMOND/DESMOND.py --exprs $dataDir/${datasets[$i]} \
        --network $network --basename $i --out_dir $baseDir/$i --alpha $a \
        --p_val $p_val -q $q --direction DOWN \
        --verbose >$baseDir/$i/$i.DOWN.LOG 2>$baseDir/$i/$i.DOWN.ERR

    python3 $HOME/Dropbox/Doutorado/Coop/hamburgo/OlgaZ/DESMOND/post-processing.py --up $baseDir/$i/$i.'alpha='$a',beta_K='$b',direction=UP,p_val='$p_val',q='$q.biclusters.tsv \
        --down $baseDir/$i/$i.'alpha='$a',beta_K='$b',direction=DOWN,p_val='$p_val',q='$q.biclusters.tsv \
        --exprs $dataDir/${datasets[$i]} \
        --network $network -s $baseDir/$i/$i',q='$q'.SNR_threshold.txt' \
        --out $baseDir/$i/$i.'alpha='$a',beta_K='$b',p_val='$p_val',q='$q.biclusters.permutations.tsv \
        --verbose >$baseDir/$i/$i.permutations.LOG 2>$baseDir/$i/$i.permutations.ERR

done
