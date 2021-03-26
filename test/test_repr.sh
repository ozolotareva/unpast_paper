#!/usr/bin/env bash
mkdir wd2
python3 ../run_desmond.py rand_matrix.tsv --out_dir wd2/ > /dev/null
d=$(diff wd/*1.5.biclusters.tsv wd2/*1.5.biclusters.tsv)
if [[ -z "$d" ]]; then
    echo "Code works the same"
else 
    echo "Something changed"
fi


