# UnPaSt

UnPaSt is a novel method for identification of differentially expressed biclusters in gene expression matrix. It searches for gene sets up- or down-regulated in subsets of samples:

![alt text](./poster/DESMOND2_steps2.png)


Webserver: https://unpast.zbh.uni-hamburg.de/


## Requirements:
<pre>
Python:
    fisher==0.1.9
    jenkspy==0.2.0
    pandas==1.4.2
    python-louvain==0.15
    matplotlib-venn==0.11.6
    numba==0.51.2
    numpy==1.22.3
    scikit-learn==0.23.1
    scikit-network==0.24.0
    scipy==1.7.1
    statsmodels==0.13.2

R:
    WGCNA==1.70-3
</pre>

## Examples
* UnPaSt requires a tab-separated file with standardized expressions of genes (or transcripts) in rows, and samples in columns. Gene and sample names must be unique. 
* A subset of 200 randomly chosen samples from TCGA-BRCA and UnPaSt output:
[test data](https://drive.google.com/file/d/1GXR_1ErIPtQkEOxE66at0uqQN76qNG7a/view?usp=sharing)

<pre>
# running UnPaSt with default parameters and example data
python run_unpast.py --exprs TCGA_200.exprs_z.tsv --basename TCGA_200_results

# with different binarization and clustering methods
python run_unpast.py --exprs TCGA_200.exprs_z.tsv --basename results --binarization ward --clustering WGCNA

# help
python run_unpast.py -h
</pre>

## Outputs
* \<basename\>.bin=[GMM|Jenks],clust=[Louvain|WGCNA|DESMOND].biclusters.tsv - a .tsv table with found biclsuters, where 
    - avgSNR is average SNR over all genes in the biclusters
    - columns "n_genes" and "n_samples" provide the numbers of genes and samples, respectively 
    - "gene","sample" contain gene and sample names respectively
    - "gene_indexes" and  "sample_indexes" - 0-based gene and sample indexes in the input matrix.
* binarized expressions, background distributions of SNR for each bicluster size and binarization statistics [if clustering is WGCNA,  or  '--save_binary' flag is added]
* modules found by WGCNA [if clustering is WGCNA]

## About 
UnPaSt is an unconstrained version of DESMOND method ([repository](https://github.com/ozolotareva/DESMOND), [publication](https://academic.oup.com/bioinformatics/article/37/12/1691/6039116?login=true))

Major modifications:
 * it does not require the network of gene interactions 
 * UnPaSt clusters individual genes instead of gene pairs
 * uses Gaussian mixture models or Jenks method for binarization of individual gene expressions
 * SNR threshold is authomatically determined; it depends on bicluster size in samples and user-defined p-value cutoff
 
## License
Free for non-for-profit use. For commercial use please contact the developers. 

### Poster CDCS workshop'22
![./poster/DESMOND2_poster_v5.png](./poster/DESMOND2_poster_v5.png)
### Poster ISMB and MCCMB'21
![./poster/DESMOND2.pdf](./poster/DESMOND2.png)
