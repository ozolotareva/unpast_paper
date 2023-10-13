# UnPaSt

UnPaSt is a novel method for identification of differentially expressed biclusters.

![alt text](./poster/DESMOND2_steps2.png)


## Requirements:
<pre>
Python:
    fisher==0.1.10
    jenkspy==0.2.0
    matplotlib-venn==0.11.6
    numba==0.55.2
    numpy==1.22.3
    scikit-learn==1.1.0
    scikit-network==0.25.0
    scipy==1.7.3
    statsmodels==0.13.2
    pandas==1.4.2
    python-louvain==0.15
    statsmodels==0.13.2

R:
    WGCNA>=1.70-3
    limma>=3.42.2
</pre>


## Installation
* UnPaSt can be installed using `pip`, `poetry`, or run using `Docker`, or as a script (see examples section). Follow the appropriate instructions below for your preferred method. You need to have R and **Python 3.8-3.10** installed.

1. Using **pip**: \
    To install the project using `pip`, first make sure you have `pip` installed on your system. If you haven't installed it already, you can find the installation instructions [here](https://pip.pypa.io/en/stable/installation/). \
    Once `pip` is installed, you can install UnPaSt by running the following command:

    ```bash
    pip install unpast
    ```
    Run it:
    ```bash
    run_unpast -h
    ```
    **Dependencies**. To use this package, you will need to have R and the [WGCNA library](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) installed. You can easily install these dependencies by running the following command after installing unpast:
    ```bash
    python -m unpast.install_r_dependencies

    # or you can install it directly
    R -e "install.packages('BiocManager'); BiocManager::install(c('WGCNA', 'limma'))"
    ```

2. Installation using **Poetry**: \
    To install the package using Poetry, first make sure you have Poetry installed, clone the repo and run:
    ```bash
    poetry add unpast
    ```
    Run it:
    ```bash
    poetry run run_unpast -h
    ```
    **Dependencies**. To use this package, you will need to have R and the [WGCNA library](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) installed. You can easily install these dependencies by running the following command after installing unpast:
    ```bash
    poetry run python -m unpast.install_r_dependencies

    # or you can install it directly
    R -e "install.packages('BiocManager'); BiocManager::install(c('WGCNA', 'limma'))"
    ```
3. Running with **Docker**: \
    You can also run the package using Docker. First, pull the Docker image:
    ```bash
    docker pull freddsle/unpast:latest
    ```
    Next, run the UnPaSt:
    ```bash
    docker run -v /your/data/path/:/user_data/ freddsle/unpast:latest --exprs /user_data/exprs.tsv --out_dir /user_data/out_dir/
    ```


## Examples
* UnPaSt requires a tab-separated file with features (e.g. genes) in rows, and samples in columns. Feature and sample names must be unique. 

<pre>

cd test;
mkdir -p results;

# running UnPaSt with default parameters and example data
python ../run_unpast.py --exprs scenario_B500.exprs.tsv.gz --basename results/scenario_B500

# with different binarization and clustering methods
python ../run_unpast.py --exprs scenario_B500.exprs.tsv.gz --basename results/scenario_B500 --binarization ward --clustering Louvain

# help
python run_unpast.py -h
</pre>

## Outputs
* \<basename\>.[parameters].biclusters.tsv - a .tsv table with found biclsuters, where 
    - the first line starts from '#' and stores parameters
    - each following line represents a bicluster
    - SNR column contains SNR of a bicluster 
    - columns "n_genes" and "n_samples" provide the numbers of genes and samples, respectively 
    - "gene","sample" contain gene and sample names respectively
    - "gene_indexes" and  "sample_indexes" - 0-based gene and sample indexes in the input matrix.
* binarized expressions, background distributions of SNR for each bicluster size and binarization statistics [if clustering is WGCNA,  or  '--save_binary' flag is added]

## About 
UnPaSt is an unconstrained version of DESMOND method ([repository](https://github.com/ozolotareva/DESMOND), [publication](https://academic.oup.com/bioinformatics/article/37/12/1691/6039116?login=true))

Major modifications:
 * it does not require the network of feature interactions 
 * UnPaSt clusters individual features instead of pairs of features
 * uses 2-means, hierarchicla clustering or GMM for binarization of individual gene expressions
 * SNR threshold for featuer selection is authomatically determined; it depends on bicluster size in samples and user-defined p-value cutoff
 
## License
Free for non-for-profit use. For commercial use please contact the developers. 
