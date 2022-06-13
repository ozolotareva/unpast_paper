# DESMOND2

DESMOND2 is a novel method for identification of differentially expressed biclusters.
It is an unconstrained version of DESMOND: https://github.com/ozolotareva/DESMOND

Major modifications:
 * does not require the network 
 * DESMOND2 clusters individual genes instead of gene pairs
 * uses Gaussian mixture models for binarization of gene expressions
 * a gene may be assigned to multiple biclusters where it spent more than f time during the sampling phase
 * SNR threshold is authomatically determined based on bicluster size and user-defined p-value cutoff

## Requirements:
 * tested for Python 3.8.3
<pre>
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
</pre>
## Jupyter notebooks
* Easy-to-run random data example [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ozolotareva/DESMOND2/blob/main/random_data_example.ipynb)
* Real data example [DESMOND2_step_by_step.ipynb](DESMOND2_step_by_step.ipynb)


## Poster
![./poster/DESMOND2.pdf](./poster/DESMOND2.png)
