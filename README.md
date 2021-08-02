# DESMOND2

DESMOND2 is a novel method for identification of differentially expressed biclusters.
It is an unconstrained version of DESMOND: https://github.com/ozolotareva/DESMOND

Major modifications:
 * does not require the network 
 * DESMOND2 clusters individual genes instead of gene pairs
 * uses Gaussian mixture models for binarization of gene expressions
 * a gene may be assigned to multiple biclusters where it spent more than f time during the sampling phase
 * SNR threshold is authomatically determined based on bicluster size and user-defined p-value cutoff



## Jupyter notebooks
* Easy-to-run random data example [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ozolotareva/DESMOND2/blob/main/random_data_example.ipynb)
* Real data example [DESMOND2_step_by_step.ipynb](DESMOND2_step_by_step.ipynb)