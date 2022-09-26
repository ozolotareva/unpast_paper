import pandas as pd

metabric_gt = pd.read_csv('/Users/fernando/Documents/Research/DESMOND2_data_simulated/real_data/METABRIC_1904.annotation_v6.tsv', sep='\t')
TCGA_gt = pd.read_csv('/Users/fernando/Documents/Research/DESMOND2_data_simulated/real_data/TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv', sep='\t')

metabric_gt.groupby('mol_subt').count()
TCGA_gt.groupby('PAM50').count()