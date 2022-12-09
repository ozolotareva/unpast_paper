import os


MIN_SAMPLES_PER_COHORT = 10
SUBSAMPLE_FACTOR = 0.05
OUTPUT_FOLDER = '/local/DESMOND2_data/v6/preprocessed_v6/subsampled'

DATASETS = {
    'METABRIC': {
        'annotation': '/local/DESMOND2_data/v6/preprocessed_v6/METABRIC_1904.annotation_v6.tsv',
        'expression': '/local/DESMOND2_data/v6/preprocessed_v6/METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv',
        'subtypes': '/local/DESMOND2_data/v6/preprocessed_v6/METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv'
    },
    'TCGA': {
        'annotation': '/local/DESMOND2_data/v6/preprocessed_v6/TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv',
        'expression': '/local/DESMOND2_data/v6/preprocessed_v6/TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv',
        'subtypes': '/local/DESMOND2_data/v6/preprocessed_v6/TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv'
    }
}

GENE_MAPPING_FILE = '/local/DESMOND2_data/v6/preprocessed_v6/gene_id_mapping.tsv'

SUBSAMPLED_SUBTYPE_FILE = os.path.join('subtypes.tsv')
SUBSAMPLED_ANNOTATION_FILE = os.path.join('annotation.tsv')
SUBSAMPLED_EXPRESSION_FILE = os.path.join('expression.tsv')
SUBSAMPLED_SIZES_FILE = os.path.join('sizes.json')

UNPAST_BEST_PARAMS = {
    'METABRIC': {
        'bin_method': 'kmeans',
        'pval': .05,
        'clust_method': 'Louvain',
        'modularity': .7
    },
    'TCGA': {
        'bin_method': 'kmeans',
        'pval': .01,
        'clust_method': 'Louvain',
        'modularity': .8
    },
    'OVERALL': {
        'bin_method': 'kmeans',
        'pval': .01,
        'clust_method': 'WGCNA',
        'ds': 0
    }
}