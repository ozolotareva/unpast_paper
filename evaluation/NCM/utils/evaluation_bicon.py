from joblib import Parallel, delayed
from os.path import join
from os import listdir

import pandas as pd

from evaluation_desmond import (calculate_perfromance, make_ref_groups)


def prepar_results(pasta, n):

    file_name = [x for x in listdir(join(pasta, n)) if "results.csv" in x][0]
    df = pd.read_csv(join(pasta, n, file_name), index_col=0)
    df_neu = pd.DataFrame([[df.genes1.values, df.patients1.values], [
        df.genes2.values, df.patients2.values]],
        columns=["genes", "samples"])
    df_neu.genes = df_neu.genes.apply(
        lambda x: str(x[0]).split("|"))
    df_neu.samples = df_neu.samples.apply(
        lambda x: str(x[0]).split("|"))
    return df_neu


classifications = {"Intrinsic": ["Luminal", "Basal", "Her2", "Normal",
                                 "Claudin-low"],
                   "SCMOD2": ["ER-/HER2-", "ER+/HER2- Low Prolif",
                              "ER+/HER2- High Prolif", "HER2+"],
                   "IHC": ["IHC_TNBC", "IHC_ER", "IHC_HER2", "IHC_PR"]}

base_data_dir = "/home/fabio/Downloads/unpast_trans/data"
out_dir = "/home/Desktop/"
exprs_file_t = join(
    base_data_dir,
    "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv")
basename_t = "TCGA"

exprs_file_m = join(
    base_data_dir, "METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv")
basename_m = "METABRIC"

m_subtypes = pd.read_csv(
    join(base_data_dir,
         "METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv"),
    sep="\t", index_col=0)
m_annotation = pd.read_csv(
    join(base_data_dir, "METABRIC_1904.annotation_v6.tsv"),
    sep="\t", index_col=0)

t_subtypes = pd.read_csv(
    join(base_data_dir,
         "TCGA-BRCA_1079_17Kgenes.Xena_" +
         "TCGA_PanCan.subtypes_and_signatures_v6.tsv"),
    sep="\t", index_col=0)
t_annotation = pd.read_csv(
    join(base_data_dir, "TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv"),
    sep="\t", index_col=0)

exprs_t = pd.read_csv(exprs_file_t, sep="\t", index_col=0, nrows=1)
exprs_m = pd.read_csv(exprs_file_m, sep="\t", index_col=0, nrows=1)

known_groups_t, freqs_t = make_ref_groups(t_subtypes, t_annotation, exprs_t)
known_groups_m, freqs_m = make_ref_groups(m_subtypes, m_annotation, exprs_m)


gdc_path = "/home/fabio/Downloads/desmod_run/Bicon/new/bicon/GDC"
mbr_path = "/home/fabio/Downloads/desmod_run/Bicon/new/bicon/Mbr"

for i in range(len([x for x in listdir(gdc_path) if ".tsv" not in x])):
    subt_t = []
    subt_m = []

    for n in range(len([x for x in listdir(
            join(gdc_path, str(i))) if ".tsv" not in x])):

        result_t = prepar_results(gdc_path, join(str(i), str(n)))
        performance_t = calculate_perfromance(result_t, known_groups_t,
                                              freqs_t, set(
                                                  exprs_t.columns.values),
                                              classifications=classifications)
        subt_t.append(performance_t)

        result_m = prepar_results(mbr_path, join(str(i), str(n)))
        performance_m = calculate_perfromance(result_m, known_groups_m,
                                              freqs_m, set(
                                                  exprs_m.columns.values),
                                              classifications=classifications)
        subt_m.append(performance_m)

    pd.DataFrame.from_records(subt_t).sort_index(axis=1).to_csv(
        join(gdc_path, str(i), "BICON_TCGA.tsv"), sep="\t")

    pd.DataFrame.from_records(subt_m).sort_index(axis=1).to_csv(
        join(mbr_path, str(i), "BICON_METABRIC.tsv"), sep="\t")
