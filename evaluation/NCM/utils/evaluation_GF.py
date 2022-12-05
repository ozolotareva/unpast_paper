import pandas as pd

from os.path import join
from re import sub

from evaluation_desmond import (calculate_perfromance, make_ref_groups)


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


result_t = pd.read_csv(
    "/home/fabio/Downloads/desmod_run/GF/clusters_tcga.tsv", sep="\t")
result_t.index = [sub("\\.", "-", x, ) for x in result_t.index]

result_m = pd.read_csv(
    "/home/fabio/Downloads/desmod_run/GF/clusters_mbr.tsv", sep="\t")
result_m.index = [sub("\\.", "-", x, ) for x in result_m.index]

subt_t = []
subt_m = []

for i, run in enumerate(result_m.columns):

    # k = 5, hence range(1, 6)
    new_result_t = pd.DataFrame([], columns=["samples"], index=range(1, 6))
    new_result_t.loc[:, "samples"] = [
        result_t.index.values[result_t.loc[:, run] == x] for x in range(1, 6)]

    performance_t = calculate_perfromance(new_result_t, known_groups_t,
                                          freqs_t, set(exprs_t.columns.values),
                                          classifications=classifications)
    subt_t.append(performance_t)

    new_result_m = pd.DataFrame([], columns=["samples"], index=range(1, 6))
    new_result_m.loc[:, "samples"] = [
        result_m.index.values[result_m.loc[:, run] == x] for x in range(1, 6)]
    performance_m = calculate_perfromance(new_result_m, known_groups_m,
                                          freqs_m, set(exprs_m.columns.values),
                                          classifications=classifications)
    subt_m.append(performance_m)

pd.DataFrame.from_records(subt_t).to_csv(
    "/home/fabio/Downloads/desmod_run/GF/eval_gdc.tsv", sep="\t")
pd.DataFrame.from_records(subt_m).to_csv(
    "/home/fabio/Downloads/desmod_run/GF/eval_mbr.tsv", sep="\t")
