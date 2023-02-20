from itertools import product
from joblib import Parallel, delayed
from os.path import join
from re import sub

import pandas as pd
import numpy as np

from evaluation_desmond import (calculate_perfromance, make_ref_groups)


def eval_perf(n):

    result_t = pd.read_csv(join(
        "/home/fabio/Downloads/desmod_run/GF/GDC", n, "clusters_tcga.tsv"),
        sep="\t")
    result_t.index = [sub("\\.", "-", x, ) for x in result_t.index]

    result_m = pd.read_csv(join(
        "/home/fabio/Downloads/desmod_run/GF/Mbr", n, "clusters_mbr.tsv"),
        sep="\t")
    result_m.index = [sub("\\.", "-", x, ) for x in result_m.index]

    subt_t = []
    subt_m = []

    for _, run in enumerate(result_m.columns):

        # k = 5, hence range(1, 6)
        new_result_t = pd.DataFrame([], columns=["samples"], index=range(1, 6))
        new_result_t.loc[:, "samples"] = [result_t.index.values[
            result_t.loc[:, run] == x] for x in range(1, 6)]

        performance_t = calculate_perfromance(new_result_t, known_groups_t,
                                              freqs_t, set(
                                                  exprs_t.columns.values),
                                              classifications=classifications)
        subt_t.append(performance_t)

        new_result_m = pd.DataFrame([], columns=["samples"], index=range(1, 6))
        new_result_m.loc[:, "samples"] = [result_m.index.values[result_m.loc[
            :, run] == x] for x in range(1, 6)]
        performance_m = calculate_perfromance(new_result_m, known_groups_m,
                                              freqs_m, set(
                                                  exprs_m.columns.values),
                                              classifications=classifications)
        subt_m.append(performance_m)

    results_t = pd.DataFrame.from_records(subt_t)
    results_t.to_csv(join(
        "/home/fabio/Downloads/desmod_run/GF/GDC", n, "eval_gdc.tsv"),
        sep="\t")

    results_m = pd.DataFrame.from_records(subt_m)
    results_m.to_csv(join(
        "/home/fabio/Downloads/desmod_run/GF/Mbr", n, "eval_mbr.tsv"),
        sep="\t")

    return (results_t, results_m)


def expansor(df, n_rep):
    final_list = np.ones((n_rep*df.shape[0], 1), dtype="object")
    for n in range(df.shape[0]):
        first = [";".join(["=".join(x) for x in zip(
            df.columns[:-1], df.loc[n, "num.trees":"subgraph"])])]*n_rep
        last = ["=".join(x) for x in product(
            [df.columns[-1]], df.loc[n, "seeds"].split("|"))]
        final_list[n*n_rep:(n*n_rep)+n_rep, 0] = [";".join(i)
                                                  for i in zip(first, last)]
    return final_list


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

results = Parallel(n_jobs=-1)(delayed(eval_perf)(str(n)) for n in range(1, 28))

# NOTE merge evals
for i, (subt_t, subt_m) in enumerate(results):
    if i > 0:
        results_t = pd.concat([results_t, subt_t], ignore_index=True)
        results_m = pd.concat([results_m, subt_m], ignore_index=True)
    else:
        results_t = subt_t
        results_m = subt_m


# NOTE expanding parameters
par_gdc = pd.read_csv(
    "/home/fabio/Downloads/desmod_run/GF/GDC/par_df.tsv",
    sep="\t", dtype="object")

par_mbr = pd.read_csv(
    "/home/fabio/Downloads/desmod_run/GF/Mbr/par_df.tsv",
    sep="\t", dtype="object")

results_t["param"] = expansor(par_gdc, 5)
results_m["param"] = expansor(par_mbr, 5)

results_t.to_csv(
    "/home/fabio/Downloads/desmod_run/GF/GrandForest_TCGA.tsv", sep="\t")
results_m.to_csv(
    "/home/fabio/Downloads/desmod_run/GF/GrandForest_METABRIC.tsv", sep="\t")
