from os.path import join
from re import sub
import pandas as pd


from eval import find_best_matches, make_known_groups
from eval import find_best_matching_biclusters


def make_ref_groups(subtypes, annotation, exprs):
    # prepared a dict of subtype classifications
    # {"class1":{"subt1":[],"subt2":[]},"class2":{"subtA":[],"subtB":[]}}
    all_samples = set(exprs.columns.values)
    pam50 = make_known_groups(
        subtypes, exprs, target_col="PAM50", verbose=False)
    lum = {}
    lum["Luminal"] = pam50["LumA"].union(pam50["LumB"])
    scmod2 = make_known_groups(
        subtypes, exprs, target_col='SCMOD2', verbose=False)
    claudin = {}
    claudin["Claudin-low"] = set(
        subtypes.loc[subtypes['claudin_low']
                     == 1, :].index.values).intersection(all_samples)

    ihc = {}
    for x in ["IHC_HER2", "IHC_ER", "IHC_PR"]:
        ihc[x] = set(annotation.loc[annotation[x]
                     == "Positive", :].index.values)
    ihc["IHC_TNBC"] = set(
        annotation.loc[annotation["IHC_TNBC"] == 1, :].index.values)

    known_groups = {"PAM50": pam50, "Luminal": lum,
                    "Claudin-low": claudin, "SCMOD2": scmod2, "IHC": ihc}

    freqs = {}
    N = exprs.shape[1]
    for classification in known_groups.keys():
        for group in known_groups[classification].keys():
            n = len(known_groups[classification][group])
            freqs[group] = n/N

    return known_groups, freqs


def calculate_perfromance(results, known_groups, freqs, all_samples,
                          classifications={"Intrinsic": ["Luminal", "Basal",
                                                         "Her2", "Normal",
                                                         "Claudin-low"]}):
    # finds best matches for each subtype, calcuates J per subtype and overall
    # performance
    N = len(all_samples)
    best_matches = []

    for classification in known_groups.keys():
        bm = find_best_matches(
            results, known_groups[classification], all_samples, FDR=0.05,
            verbose=False)
        best_matches.append(bm)

    best_matches = pd.concat(best_matches, axis=0)
    best_matches = best_matches["J"].to_dict()

    for cl_name in classifications.keys():
        overall_performance = 0
        norm_factor = 0
        for group in classifications[cl_name]:
            overall_performance += best_matches[group]*freqs[group]
            norm_factor += freqs[group]
        overall_performance = overall_performance/norm_factor
        best_matches["overall_performance_"+cl_name] = overall_performance
    return best_matches


def compare_gene_clusters(tcga_result, metabric_result, N):
    # N - total number of genes
    # finds best matched TCGA -> METABRIC and METABRIC -> TCGA
    # calculates % of matched clusterst, number of genes in matched cluster,
    # and the average J index for best matches
    bm = find_best_matching_biclusters(tcga_result, metabric_result, N)
    bm = bm.dropna()
    bm2 = find_best_matching_biclusters(metabric_result, tcga_result, N)
    bm2 = bm2.dropna()

    bm = bm.loc[bm["n_shared"] > 1, :].sort_values(
        by="n_shared", ascending=False)
    bm2 = bm2.loc[bm2["n_shared"] > 1, :].sort_values(
        by="n_shared", ascending=False)

    clust_similarity = {}
    # number of biclusters
    clust_similarity["n_1"] = tcga_result.shape[0]
    clust_similarity["n_2"] = metabric_result.shape[0]
    clust_similarity["percent_matched_1"] = bm.shape[0]/tcga_result.shape[0]
    clust_similarity["percent_matched_2"] = bm2.shape[0] / \
        metabric_result.shape[0]
    clust_similarity["n_shared_genes_1"] = bm.loc[:, "n_shared"].sum()
    clust_similarity["n_shared_genes_2"] = bm2.loc[:, "n_shared"].sum()
    clust_similarity["avg_bm_J_1"] = bm.loc[:, "J"].mean()
    clust_similarity["avg_bm_J_2"] = bm2.loc[:, "J"].mean()

    return clust_similarity, bm, bm2


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
