from os.path import join
from os.path import exists

import pandas as pd

from method import read_bic_table
from eval import (find_best_matches, make_known_groups)


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


def adapt_bicluster(exprs, fname):

    biclusters = pd.read_csv(fname, sep="\t", index_col=0, comment="#")

    biclusters["genes_up"] = pd.NA
    biclusters["genes_down"] = pd.NA
    for i, direction in enumerate(biclusters.direction):
        if direction == "UP":
            biclusters.loc[i, "genes_up"] = biclusters.loc[i, "genes"]
        elif direction == "DOWN":
            biclusters.loc[i, "genes_down"] = biclusters.loc[i, "genes"]

    genes_list = list(exprs.index)
    patients_list = list(exprs.columns)
    biclusters["gene_indexes"] = biclusters.genes.apply(
        lambda x: " ".join([str(genes_list.index(w)) for w in x.split(" ")]))
    biclusters["sample_indexes"] = biclusters.samples.apply(
        lambda x: " ".join(
            [str(patients_list.index(w)) for w in x.split(" ")]))

    biclusters.to_csv(fname[:-4] + "_edited.tsv", sep="\t")

    return True


if __name__ == "__main__":
    classifications = {"Intrinsic": ["Luminal", "Basal", "Her2", "Normal",
                                     "Claudin-low"],
                       "SCMOD2": ["ER-/HER2-", "ER+/HER2- Low Prolif",
                                  "ER+/HER2- High Prolif", "HER2+"],
                       "IHC": ["IHC_TNBC", "IHC_ER", "IHC_HER2", "IHC_PR"]}

    base_data_dir = "/home/fabio/Downloads/unpast_trans/data"
    out_dir = "/home/Desktop/"

    # NOTE GDC
    exprs_file_t = join(
        base_data_dir,
        "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv")

    t_subtypes = pd.read_csv(
        join(base_data_dir,
             "TCGA-BRCA_1079_17Kgenes.Xena_" +
             "TCGA_PanCan.subtypes_and_signatures_v6.tsv"),
        sep="\t", index_col=0)
    t_annotation = pd.read_csv(
        join(base_data_dir,
             "TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv"),
        sep="\t", index_col=0)

    exprs_t = pd.read_csv(exprs_file_t, sep="\t", index_col=0)
    known_groups_t, freqs_t = make_ref_groups(
        t_subtypes, t_annotation, exprs_t)

    # NOTE MBR
    exprs_file_m = join(
        base_data_dir, "METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv")

    m_subtypes = pd.read_csv(
        join(base_data_dir,
             "METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv"),
        sep="\t", index_col=0)
    m_annotation = pd.read_csv(
        join(base_data_dir, "METABRIC_1904.annotation_v6.tsv"),
        sep="\t", index_col=0)

    exprs_m = pd.read_csv(exprs_file_m, sep="\t", index_col=0)
    known_groups_m, freqs_m = make_ref_groups(
        m_subtypes, m_annotation, exprs_m)

    base_path = "/home/fabio/Downloads/desmod_run/D1"

    for folder, alpha in zip(["0_5", "2_5"], ["0.5", "2.5"]):
        subt_t = []
        subt_m = []

        for n in range(5):

            # NOTE GDC
            fname = f"/home/fabio/Downloads/desmod_run/D1/{folder}/GDC/{n}" + \
                    f"/GDC.alpha={alpha},beta_K=1.0,p_val=0.01," + \
                    "q=0.5.biclusters.permutations.tsv"

            if not exists(fname[:-4] + "_edited.tsv"):
                adapt_bicluster(exprs_t, fname)

            result_t = read_bic_table(fname[:-4] + "_edited.tsv")

            # NOTE MBR
            fname = f"/home/fabio/Downloads/desmod_run/D1/{folder}/Mbr/{n}" + \
                    f"/Mbr.alpha={alpha},beta_K=1.0,p_val=0.01," + \
                    "q=0.5.biclusters.permutations.tsv"

            if not exists(fname[:-4] + "_edited.tsv"):
                adapt_bicluster(exprs_m, fname)

            result_m = read_bic_table(fname[:-4] + "_edited.tsv")

            # NOTE calculations
            performance_t = calculate_perfromance(
                result_t, known_groups_t,
                freqs_t, set(
                    exprs_t.columns.values),
                classifications=classifications)
            subt_t.append(performance_t)

            performance_m = calculate_perfromance(
                result_m, known_groups_m,
                freqs_m, set(
                    exprs_m.columns.values),
                classifications=classifications)
            subt_m.append(performance_m)

        pd.DataFrame.from_records(subt_t).sort_index(axis=1).to_csv(
            f"/home/fabio/Downloads/desmod_run/D1/{folder}/GDC/eval_gdc.tsv",
            sep="\t")

        pd.DataFrame.from_records(subt_m).sort_index(axis=1).to_csv(
            f"/home/fabio/Downloads/desmod_run/D1/{folder}/Mbr/eval_mbr.tsv",
            sep="\t")
