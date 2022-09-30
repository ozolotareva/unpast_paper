import os

import pandas as pd
from utils.eval import find_best_matches, make_known_groups


def match_known_subtypes(results, subtypes, annotation, exprs):
    all_samples = set(exprs.columns.values)
    pam50 = make_known_groups(subtypes, exprs, target_col="PAM50", verbose=False)
    lum = {}
    lum["Luminal"] = pam50["LumA"].union(pam50["LumB"])
    scmod2 = make_known_groups(subtypes, exprs, target_col='SCMOD2', verbose=False)
    claudin = {}
    claudin["Claudin-low"] = set(subtypes.loc[subtypes['claudin_low'] == 1, :].index.values).intersection(all_samples)

    ihc = {}
    for x in ["IHC_HER2", "IHC_ER", "IHC_PR", "IHC_TNBC"]:
        ihc[x] = set(annotation.loc[annotation[x] == "Positive", :].index.values)

    known_groups = [pam50, lum, claudin, scmod2, ihc]
    best_matches = []
    for group in known_groups:
        bm = find_best_matches(results, group, all_samples, FDR=0.05, verbose=False)
        best_matches.append(bm)
    best_matches = pd.concat(best_matches, axis=0)
    return best_matches


from utils.eval import find_best_matching_biclusters


def compare_gene_clusters(tcga_result, metabric_result, N):
    # N - total number of genes
    # finds best matched TCGA -> METABRIC and METABRIC -> TCGA
    # calculates % of matched clusterst, number of genes in matched cluster,
    # and the average J index for best matches
    bm = find_best_matching_biclusters(tcga_result, metabric_result, N)
    bm = bm.dropna()
    bm2 = find_best_matching_biclusters(metabric_result, tcga_result, N)
    bm2 = bm2.dropna()

    bm = bm.loc[bm["n_shared"] > 1, :].sort_values(by="n_shared", ascending=False)
    bm2 = bm2.loc[bm2["n_shared"] > 1, :].sort_values(by="n_shared", ascending=False)

    clust_similarity = {}
    # number of biclusters
    clust_similarity["n_1"] = tcga_result.shape[0]
    clust_similarity["n_2"] = metabric_result.shape[0]
    # print("% matched biclusters:",bm.shape[0]/tcga_result.shape[0],bm2.shape[0]/metabric_result.shape[0])
    clust_similarity["percent_matched_1"] = bm.shape[0] / tcga_result.shape[0]
    clust_similarity["percent_matched_2"] = bm2.shape[0] / metabric_result.shape[0]
    # print("n matched genes:",bm.loc[:,"n_shared"].sum(),bm2.loc[:,"n_shared"].sum())
    clust_similarity["n_shared_genes_1"] = bm.loc[:, "n_shared"].sum()
    clust_similarity["n_shared_genes_2"] = bm2.loc[:, "n_shared"].sum()
    # print("avg. J:",bm.loc[:,"J"].mean(),bm2.loc[:,"J"].mean())
    clust_similarity["avg_bm_J_1"] = bm.loc[:, "J"].mean()
    clust_similarity["avg_bm_J_2"] = bm2.loc[:, "J"].mean()

    return clust_similarity, bm, bm2


def read_results(tool_name, result_file):
    if tool_name in ['isa2', 'fabia', 'qubic']:
        result = pd.read_csv(result_file, names=['samples', 'genes'], sep="\t")
        result["samples"] = result["samples"].apply(lambda x: set(x.strip().split(" ")))
        result["genes"] = result["genes"].apply(lambda x: set(x.strip().split(" ")))
        return result
    elif tool_name == 'debi':
        samples = []
        genes = []
        line_nr = -1
        with open(result_file, 'r') as fh:
            for line in fh.readlines():
                line_nr = (line_nr + 1) % 3
                if line_nr == 1:
                    genes.append(set(line.strip().split(" ")))
                elif line_nr == 2:
                    samples.append(set(line.strip().split(" ")))
        return pd.DataFrame({'samples': samples, 'genes': genes})
    elif tool_name == 'qubic2':
        samples = []
        genes = []
        gene_lines = []
        sample_lines = []
        with open(result_file, 'r') as fh:
            for line in fh.readlines():
                if line.startswith(" Genes"):
                    gene_lines.append(line.strip())
                elif line.startswith(" Conds"):
                    sample_lines.append(line.strip())
        for entry in range(0, len(gene_lines)):
            genes.append(set() if '[0]' in gene_lines[entry] else set(gene_lines[entry].split(": ")[1].strip().split(" ")))
            samples.append(
                set() if '[0]' in sample_lines[entry] else set(sample_lines[entry].split(": ")[1].strip().split(" ")))

        return pd.DataFrame({'samples': samples, 'genes': genes})
    elif tool_name == 'dataframe':
        result = pd.read_csv(result_file, delimiter="\t", index_col=0)
        result["samples"] = result["samples"].apply(eval)
        result["genes"] = result["genes"].apply(eval)
        return result


def run_eval(tool_name, expr_file, subtype_file, annot_file,  result_file):
    # expression file

    exprs = pd.read_csv(expr_file, sep="\t", index_col=0)
    exprs[exprs > 3] = 3
    exprs[exprs < -3] = -3
    samples = list(exprs.columns)

    # read subtypes and annotations truth from file
    subtypes = pd.read_csv(subtype_file, sep="\t",index_col=0)
    annotation = pd.read_csv(annot_file, sep="\t", index_col=0)


    result = read_results(tool_name, result_file)
    # print(result)
    # result.to_csv("/home/andim/Downloads/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv.qubic2_df.tsv", sep="\t")
    # print(samples)
    best_matches = match_known_subtypes(result, subtypes, annotation, exprs)
    best_matches = best_matches["J"].to_dict()
    # t_best_matches.update(params_dict)
    # t_best_matches["time"] = time_t
    # best_matches = find_best_matches(result, known_groups, samples, FDR=0.05)
    # try:
    return (best_matches, result)
    # except:
    #     return (0.0,result)


#
#
# name = 'dataframe'
# expr_file = '/home/andim/Downloads/desmond2-eval/C.n_genes=500,m=4,std=1,overlap=yes.exprs_z.tsv'
# truth = '/home/andim/Downloads/desmond2-eval/C.n_genes=500,m=4,std=1,overlap=yes.biclusters.tsv'
# result = '/home/andim/Downloads/desmond2-eval/C.n_genes=500,m=4,std=1,overlap=yes_r=1-q=0.1-c=0.51-f=1-P=F-C=T-type=area-biclusters_df.tsv'
# #
# (scores, result) = run_eval(name, expr_file, truth, result)
# print(scores["J_weighted"].sum())
