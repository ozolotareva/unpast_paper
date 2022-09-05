import pandas as pd
from fisher import pvalue
from statsmodels.stats.multitest import fdrcorrection


def apply_fdr(df_pval,known_groups):
    df_fdr = {}
    for group in list(known_groups.keys()):
        bh_res, adj_pval = fdrcorrection(df_pval[group].values,alpha=0.05)
        df_fdr[group] =  adj_pval
    df_fdr = pd.DataFrame.from_dict(df_fdr)
    df_fdr.index = df_pval.index
    # df_fdr["associated"] = df_fdr.apply(lambda row: row[row<0.05].index.values,axis=1)
    return df_fdr


def evaluate_overlaps(biclusters, known_groups, N):
    # compute exact Fisher's p-values and Jaccard overlaps
    pvals = {}
    is_enriched = {}
    jaccards = {}
    for i in biclusters.index.values:
        #print(i)
        bic = biclusters.loc[i, :]
        pvals[i] = {}
        is_enriched[i] = {}
        jaccards[i] = {}
        bic_samples = bic["samples"]
        for group in list(known_groups.keys()):
            group_samples = known_groups[group]
            shared = len(bic_samples.intersection(group_samples))
            bic_only = len(bic_samples.difference(group_samples))
            group_only = len(group_samples.difference(bic_samples))
            union = shared + bic_only + group_only
            pval = pvalue(shared, bic_only, group_only, N - union)
            if pval.right_tail < pval.left_tail:
                pvals[i][group] = pval.right_tail
                is_enriched[i][group] = True
                jaccards[i][group] = shared / union
            else:
                # save left-tail p-value and record that this is not enrichment
                pvals[i][group] = pval.left_tail
                is_enriched[i][group] = False
                # compute J for bicluster and group complement, (e.g. not_LumA instead of LumA)
                shared_compement = len(bic_samples) - shared
                union_compement = N - group_only
                jaccards[i][group] = shared_compement / union_compement

    pvals = pd.DataFrame.from_dict(pvals).T
    is_enriched = pd.DataFrame.from_dict(is_enriched).T
    jaccards = pd.DataFrame.from_dict(jaccards).T
    return pvals, is_enriched, jaccards


def find_best_matches(biclusters, known_groups, N, FDR=0.05,  min_SNR=False, min_n_genes=False):
    # for each known group identifies the best matching bicluster - a significant match with max.Jaccard
    # if bicluster is significantly under-represented in a known group,
    # compares bicluster to group complement (e.g. not LumA instead of LumA) and sets is_enriched = False
    # returns all
    results = {}
    if min_SNR:
        biclusters = biclusters.loc[biclusters["avgSNR"] >= min_SNR, :]
    if min_n_genes:
        biclusters = biclusters.loc[biclusters["n_genes"] >= min_n_genes, :]
    if biclusters.shape[0] == 0:
        for group in sorted(list(known_groups.keys())):
            results[group] = {"group_size": group_size, "J": 0, "J_weighted": 0}
        return results

    # calculate overlap p-vals
    df_pval, is_enriched, df_jaccard = evaluate_overlaps(biclusters, known_groups, N)

    # BH-adjust for multiple testing
    df_fdr = apply_fdr(df_pval, known_groups)
    # print(df_pval)
    # print(df_fdr)

    for group in sorted(list(known_groups.keys())):
        significant_matches_j = df_jaccard.loc[df_fdr[group] <= FDR, :]
        group_size = len(known_groups[group])
        if significant_matches_j.shape[0] > 0:
            significant_matches_j = significant_matches_j.sort_values(by=group, ascending=False).head(1)

            bm_id = significant_matches_j.index[0]

            bm = biclusters.loc[bm_id, :]
            j = df_jaccard.loc[bm_id, group]
            results[group] = {"group_size": group_size, "J": j,  # "J_weighted": j*len(known_groups[group])/N,
                              "is_enriched": is_enriched.loc[bm_id, group],
                              "best_match_id": bm_id}
            results[group].update(bm.to_dict())
        else:
            results[group] = {"group_size": group_size, "J": 0}
    results = pd.DataFrame.from_dict(results).T
    results.index.name = "known_group"
    total_bicluster_samples = results["group_size"].sum()
    results["J_weighted"] = results["J"] * results["group_size"] / total_bicluster_samples
    return results


def run_eval(expr_file, ground_truth_file, result_file):
    # expression file
    exprs = pd.read_csv(expr_file, sep="\t", index_col=0, header=0)
    N = len(exprs.columns)
    # read ground truth from file
    ground_truth = pd.read_csv(ground_truth_file, sep ="\t", index_col=0)
    ground_truth["samples"] = ground_truth["samples"].apply(lambda x: set(x.split(" ")))
    if "genes" in ground_truth.columns.values:
        ground_truth["genes"] = ground_truth["genes"].apply(lambda x: set(x.split(" ")))

    # ground_truth

    # prepare a dict with sample groups corresponding to known bicluster
    known_groups = {}
    for group in ground_truth.index.values:
        known_groups[group] = ground_truth.loc[group, "samples"]

    result = pd.read_csv(result_file, sep="\t", index_col=0)
    result['samples'] = result['samples'].apply(eval)

    best_matches = find_best_matches(biclusters=result, known_groups=known_groups, N=N, FDR=0.05)
    return best_matches["J_weighted"].sum()

