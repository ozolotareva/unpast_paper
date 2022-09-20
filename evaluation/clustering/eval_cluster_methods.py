import pandas as pd
from fisher import pvalue
from statsmodels.stats.multitest import fdrcorrection
import sys


def apply_fdr(df_pval):
    df_fdr = {}
    for group in df_pval.columns.values:
        bh_res, adj_pval = fdrcorrection(df_pval[group].values,alpha=0.05)
        df_fdr[group] =  adj_pval
    df_fdr = pd.DataFrame.from_dict(df_fdr)
    df_fdr.index = df_pval.index
    #df_fdr["associated"] = df_fdr.apply(lambda row: row[row<0.05].index.values,axis=1)
    return df_fdr


def evaluate_overlaps(biclusters, known_groups, all_elements, dimension="samples"):
    # compute exact Fisher's p-values and Jaccard overlaps for samples
    pvals = {}
    is_enriched = {}
    jaccards = {}
    N = len(all_elements)
    # sanity check - biclusters
    for i in biclusters.index.values:
        bic_members = biclusters.loc[i, dimension]
        if not bic_members.intersection(all_elements) == bic_members:
            print("bicluster {} elements {} are not in 'all_elements'".format(i, " ".join(
                bic_members.difference(all_elements))), file=sys.stderr)
            return
            # sanity check and sorting
    group_names = list(known_groups.keys())
    sorted_group_names = [group_names[0]]  # group names ordered by group size
    for group in group_names:
        group_members = known_groups[group]
        if not group_members.intersection(all_elements) == group_members:
            print(group, "elements are not in 'all_elements'", file=sys.stderr)
            return

        if group != group_names[0]:
            for i in range(len(sorted_group_names)):
                if len(group_members) < len(known_groups[sorted_group_names[i]]):
                    sorted_group_names = sorted_group_names[:i] + [group] + sorted_group_names[i:]
                    break
                elif i == len(sorted_group_names) - 1:
                    sorted_group_names = [group] + sorted_group_names
    # print(sorted_group_names)
    for group in sorted_group_names:
        group_members = known_groups[group]
        pvals[group] = {}
        is_enriched[group] = {}
        jaccards[group] = {}
        for i in biclusters.index.values:
            bic = biclusters.loc[i, :]
            bic_members = bic[dimension]
            shared = len(bic_members.intersection(group_members))
            bic_only = len(bic_members.difference(group_members))
            group_only = len(group_members.difference(bic_members))
            union = shared + bic_only + group_only
            pval = pvalue(shared, bic_only, group_only, N - union)
            if pval.right_tail < pval.left_tail:
                pvals[group][i] = pval.right_tail
                is_enriched[group][i] = True
                jaccards[group][i] = shared / union
            else:
                # save left-tail p-value and record that this is not enrichment
                pvals[group][i] = pval.left_tail
                is_enriched[group][i] = False
                # take complement for the biggest group
                if len(bic_members) > len(group_members):

                    # compute J for bicluster and group complement, (e.g. not_LumA instead of LumA)
                    shared_complement = len(group_members) - shared
                    union_complement = N - union + group_only
                else:
                    shared_complement = len(bic_members) - shared
                    union_complement = N - union + bic_only
                jaccards[group][i] = shared_complement / union_complement

                # print(group,jaccards[group])

    pvals = pd.DataFrame.from_dict(pvals).loc[:, sorted_group_names]
    is_enriched = pd.DataFrame.from_dict(is_enriched).loc[:, sorted_group_names]
    jaccards = pd.DataFrame.from_dict(jaccards).loc[:, sorted_group_names]
    return pvals, is_enriched, jaccards


def find_best_matches(biclusters, known_groups, all_elements, FDR=0.05,
                      min_SNR=False, min_n_genes=False,
                      dimension="samples", verbose=False,
                      match_unique=True):
    # for each known group starting from the largest one,
    # identifies the best matching bicluster - a significant match with max.Jaccard
    # matched biclusters are removed from comparizon
    # if a bicluster is significantly under-represented in a known group,
    # compares bicluster to group complement (e.g. not LumA instead of LumA) and sets is_enriched = False
    # returns all
    results = {}
    if min_SNR:
        biclusters = biclusters.loc[biclusters["avgSNR"] >= min_SNR, :]
    if min_n_genes:
        biclusters = biclusters.loc[biclusters["n_genes"] >= min_n_genes, :]
    if biclusters is None or biclusters.shape[0] == 0:
        for group in sorted(list(known_groups.keys())):
            results[group] = {"J": 0, "J_weighted": 0}
        return results

    # calculate overlap p-vals
    try:
        df_pval, is_enriched, df_jaccard = evaluate_overlaps(biclusters, known_groups,
                                                             all_elements, dimension=dimension)
    except:
        print("failed to calculate overlap p-values", file=sys.stderr)
        out = evaluate_overlaps(biclusters, known_groups, all_elements, dimension=dimension)
        return

    # BH-adjust for multiple testing
    df_fdr = apply_fdr(df_pval)

    not_matched_biclusters = df_fdr.index.values
    for group in df_fdr.columns.values:
        # choose biclusters with significant overlaps and not overlapping a bigger set
        passed_biclusters = set(df_fdr.loc[df_fdr[group] <= FDR, :].index.values).intersection(not_matched_biclusters)
        significant_matches_j = df_jaccard.loc[list(passed_biclusters), group]

        group_size = len(known_groups[group])
        if verbose:
            print("\t", group, "significant matches:", significant_matches_j.shape[0], file=sys.stdout)
        if significant_matches_j.shape[0] > 0:
            significant_matches_j = significant_matches_j.sort_values().tail(1)

            bm_id = significant_matches_j.index[0]

            # exclude best match
            if match_unique:
                not_matched_biclusters = [x for x in not_matched_biclusters if x != bm_id]

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
    total_bicluster_members = results["group_size"].sum()
    results["J_weighted"] = results["J"] * results["group_size"] / total_bicluster_members
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

    best_matches = find_best_matches(biclusters=result, known_groups=known_groups, FDR=0.05)
    return best_matches["J_weighted"].sum()

