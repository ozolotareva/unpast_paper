import sys
import numpy as np
import pandas as pd
from fisher import pvalue
from statsmodels.stats.multitest import fdrcorrection

from utils.method import zscore
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

def generate_exprs(data_sizes, g_size=5, frac_samples=[0.05, 0.1, 0.25, 0.5], m=2.0, std=1,
                   z=True,
                   outdir="./", outfile_basename="",
                   g_overlap=False, s_overlap=True,
                   seed=42, add_coexpressed=[]):
    n_genes, N = data_sizes
    biclusters = {}
    coexpressed_modules = []

    # generate background model
    np.random.seed(seed)
    exprs = pd.DataFrame(np.random.normal(loc=0, scale=1.0, size=(n_genes, N)))
    exprs.columns = ["s_" + str(x) for x in exprs.columns.values]
    exprs.index = ["g_" + str(x) for x in exprs.index.values]

    # implant bicluster
    N = exprs.shape[1]
    bic_g = []
    bic_s = []
    bg_g = set(exprs.index.values).difference(set(bic_g))
    bg_s = set(exprs.columns.values).difference(set(bic_s))
    bicluster_genes = []
    for s_frac in frac_samples:
        s_size = int(s_frac * N)
        # select random sets of samples and genes from the background
        bic_genes = list(np.random.choice(list(bg_g), size=g_size, replace=False))
        bic_samples = list(np.random.choice(list(bg_s), size=s_size, replace=False))
        bic_g += bic_genes
        bic_s += bic_samples
        # identify samples outside the bicluster
        if not g_overlap:
            bg_g = bg_g.difference(set(bic_g))
        if not s_overlap:
            bg_s = bg_s.difference(set(bic_s))
        # generate bicluster
        biclusters[s_frac] = {"genes": set(bic_genes), "samples": set(bic_samples), "frac": s_frac,
                              "n_genes": len(bic_genes), "n_samples:": len(bic_samples)}
        bic_exprs = np.random.normal(loc=m, scale=std, size=(g_size, s_size))
        # implant biclusters
        exprs.loc[bic_genes, bic_samples] += bic_exprs
        bicluster_genes += bic_genes

    # add modules of co-expressed genes 
    bg_g = set(exprs.index.values).difference(set(bicluster_genes))
    r = 0.5
    for module in add_coexpressed:
        module_genes = list(np.random.choice(list(bg_g), size=module, replace=False))
        n = exprs.loc[module_genes[0], :]
        for i in range(1, module):
            n_i = n * r + np.sqrt(1 - r ** 2) * exprs.loc[module_genes[i], :]
            exprs.loc[module_genes[i], :] = n_i
        print("\tco-exprs. module ", module, "r=",
              (exprs.loc[module_genes, :].T.corr().sum().sum() - module) / (module ** 2 / 2 - module))
        coexpressed_modules.append(module_genes)

    if z:
        # center to 0 and scale std to 1
        exprs = zscore(exprs)
    biclusters = pd.DataFrame.from_dict(biclusters).T
    # biclusters.set_index("frac",inplace = True,drop=True)
    biclusters_ = biclusters.copy()

    if outfile_basename:
        # overlap extension
        if s_overlap:
            if g_overlap:
                overlap_ext = ",overlap=yes"
            else:
                overlap_ext = ",overlap=s"
        elif g_overlap:
            overlap_ext = ",overlap=g"
        else:
            overlap_ext = ",overlap=no"
        # save expressions 
        exprs_file = outdir + "/" + outfile_basename + ".n_genes=" + str(g_size) + ",m=" + str(m) + ",std=" + str(std)
        exprs_file += overlap_ext + ".exprs_z.tsv"
        print("expressions:", exprs_file)
        exprs.to_csv(exprs_file, sep="\t")

        # save ground truth 
        biclusters["n_genes"] = biclusters["genes"].apply(lambda x: len(x))
        biclusters["n_samples"] = biclusters["samples"].apply(lambda x: len(x))
        biclusters["genes"] = biclusters["genes"].apply(lambda x: " ".join((map(str, sorted(x)))))
        biclusters["samples"] = biclusters["samples"].apply(lambda x: " ".join((map(str, sorted(x)))))

        biclusters_file = outdir + "/" + outfile_basename + ".n_genes=" + str(g_size) + ",m=" + str(m) + ",std=" + str(
            std)
        biclusters_file += overlap_ext + ".biclusters.tsv"
        biclusters.to_csv(biclusters_file, sep="\t")
        print("true bilusters:", biclusters_file)
        biclusters.to_csv(biclusters_file, sep="\t")

    return exprs, biclusters_, coexpressed_modules


def make_known_groups(annot, exprs, target_col="genefu_z", verbose=False):
    samples = set(exprs.columns.values).intersection(set(annot.index.values))
    if verbose:
        print("Total samples:", len(samples), file=sys.stdout)
    annot = annot.loc[list(samples), :]
    groups = set(annot.loc[:, target_col].values)

    known_groups = {}
    for group in groups:
        if group == group:
            group_samples = set(annot.loc[annot[target_col] == group, :].index.values)
            group_samples = group_samples.intersection(samples)
            if len(group_samples) > int(len(samples) / 2):
                print("take complement of ", group, file=sys.stderr)
                group_samples = samples.difference(group_samples)
            known_groups[group] = group_samples  # {"set":group_samples,"complement": samples.difference(group_samples)}
            if verbose:
                print(group, round(len(group_samples) / len(samples), 2),
                      len(group_samples), len(samples.difference(group_samples)))
    return known_groups

def make_ref_groups(subtypes, annotation,exprs):
    # prepared a dict of subtype classifications {"class1":{"subt1":[],"subt2":[]},"class2":{"subtA":[],"subtB":[]}}
    all_samples = set(exprs.columns.values)
    pam50 = make_known_groups(subtypes, exprs,target_col = "PAM50",verbose=False)
    lum = {}
    lum["Luminal"] = pam50["LumA"].union(pam50["LumB"])
    scmod2 = make_known_groups(subtypes, exprs,target_col = 'SCMOD2',verbose=False)
    claudin = {} 
    claudin["Claudin-low"] = set(subtypes.loc[subtypes['claudin_low']==1,:].index.values).intersection(all_samples)
    
    ihc = {}
    for x in ["IHC_HER2","IHC_ER","IHC_PR"]:
        ihc[x] = set(annotation.loc[annotation[x]=="Positive",:].index.values)
    ihc["IHC_TNBC"] = set(annotation.loc[annotation["IHC_TNBC"]==1,:].index.values)
    
    known_groups = {"PAM50":pam50,"Luminal":lum,"Claudin-low":claudin,"SCMOD2":scmod2,"IHC":ihc}
    
    freqs = {}
    N =  exprs.shape[1]
    for classification in known_groups.keys():
        for group in known_groups[classification].keys():
            n = len(known_groups[classification][group])
            freqs[group] = n/N
            
    return known_groups, freqs

def compare_gene_clusters(tcga_result,metabric_result, N):
    # N - total number of genes
    # finds best matched TCGA -> METABRIC and METABRIC -> TCGA
    # calculates % of matched clusterst, number of genes in matched cluster, 
    # and the average J index for best matches 
    bm = find_best_matching_biclusters(tcga_result,metabric_result, N)
    bm = bm.dropna()
    bm2 = find_best_matching_biclusters(metabric_result, tcga_result, N)
    bm2 = bm2.dropna()
    
    bm = bm.loc[bm["n_shared"]>1,:].sort_values(by="n_shared",ascending = False)
    bm2 = bm2.loc[bm2["n_shared"]>1,:].sort_values(by="n_shared",ascending = False)
    
    
    clust_similarity = {}
    # number of biclusters 
    clust_similarity["n_1"] = tcga_result.shape[0]
    clust_similarity["n_2"] = metabric_result.shape[0]
    #print("% matched biclusters:",bm.shape[0]/tcga_result.shape[0],bm2.shape[0]/metabric_result.shape[0])
    clust_similarity["percent_matched_1"] = bm.shape[0]/tcga_result.shape[0]
    clust_similarity["percent_matched_2"] = bm2.shape[0]/metabric_result.shape[0]
    #print("n matched genes:",bm.loc[:,"n_shared"].sum(),bm2.loc[:,"n_shared"].sum())
    clust_similarity["n_shared_genes_1"] = bm.loc[:,"n_shared"].sum()
    clust_similarity["n_shared_genes_2"] = bm2.loc[:,"n_shared"].sum()
    #print("avg. J:",bm.loc[:,"J"].mean(),bm2.loc[:,"J"].mean())
    clust_similarity["avg_bm_J_1"] = bm.loc[:,"J"].mean()
    clust_similarity["avg_bm_J_2"] = bm2.loc[:,"J"].mean()
    
    
    return clust_similarity, bm, bm2

def calculate_perfromance(results, known_groups, freqs, all_samples,
                          classifications={"Intrinsic":["Luminal","Basal","Her2","Normal","Claudin-low"]}):
    # finds best matches for each subtype, calcuates J per subtype and overall performance
    N = len(all_samples)
    best_matches = []
    
    for classification in known_groups.keys():
        bm = find_best_matches(results,known_groups[classification],all_samples,FDR=0.05,verbose = False)
        best_matches.append(bm)
            
    best_matches = pd.concat(best_matches, axis=0)
    best_matches = best_matches["J"].to_dict()
    
    for cl_name in classifications.keys():
        overall_performance = 0
        norm_factor = 0
        for group in classifications[cl_name]:
            overall_performance += best_matches[group]*freqs[group]
            norm_factor +=freqs[group]
        overall_performance = overall_performance/norm_factor 
        best_matches["overall_performance_"+cl_name] = overall_performance
    return best_matches

def apply_fdr(df_pval):
    df_fdr = {}
    for group in df_pval.columns.values:
        bh_res, adj_pval = fdrcorrection(df_pval[group].values, alpha=0.05)
        df_fdr[group] = adj_pval
    df_fdr = pd.DataFrame.from_dict(df_fdr)
    df_fdr.index = df_pval.index
    # df_fdr["associated"] = df_fdr.apply(lambda row: row[row<0.05].index.values,axis=1)
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
            bic_members = bic_members.intersection(all_elements)
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
                #if len(bic_members) > len(group_members):
                    # compute J for bicluster and group complement, (e.g. not_LumA instead of LumA)
               #     shared_complement = len(group_members) - shared
               #     union_complement = N - union + group_only
                #else:
                bic_members = all_elements.difference(bic_members)
                shared_complement = len(bic_members.intersection(group_members))
                union_complement = len(bic_members.union(group_members))
                jaccards[group][i] = 0 if union_complement == 0 else shared_complement / union_complement

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
        df_pval, is_enriched, df_jaccard = evaluate_overlaps(biclusters, known_groups, all_elements, dimension=dimension)
    except:
        print("failed to calculate overlap p-values",file=sys.stderr)
        out = evaluate_overlaps(biclusters, known_groups, all_elements, dimension=dimension)

    # BH-adjust for multiple testing
    df_fdr = apply_fdr(df_pval)

    not_matched_biclusters = list(df_fdr.index.values)
    groups = list(df_fdr.columns.values)
    while len(groups) > 0:
        best_group = None
        best_j = -1
        best_bm = None
        best_bm_id = -1
        best_results = None
        for group in groups:
            # choose biclusters with significant overlaps and not overlapping a bigger set
            passed_biclusters = set(df_fdr.loc[df_fdr[group] <= FDR, :].index.values).intersection(
                not_matched_biclusters)
            significant_matches_j = df_jaccard.loc[list(passed_biclusters), group]

            group_size = len(known_groups[group])
            if verbose:
                print("\t", group, "significant matches:", significant_matches_j.shape[0], file=sys.stdout)
            if significant_matches_j.shape[0] > 0:
                significant_matches_j = significant_matches_j.sort_values().tail(1)

                bm_id = significant_matches_j.index[0]

                bm = biclusters.loc[bm_id, :]
                j = df_jaccard.loc[bm_id, group]
                if best_j < j:
                    best_j = j
                    best_group = group
                    best_bm = bm
                    best_bm_id = bm_id
                    best_results = {"group_size": group_size, "J": j,  # "J_weighted": j*len(known_groups[group])/N,
                                    "is_enriched": is_enriched.loc[bm_id, group],
                                    "best_match_id": bm_id}
        if best_group is not None:
            results[best_group] = best_results
            results[best_group].update(best_bm.to_dict())
            # exclude best group
            groups.remove(best_group)
            # exclude best match
            if match_unique:
                not_matched_biclusters = [x for x in not_matched_biclusters if x != best_bm_id]
        else:
            for group in groups:
                # fill group_size clumns for unmatched groups
                group_size = len(known_groups[group])
                results[group] = {"group_size": group_size, "J": 0}
            break
    results = pd.DataFrame.from_dict(results).T
    results.index.name = "known_group"
    total_bicluster_members = results["group_size"].sum()
    results["J_weighted"] = results["J"] * results["group_size"] / total_bicluster_members
    return results


def find_best_matching_biclusters(bics1, bics2, N, by="genes", adj_pval_thr=0.05):
    # takes two biluster dafaframes from read_bic_table
    # by = "genes" or "samples"
    # N - total number of objects (genes or samples)
    # finds best matches of bics1 biclusters among bics2 biclusters
    n_bics1 = bics1.shape[0]
    n_bics2 = bics2.shape[0]
    best_matches = {}  # OrderedDict({})
    for row1 in bics1.iterrows():
        bic1 = row1[1]
        i1 = row1[0]
        set1 = bic1[by]
        s1 = len(set1)

        best_matches[i1] = {}
        bm_J = 0
        bm_o = 0
        bm_adj_pval = 1
        bm_id = -1

        for row2 in bics2.iterrows():
            bic2 = row2[1]
            i2 = row2[0]
            set2 = bic2[by]
            o = len(set1.intersection(set2))
            J = 0
            adj_pval = 1
            if o > 0:
                s2 = len(set2)
                s2_ = s2 - o
                s1_ = s1 - o
                u = s1_ + s2_ + o
                bg = N - u
                if by == "genes":
                    pval = pvalue(o, s1_, s2_, bg).right_tail
                    J = o * 1.0 / u
                elif by == "samples":
                    # two-tailed pvalue is computed for samples
                    # because bicluster and background sets can be flipped 
                    pval = pvalue(o, s1_, s2_, bg)
                    if pval.right_tail <= pval.left_tail:
                        pval = pval.right_tail
                    elif max(s1,
                             s2) > 0.4 * N:  # try flipping the largest bicluster if it is close to 50% of the cohort
                        pval = pval.left_tail
                        if s1 > s2:  # flip s1
                            u = bg + s2
                            o = s2_
                        else:
                            u = bg + s1
                            o = s1_
                    J = o * 1.0 / u

                adj_pval = pval * n_bics2 * n_bics1 / 2
                if adj_pval < adj_pval_thr:

                    if J > bm_J or (J == bm_J and adj_pval < bm_adj_pval):
                        bm_J = J
                        bm_adj_pval = adj_pval
                        bm_id = i2
                        bm_o = o

        best_matches[i1][by] = set1
        best_matches[i1]["n_" + by] = s1
        if bm_id != -1:
            best_matches[i1]["bm_" + by] = bics2.loc[bm_id, by]
            best_matches[i1]["bm_n_" + by] = bics2.loc[bm_id, "n_" + by]
            best_matches[i1]["n_shared"] = bm_o
            best_matches[i1]["J"] = bm_J
            best_matches[i1]["adj_pval"] = bm_adj_pval
            best_matches[i1]["shared_" + by] = bics2.loc[bm_id, by].intersection(set1)
            # best_matches[i1]["bic_n_"+by] = bics2.loc[bm_id,"n_"+by]
            best_matches[i1]["bm_id"] = bm_id
    best_matches = pd.DataFrame.from_dict(best_matches).T
    return best_matches


def bic_survival(surv_anno,samples,event = "OS",surv_time = "",
                 lr = True, verbose = True):
    # check  complete separation
    # if all events are either inside or outside sample group
    if not surv_time:
        surv_time= event+".time"
    surv_data = surv_anno.copy()
    surv_data = surv_data.dropna(axis=0)
    
    # check zero variance columns:
    v =  surv_data.var()
    for col in v.index:
        if v[col] == 0:
            if verbose:
                print(col,"with zero variance excluded",file =sys.stderr)
            surv_data = surv_data.drop(col,axis=1)
    
    surv_data.loc[:,"x"] = 0
    surv_data.loc[list(set(samples).intersection(set(surv_data.index.values))),"x"] = 1

    pval = np.nan
    hr, upper_95CI, lower_95CI =  np.nan,np.nan,np.nan
    results = {}
    
    events = surv_data[event].astype(bool)
    v1 = surv_data.loc[events, 'x'].var()
    v2 = surv_data.loc[~events, 'x'].var()
    
    if v1 ==0 or v2 ==0:
        if verbose:
            in_bic = surv_data.loc[surv_data["x"]==1,:].shape[0]
            in_bg = surv_data.loc[surv_data["x"]==0,:].shape[0]
            print("perfect separation for biclsuter of  %s/%s samples"%(in_bic,in_bg),
                 "variances: {:.2f} {:.2f}".format(v1, v2), file = sys.stderr)
    else:
        try:
            cph = CoxPHFitter()
            res = cph.fit(surv_data, duration_col=surv_time, event_col= event, show_progress=False)
            res_table = res.summary
            res_table  = res_table#.sort_values("p")    
            pval = res_table.loc["x","p"]
            hr = res_table.loc["x","exp(coef)"]
            upper_95CI = res_table.loc["x","exp(coef) upper 95%"]
            lower_95CI = res_table.loc["x","exp(coef) lower 95%"]
        except:
            pass
    
    results = {"p_value":pval,"HR":hr,
                  "upper_95CI":upper_95CI,
                  "lower_95CI":lower_95CI}
    # Log-rank test
    if lr:
        bic = surv_data.loc[surv_data["x"]==1,:]
        bg = surv_data.loc[surv_data["x"]==0,:]
        
        lr_result = logrank_test(bic.loc[:,surv_time], bg.loc[:,surv_time],
                                 event_observed_A=bic.loc[:,event],
                                 event_observed_B=bg.loc[:,event])
        results["LogR_p_value"] = lr_result.p_value
    
    
    return results

def add_survival(biclusters, # dataframes with biclustes
                 sample_data, # sample annotation
                 event= "OS", surv_time = "", # event and time column names
                 covariates=[],
                 min_n_events=5,
                 verbose = True):
    # if too few events, add na columns
    if sample_data[event].sum() < min_n_events:
        df = biclusters.copy()
        for col in [".p_value",".p_value_BH",
                    ".HR",".upper_95CI",".lower_95CI",
                    ".LogR_p_value",".LogR_p_value_BH"]:
            df[event+col] = np.nan
        return df
    if not surv_time:
        surv_time= event+".time"
    surv_results = {}
    for bic in biclusters.iterrows():
        sample_set = bic[1]["samples"]
        surv_data = sample_data.loc[:,covariates +[event,surv_time]]
        
        surv_results[bic[0]] = bic_survival(surv_data,
                                   sample_set,
                                   event = event,
                                   surv_time = surv_time,
                                   verbose = verbose)
        if "pval" in surv_results[bic[0]].keys():
            if np.isnan(surv_results[bic[0]]["pval"]):
                print("failed to fit CPH model for %s ~ bicluster %s"%(event,bic[0]),
                  file = sys.stderr)
        
    surv_results = pd.DataFrame.from_dict(surv_results).T
    surv_results.columns = [event+"."+x for x in surv_results.columns]
    
    pvals = surv_results.loc[~surv_results[event+".p_value"].isna(),event+".p_value"].values
    bh_res, adj_pval = fdrcorrection(pvals, alpha=0.05)
    surv_results.loc[~surv_results[event+".p_value"].isna(),event+".p_value_BH"] = adj_pval
    
    pvals = surv_results.loc[~surv_results[event+".LogR_p_value"].isna(),event+".LogR_p_value"].values
    bh_res, adj_pval = fdrcorrection(pvals, alpha=0.05)
    surv_results.loc[~surv_results[event+".LogR_p_value"].isna(),event+".LogR_p_value_BH"] = adj_pval
    
    return pd.concat([biclusters, surv_results],axis=1)


def test_sample_overlap(row, sample_set, N):
    # usage: 
    # biclusters_df.apply(lambda row: test_sample_overlap(row, sample_set, N),axis=1)
    # N - total number of samples in dataset
    bic_samples = row["samples"]
    o = len(sample_set.intersection(bic_samples))
    bic_only = len(bic_samples)-o
    sample_set_only = len(sample_set)-o
    bg = N-o-bic_only-sample_set_only
    p= pvalue(o,bic_only,sample_set_only,bg ).right_tail
    #if p<0.001:
    #    print(p,(o,bic_only,sample_set_only,bg),row["genes"])
    return pd.Series({"pval":p,"counts":(o,bic_only,sample_set_only,bg)})

def add_sex(biclusters,males = [],females=[]):
    sample_sets = {}
    #if len(males)>0:
    sample_sets["male"] = set(males)
    #if len(females)>0:
    sample_sets["female"] = set(females)
    
    N = len(males)+len(females)
    dfs =[]
    for sex in sample_sets.keys():
        sample_set = sample_sets[sex]
        df = biclusters.apply(lambda row: test_sample_overlap(row, sample_set, N),axis=1)
        df.columns = [sex+"."+x for x in df.columns]
        bh_res, adj_pval = fdrcorrection(df[sex+".pval"].values, alpha=0.05)
        df[sex+".pval_BH"] =  adj_pval
        dfs.append(df)
    dfs = pd.concat(dfs,axis=1)
    dfs["sex.pval_BH"] = dfs.loc[:,["male.pval_BH","female.pval_BH"]].min(axis=1)     
    dfs["sex"] = ""
    try:
        dfs.loc[dfs["male.pval_BH"]<0.05,"sex"] = "male"
    except:
        pass
    try:
        dfs.loc[dfs["female.pval_BH"]<0.05,"sex"] = "female"
    except:
        pass
    return pd.concat([biclusters,dfs],axis=1)