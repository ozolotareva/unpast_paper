import os, sys

import pandas as pd
from utils.eval import find_best_matches, make_known_groups
from utils.eval import find_best_matching_biclusters


def make_ref_groups(subtypes, annotation, exprs):
    # prepared a dict of subtype classifications {"class1":{"subt1":[],"subt2":[]},"class2":{"subtA":[],"subtB":[]}}
    all_samples = set(exprs.columns.values)
    pam50 = make_known_groups(subtypes, exprs, target_col="PAM50", verbose=False)
    lum = {}
    lum["Luminal"] = pam50["LumA"].union(pam50["LumB"])
    scmod2 = make_known_groups(subtypes, exprs, target_col='SCMOD2', verbose=False)
    claudin = {}
    claudin["Claudin-low"] = set(subtypes.loc[subtypes['claudin_low'] == 1, :].index.values).intersection(all_samples)

    ihc = {}
    for x in ["IHC_HER2", "IHC_ER", "IHC_PR"]:
        ihc[x] = set(annotation.loc[annotation[x] == "Positive", :].index.values)
    ihc["IHC_TNBC"] = set(annotation.loc[annotation["IHC_TNBC"] == 1, :].index.values)

    known_groups = {"PAM50": pam50, "Luminal": lum, "Claudin-low": claudin, "SCMOD2": scmod2, "IHC": ihc}

    freqs = {}
    N = exprs.shape[1]
    for classification in known_groups.keys():
        for group in known_groups[classification].keys():
            n = len(known_groups[classification][group])
            freqs[group] = n / N

    return known_groups, freqs


def calculate_perfromance(results, known_groups, freqs, all_samples,
                          classifications={"Intrinsic": ["Luminal", "Basal", "Her2", "Normal", "Claudin-low"]}):
    # finds best matches for each subtype, calcuates J per subtype and overall performance
    N = len(all_samples)
    best_matches = []

    for classification in known_groups.keys():
        bm = find_best_matches(results, known_groups[classification], all_samples, FDR=0.05, verbose=False)
        best_matches.append(bm)

    best_matches = pd.concat(best_matches, axis=0)
    best_matches = best_matches["J"].to_dict()

    for cl_name in classifications.keys():
        overall_performance = 0
        norm_factor = 0
        for group in classifications[cl_name]:
            overall_performance += best_matches[group] * freqs[group]
            norm_factor += freqs[group]
        overall_performance = overall_performance / norm_factor
        best_matches["overall_performance_" + cl_name] = overall_performance
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
    elif tool_name == 'coalesce':
        samples = []
        genes = []
        gene_lines = []
        sample_lines = []
        with open(result_file, 'r') as fh:
            for line in fh.readlines():
                if line.startswith("Gene_names"):
                    gene_lines.append(line.strip())
                elif line.startswith("Condition_names"):
                    sample_lines.append(line.strip())
        for entry in range(0, len(gene_lines)):
            genes.append(set(gene_lines[entry].split(":\t")[1].strip().split("\t")))
            samples.append(set(sample_lines[entry].split(":\t")[1].strip().split("\t")))
        return pd.DataFrame({'samples': samples, 'genes': genes})
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
            genes.append(
                set() if '[0]' in gene_lines[entry] else set(gene_lines[entry].split(": ")[1].strip().split(" ")))
            samples.append(
                set() if '[0]' in sample_lines[entry] else set(sample_lines[entry].split(": ")[1].strip().split(" ")))

        return pd.DataFrame({'samples': samples, 'genes': genes})
    elif tool_name == 'dataframe':
        try:
            result = pd.read_csv(result_file, delimiter="\t", index_col=0)
            result["samples"] = result["samples"].map(lambda a: a.replace('.', '-')).apply(eval)
            result["n_samples"] = result["samples"].apply(len)
            result["genes"] = result["genes"].apply(eval)
            result["n_genes"] = result["genes"].apply(len)
        except:
            return pd.DataFrame()
        return result


def get_file_pairs(directory):
    tcga_files = dict()
    metabric_files = dict()
    for file in os.listdir(directory):
        if not file.endswith("-biclusters_df.tsv"):
            continue
        params = file.split("-biclusters_df.tsv")[0].split('log2_')[1]
        if file.startswith('TCGA'):
            tcga_files[params] = os.path.join(directory, file)
        elif file.startswith('METABRIC'):
            metabric_files[params] = os.path.join(directory, file)
    file_pairs = []
    for params in tcga_files.keys():
        if params in metabric_files:
            file_pairs.append([params, tcga_files[params], metabric_files[params]])
    return file_pairs


def evaluate(wd, methods=None):
    data_path = "/local/DESMOND2_data/v6/preprocessed_v6/"
    classifications = {"Intrinsic": ["Luminal", "Basal", "Her2", "Normal", "Claudin-low"],
                       "SCMOD2": ["ER-/HER2-", "ER+/HER2- Low Prolif", "ER+/HER2- High Prolif", "HER2+"],
                       "IHC": ["IHC_TNBC", "IHC_ER", "IHC_HER2", "IHC_PR"]}
    exprs_file_t = os.path.join(data_path, "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv")

    exprs_file_m = os.path.join(data_path, "METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv")

    m_subtypes = pd.read_csv(
        os.path.join(data_path, "METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv"),
        sep="\t", index_col=0)
    m_annotation = pd.read_csv(os.path.join(data_path, "METABRIC_1904.annotation_v6.tsv"), sep="\t", index_col=0)

    t_subtypes = pd.read_csv(
        os.path.join(data_path, "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv"),
        sep="\t", index_col=0)
    t_annotation = pd.read_csv(os.path.join(data_path, "TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv"),
                               sep="\t", index_col=0)

    exprs_t = pd.read_csv(exprs_file_t, sep="\t", index_col=0)
    exprs_t[exprs_t > 3] = 3
    exprs_t[exprs_t < -3] = -3

    exprs_m = pd.read_csv(exprs_file_m, sep="\t", index_col=0)
    exprs_m[exprs_m > 3] = 3
    exprs_m[exprs_m < -3] = -3

    known_groups_t, freqs_t = make_ref_groups(t_subtypes, t_annotation, exprs_t)
    known_groups_m, freqs_m = make_ref_groups(m_subtypes, m_annotation, exprs_m)

    tool_name_dict = {'qubic': 'QUBIC', 'qubic2': 'QUBIC2', 'isa2': 'ISA2', 'debi': 'DeBi', 'coalesce': 'COALESCE',
                      'xmotifs': 'xMotifs', 'fabia': 'FABIA'}

    for tool_name in os.listdir(wd):
        if not os.path.isdir(os.path.join(wd, tool_name)):
            continue
        if methods is not None and tool_name not in methods:
            continue
        subt_t = []
        subt_m = []
        clustering_similarities = []
        for files in get_file_pairs(os.path.join(wd, tool_name, "default")):
            run = 0
            if "run" in files[0]:
                run = int(files[0].split("run")[1].split('.')[0])
            params_dict = {"parameters": files[0], "run": run, 'seed': None}
            print(params_dict)
            result_file_t = files[1]
            result_file_m = files[2]

            print(set(exprs_t.columns.values))
            try:
                result_t = read_results("dataframe", result_file_t)
                performance_t = calculate_perfromance(result_t, known_groups_t,
                                                      freqs_t, set(exprs_t.columns.values),
                                                      classifications=classifications)
                print(f'{tool_name} {files[0]} TCGA: {performance_t}')
                performance_t.update(params_dict)
                subt_t.append(performance_t)
                t_failed = False
            except:
                print("TCGA biclustering failed with ", files[0], file=sys.stderr)
                print(files[1])
                t_failed = True
                subt_t.append(params_dict)

            try:
                result_m = read_results("dataframe", result_file_m)
                performance_m = calculate_perfromance(result_m, known_groups_m,
                                                      freqs_m, set(exprs_m.columns.values),
                                                      classifications=classifications)
                print(f'{tool_name} {files[0]} METABRIC: {performance_m}')
                performance_m.update(params_dict)
                subt_m.append(performance_m)
                m_failed = False
            except:
                print("METABRIC biclustering failed with ", files[0], file=sys.stderr)
                print(files[2])
                m_failed = True
                subt_m.append(params_dict)

            if not (t_failed or m_failed):
                N = exprs_m.shape[0]
                try:
                    clust_sim, bm, bm2 = compare_gene_clusters(result_t, result_m, N)
                except:
                    print(f"Clustering comparison failed for {files[0]}: {files[1]}, {files[2]}", file=sys.stderr)
                    clust_sim = {}
            else:
                clust_sim = {}
            clust_sim.update(params_dict)
            clustering_similarities.append(clust_sim)

        results = pd.DataFrame.from_records(clustering_similarities)
        print(tool_name)
        print(results)
        results.to_csv(os.path.join(wd, f"{tool_name_dict[tool_name.lower()]}_similarities.tsv"), sep="\t")
        pd.DataFrame.from_records(subt_t).to_csv(
            os.path.join(wd, f"{tool_name_dict[tool_name.lower()]}_TCGA.tsv"), sep="\t")
        pd.DataFrame.from_records(subt_m).to_csv(
            os.path.join(wd, f"{tool_name_dict[tool_name.lower()]}_METABRIC.tsv"), sep="\t")


if __name__ == '__main__':
    import sys

    evaluate(sys.argv[1])
