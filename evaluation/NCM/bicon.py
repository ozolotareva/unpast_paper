from bicon import (BiCoN, data_preprocessing, results_analysis)

# Network An interaction network should be present as a table with two columns
# that represent two interacting genes. Without a header!

path_expr = "/root/projects/data/real_data/" + \
    "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv"
# "/root/projects/data/real_data/METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv"
path_net = "/root/projects/data/outputs/bicon_network.tsv"

size = range(2000, 5001, 1000)[0]
express_data, network, labels, _ = data_preprocessing(
    path_expr, path_net, log2=False, zscores=False, size=size,
    formats=["tsv", "tsv"])

# minimal solution subnetwork size
L_g_min = 10
# maximal solution subnetwork size
L_g_max = 15

model = BiCoN(express_data, network, L_g_min, L_g_max)
solution, scores = model.run_search()
results = results_analysis(solution, labels)
results.patients1
results.genes1

results = results_analysis(solution, labels, convert=True, origID='symbol')
results.save(output="results/results.csv")
results.show_networks(express_data, network, output="results/network.png")
results.show_clustermap(express_data, network, output="results/clustermap.png")
results.enrichment_analysis(
    library='GO_Biological_Process_2018', output="results")
results.convergence_plot(scores)
