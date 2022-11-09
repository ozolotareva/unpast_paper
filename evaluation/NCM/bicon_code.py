from concurrent.futures import ProcessPoolExecutor
from itertools import product
from os import makedirs

from bicon import (BiCoN, data_preprocessing, results_analysis)
import pandas as pd


def analysis(ds, path_expr, i, size, l_min, l_max):
    makedirs(f"/root/projects/data/outputs/Bicon/{ds}/{i}", exist_ok=True)

    express_data, network, labels, _ = data_preprocessing(
        path_expr, path_net, log2=False, zscores=False, size=size,
        formats=["tsv", "tsv"])

    model = BiCoN(express_data, network, l_min, l_max)
    solution, scores = model.run_search()
    results = results_analysis(solution, labels)
    # results.patients1
    # results.genes1

    results.save(
        output=f"/root/projects/data/outputs/Bicon/{ds}/{i}/results.csv")
    results.show_networks(
        express_data, network,
        output=f"/root/projects/data/outputs/Bicon/{ds}/{i}/network.png")
    results.show_clustermap(
        express_data, network,
        output=f"/root/projects/data/outputs/Bicon/{ds}/{i}/clustermap.png")
    results.enrichment_analysis(
        library='GO_Biological_Process_2018',
        output=f"/root/projects/data/outputs/Bicon/{ds}/{i}")
    results.convergence_plot(scores)


def main():

    path_gdc = "/root/projects/data/real_data/" + \
        "METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv"
    path_mbr = "/root/projects/data/real_data/" + \
        "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv"

    with ProcessPoolExecutor() as executor:
        for i, (size, l_min, l_max) in enumerate(par):
            executor.submit(analysis, "GDC", path_gdc, i, size, l_min, l_max)
            executor.submit(analysis, "Mbr", path_mbr, i, size, l_min, l_max)


if __name__ == "__main__":

    # Network An interaction network should be present as a table with two
    # columns that represent two interacting genes. Without a header!
    path_net = "/root/projects/data/outputs/bicon_network.tsv"

    size = range(2000, 5001, 1000)
    # minimal solution subnetwork size
    L_g_min = range(5, 20, 5)
    # maximal solution subnetwork size
    L_g_max = range(10, 25, 5)
    par = [x for x in product(size, L_g_min, L_g_max) if x[1] < x[2]]
    pd.DataFrame(
        par, columns=["size", "l_min", "l_max"]).to_csv(
            "/root/projects/data/outputs/Bicon/par_df.tsv", sep="\t")

    main()
