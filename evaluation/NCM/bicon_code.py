from itertools import product
from os import (makedirs, path)

from bicon import (BiCoN, data_preprocessing, results_analysis)
import pandas as pd


def analysis(ds, express_data, i, labels, network, l_min, l_max):
    makedirs(path.join(base_dir, ds, str(i)), exist_ok=True)

    model = BiCoN(express_data, network, l_min, l_max)
    solution, scores = model.run_search()
    results = results_analysis(solution, labels)

    results.save(path.join(base_dir, ds, str(i), "results.csv"))
    results.show_networks(
        express_data, network,
        output=path.join(base_dir, ds, str(i), "network.png"))
    results.show_clustermap(
        express_data, network,
        output=path.join(base_dir, ds, str(i), "clustermap.png"))
    results.convergence_plot(
        scores,
        path.join(base_dir, ds, str(i), "convergence_plot.png"))

    return True


def main():

    new_size = 0
    for i, (size, l_min, l_max) in enumerate(par):
        if size != new_size:
            express_data, network, labels, _ = data_preprocessing(
                path_gdc, path_net, log2=False, zscores=False,
                size=size, formats=["tsv", "tsv"])
            new_size = size

        analysis("GDC", express_data, i, labels, network, l_min, l_max)

    new_size = 0
    for i, (size, l_min, l_max) in enumerate(par):
        if size != new_size:
            express_data, network, labels, _ = data_preprocessing(
                path_mbr, path_net, log2=False, zscores=False,
                size=size, formats=["tsv", "tsv"])
            new_size = size
        analysis("Mbr", express_data, i, labels, network, l_min, l_max)


if __name__ == "__main__":

    path_net = "/home/fabio/Downloads/unpast_trans/data/bicon_network.tsv"
    base_dir = "/home/fabio/Downloads/desmod_run/Bicon"

    path_gdc = "/home/fabio/Downloads/unpast_trans/data/" + \
        "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv"
    path_mbr = "/home/fabio/Downloads/unpast_trans/data/" + \
        "METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv"

    size = range(2000, 5001, 1000)
    # minimal solution subnetwork size
    L_g_min = range(5, 20, 5)
    # maximal solution subnetwork size
    L_g_max = range(10, 25, 5)
    par = [x for x in product(size, L_g_min, L_g_max) if x[1] < x[2]]
    pd.DataFrame(
        par, columns=["size", "l_min", "l_max"]).to_csv(
            path.join(base_dir, "par_df.tsv"), sep="\t")

    main()
