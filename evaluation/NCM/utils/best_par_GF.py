import pandas as pd
import numpy as np


def get_best(df):
    score = df["overall_performance_Intrinsic"]
    medias = []

    for i in range(0, len(score), 5):
        medias.append(round(np.average(score[i:i+5]), 3))

    return medias.index(np.max(medias))*5


gdc = pd.read_csv(
    "/home/fabio/Downloads/desmod_run/GF/GrandForest_TCGA.tsv", sep="\t")
mbr = pd.read_csv(
    "/home/fabio/Downloads/desmod_run/GF/GrandForest_METABRIC.tsv", sep="\t")


indice_gdc = get_best(gdc)
gdc.loc[indice_gdc-5:indice_gdc-1,
        :].to_csv("/home/fabio/Downloads/desmod_run/GF/gdc_filtered.tsv",
                  sep="\t", index=False)

indice_mbr = get_best(mbr)
mbr.loc[indice_mbr-5:indice_mbr-1,
        :].to_csv("/home/fabio/Downloads/desmod_run/GF/mbr_filtered.tsv",
                  sep="\t", index=False)
