from logging import (basicConfig, INFO, getLogger)

import pandas as pd

basicConfig(filename="/root/projects/data/outputs/biogrid_log.txt",
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            filemode='a', datefmt='%H:%M:%S', level=INFO)
logger = getLogger()
logger.info("Filtering BIOGRID-MV-Physical PPI network")

# wget --retry-connrefused --waitretry=1 --read-timeout=20 \
#       --timeout=15 -t 0 -cv https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.214/BIOGRID-MV-Physical-4.4.214.tab3.zip
ppi = pd.read_csv(
    "/root/projects/data/outputs/BIOGRID-MV-Physical-4.4.214.tab3.txt",
    sep="\t", low_memory=False)

# remove non human inter
human_A = ppi["Organism Name Interactor A"] == "Homo sapiens"
human_B = ppi["Organism Name Interactor B"] == "Homo sapiens"
ppi_subset = ppi.loc[human_A & human_B, :]
logger.info(
    f"Non-human interactions removed: {ppi.shape[0] - ppi_subset.shape[0]}")
logger.info(f"N inters remaining: {ppi_subset.shape[0]}")

# remove sef-inter
ppi_filtered = ppi_subset.loc[ppi_subset["Official Symbol Interactor A"]
                              != ppi_subset["Official Symbol Interactor B"], :]
logger.info(
    f"Self-inter removed: {ppi_subset.shape[0] - ppi_filtered.shape[0]}")

# export network
ppi_filtered = ppi_filtered[[
    "Official Symbol Interactor A", "Official Symbol Interactor B"]]
ppi_filtered.rename(columns={"Official Symbol Interactor A": "source",
                             "Official Symbol Interactor B": "target"},
                    inplace=True)
ppi_filtered.to_csv(
    "/root/projects/data/outputs/BIOGRID-MV-Physical-4.4.214_filtered.tsv",
    sep="\t")

ppi_filtered.to_csv(
    "/root/projects/data/outputs/bicon_network.tsv",
    sep="\t", header=False, index=False)

logger.info(f"Final shape: {ppi_filtered.shape[0]}")
