from logging import (basicConfig, INFO, getLogger)
from re import (compile, sub)
from requests import (get, post, Session)
from requests.adapters import (HTTPAdapter, Retry)
from time import sleep
from tqdm import tqdm
from urllib.parse import (quote_plus, urlparse, parse_qs, urlencode)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def uniprot2gene_old(prot):
    query = quote_plus(
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  ' +
        'virtualSchemaName = "default" formatter = "TSV" header = "1" ' +
        'uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >' +
        '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >' +
        f'<Filter name = "uniprotswissprot" value = "{prot}"/>' +
        '<Attribute name = "uniprotswissprot" />' +
        '<Attribute name = "external_synonym" />' +
        '<Attribute name = "external_gene_name" />' +
        # '<Attribute name = "ensembl_gene_id" />' +
        '<Attribute name = "hgnc_symbol" /></Dataset></Query>')
    convertion_table = pd.read_table(
        'http://www.ensembl.org/biomart/martservice?query=' + query,
        sep="\t").dropna()
    return convertion_table


def get_next_link(headers):
    re_next_link = compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_response, session):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield batch_response.json()
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results):
    for key in ("results", "failedIds"):
        if key in batch_results and batch_results[key]:
            all_results[key] += batch_results[key]
    return all_results


def get_id_mapping_results_search(url, session):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    size = 500
    query["size"] = size
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    request = session.get(parsed.geturl())
    results = request.json()
    for batch in get_batch(request, session):
        results = combine_batches(results, batch)
    return results


def check_id_mapping_results_ready(job_id, session, api_url):
    while True:
        request = session.get(f"{api_url}/idmapping/status/{job_id}")
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                sleep(2)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def uniprot2gene_name(ids):
    api_url = "https://rest.uniprot.org"
    session = Session()
    session.mount("https://", HTTPAdapter(max_retries=Retry(total=5,
                  backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])))

    job_id = post(f"{api_url}/idmapping/run", data={"from": "UniProtKB_AC-ID",
                  "to": "Gene_Name", "ids": ",".join(ids)}).json()["jobId"]

    if check_id_mapping_results_ready(job_id, session, api_url):
        link = session.get(
            f"{api_url}/idmapping/details/{job_id}").json()["redirectURL"]
        results = get_id_mapping_results_search(link, session)

        results_df = pd.DataFrame(results['results'])
        results_df.rename(
            columns={"from": "uniprot", "to": "geneName"}, inplace=True)

        return results_df


basicConfig(filename="/root/projects/data/outputs/ppi_log.txt",
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            filemode='a', datefmt='%H:%M:%S', level=INFO)
logger = getLogger()
logger.info("Retrieving PPI network")

ppi = pd.DataFrame(
    get("https://api.nedrex.net/ppi?taxid=9606&iid_evidence=exp").json())
ppi_subset = ppi[["memberOne", "memberTwo"]].apply(
    lambda x: [sub("uniprot.", "", i) for i in x], axis=0)
ppi_subset.rename(columns={"memberOne": "source",
                           "memberTwo": "target"}, inplace=True)

# check errors self interaction
aaa = ppi_subset.apply(lambda x: x[0] != x[1], axis=1)
# remove self-interactions
ppi_subset = ppi_subset.loc[aaa, :]
logger.info(f"Self-interactions: {ppi.shape[0] - ppi_subset.shape[0]}")

unique_prot = ppi_subset.stack().unique().tolist()
size = 1000
for x in tqdm(range(0, len(unique_prot), size)):
    if x == 0:
        convertion_df = uniprot2gene_name(unique_prot[0+x:size+x])
    else:
        convertion_df = pd.concat(
            [convertion_df, uniprot2gene_name(unique_prot[0+x:size+x])])

ppi_subset = ppi_subset.merge(convertion_df.rename(columns={
    "uniprot": "source",
    "geneName": "source_symbol"}),
    how="left", on="source")

ppi_subset = ppi_subset.merge(convertion_df.rename(columns={
    "uniprot": "target",
    "geneName": "target_symbol"}),
    how="left", on="target")

# vis NaN
ppi_na = ppi_subset.isna()
sns.heatmap(ppi_na, yticklabels=False, cbar=False, cmap='viridis')
plt.ylabel("")
plt.savefig("na_network.png", dpi=600)
plt.close()

ppi_filtered = ppi_subset.dropna()
final_shape = ppi_filtered.shape[0]
shape_difference = ppi_subset.shape[0] - final_shape
fail_per = round((shape_difference/final_shape)*100, 2)
logger.info(
    f"Interactions lost due to mapID failure: {shape_difference}({fail_per}%)")
logger.info(f"N# Interactions: {final_shape}")

ppi_filtered.to_csv("/root/projects/data/outputs/ppi_network.tsv", sep="\t")
