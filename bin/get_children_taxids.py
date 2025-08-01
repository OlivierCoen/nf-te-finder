#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import requests
import sys
import argparse
import xml.etree.ElementTree as ET
import logging

from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
    before_sleep_log,
)

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

# Modern NCBI API
NCBI_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_API_HEADERS = {
    "accept": "application/json",
    "content-type": "application/json"
}

# E-UTILITIES OL API
ESEARCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESEARCH_RETMAX = 1000000000 # max retmax that worked

OUTFILE_SUFFIX = ".species_taxids.txt"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument(
        "--family", type=str, required=True, help="Family name"
    )
    return parser.parse_args()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(taxid: str):
    taxons = [taxid]
    data = {
        "taxons": taxons
    }
    response = requests.post(NCBI_API_URL, headers=NCBI_API_HEADERS, json=data)
    response.raise_for_status()
    return response.json()


@retry(
    retry=retry_if_exception_type(requests.exceptions.HTTPError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_esearch_query(query: str, database: str):
    """
    Query NCBI's db with Esearch API
    """
    params = dict(
        db=database,
        term=query,
        retmax=ESEARCH_RETMAX
    )
    response = requests.get(ESEARCH_BASE_URL, params=params)
    response.raise_for_status()
    return response.text


def parse_ids_from_xml(xml_string: str):
    """
    Parse XML string and get text in <Id>...</Id> blocks
    :param xml_string:
    :return: list of Ids
    """
    # parsing XML returned by API
    root = ET.fromstring(xml_string)
    # species taxon IDs are contained in <Id>...</Id> blocks
    return [id_element.text for id_element in root.findall('.//Id')]


def get_family_taxid(family: str):

    result = send_request_to_ncbi_taxonomy(family)

    if len(result['taxonomy_nodes']) > 1:
        raise ValueError(f"Multiple taxids for family {family}")
    node = result['taxonomy_nodes'][0]

    if "taxonomy" not in node:
        logger.info(f"Could not find taxonomy results for family {family}")
        if "errors" in node:
            for error in node["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)

    return node['taxonomy']['tax_id']


def get_children_taxids(taxid: str) -> list[str]:
    """
   Get list of children taxonomy IDs given a family taxonomy ID
   :param taxid:
   :return: list of SRA experiment IDs
   """
    xml_string = send_esearch_query(
        query=f"txid{taxid}[Subtree]",
        #query=f"txid{taxid}[Subtree] AND species[Rank]",
        database="taxonomy"
    )
    return parse_ids_from_xml(xml_string)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    family = args.family

    family_taxid = get_family_taxid(family)
    logger.info(f"Family taxid: {family_taxid}")

    logger.info(f"Getting children taxids for family {family}")
    children_taxids = get_children_taxids(family_taxid)
    logger.info(f"Obtained {len(children_taxids)} children taxids\n")

    # adding family taxid to children taxids
    children_taxids.append(family_taxid)

    outfile = f"{family}{OUTFILE_SUFFIX}"
    with open(outfile, 'w') as fout:
        for taxid in children_taxids:
            fout.write(f"{taxid}\n")

    logger.info("Done")
