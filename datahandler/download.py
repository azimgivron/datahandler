import argparse
import json
import logging
import urllib.request

import requests
from datahandler.data_structure import NetworkConfig, OMIMAPIConfig
from datahandler.const import GENE_GENE_SIMILARITY_URL

"""
Module to download gene–gene similarity networks and
retrieve gene–disease and disease–disease data via the OMIM API.
Provides CLI entry points via a main function.
"""

def download_from_url(config: NetworkConfig) -> None:
    """
    Download data from a URL specified in a NetworkConfig.

    Args:
        config: NetworkConfig with url and filename attributes.

    Raises:
        URLError: If the URL is unreachable.
        IOError: If writing to disk fails.
    """
    urllib.request.urlretrieve(config.url, filename=config.filename)
    logging.info("Downloaded %s", config.filename)


class OMIMClient:
    """
    Client for interacting with the OMIM API.
    """
    def __init__(self, config: OMIMAPIConfig = None):
        self.config = config or OMIMAPIConfig()
        if not self.config.api_key:
            raise ValueError("OMIM_API_KEY environment variable not set.")

    def _get(self, endpoint: str, params: dict) -> dict:
        url = f"{self.config.base_url}/{endpoint}"
        params.update({"apiKey": self.config.api_key, "format": "json"})
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json().get('omim', {})

    def fetch_all_gene_maps(self, batch_size: int = 500) -> list:
        """Page through the geneMap endpoint to retrieve all gene–disease associations.
        
        A geneMap (accessed via the /geneMap endpoint) is OMIM’s structured
        representation of gene–phenotype (disease) associations. Each geneMap
        entry corresponds to one OMIM “gene–phenotype relationship” and contains,
        among other fields:
        
        * mimNumber – the OMIM entry number for the phenotype (disease)
        * entrezGeneId or geneSymbols – identifiers and official symbols for 
            the gene(s) implicated
        * relationshipType – how the gene relates to the phenotype 
            (e.g. “causal,” “susceptibility,” “marker”)
        * phenotypeMapKey – a numeric code indicating the level of certainty
            or type of evidence for the association
        * cytoLocation – the chromosomal location of the gene
        * comments – free-text notes on modifiers, allelic variants, or other nuances
        
        """
        all_maps = []
        start = 0
        while True:
            data = self._get("geneMap", {"start": start, "limit": batch_size})
            chunk = data.get('geneMapList', [])
            if not chunk:
                break
            all_maps.extend(chunk)
            start += batch_size
        return all_maps


def save_json(data, out_file: str) -> None:
    """
    Save Python object as JSON to disk.

    Args:
        data: Serializable Python object.
        out_file: Path to output file.
    """
    with open(out_file, "w") as f:
        json.dump(data, f)
    logging.info("Saved JSON data to %s", out_file)


# Usage wrappers
def download_gene_gene_network() -> None:
    """Download the HumanNet v1 gene–gene similarity network."""
    config = NetworkConfig(
        url=GENE_GENE_SIMILARITY_URL,
        filename="gene-gene-similarity-network.txt",
    )
    download_from_url(config)

def retrieve_all_gene_disease_associations(
    out_file: str = "all_gene_disease.json",
) -> None:
    """Fetch and save all gene–disease associations from OMIM."""
    client = OMIMClient()
    associations = client.fetch_all_gene_maps()
    save_json(associations, out_file)

def main():
    """
    Command-line interface for network downloads and OMIM data retrieval.
    """
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(
        description="Download networks and OMIM data via CLI."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Gene-gene
    parser_gg = subparsers.add_parser(
        "download-gene-gene", help="Download gene–gene similarity network"
    )

    # Gene-disease
    parser_gd = subparsers.add_parser(
        "download-gene-disease", help="Fetch all gene–disease associations"
    )
    parser_gd.add_argument(
        "--output", "-o", default="all_gene_disease.json", help="Output JSON file"
    )

    # Disease entries
    parser_de = subparsers.add_parser(
        "download-disease-entries", help="Fetch disease entries by MIM number"
    )
    parser_de.add_argument(
        "mim_numbers", nargs="+", type=int, help="List of OMIM MIM numbers"
    )
    parser_de.add_argument(
        "--output", "-o", default="omim_disease_entries.json", help="Output JSON file"
    )

    args = parser.parse_args()
    if args.command == "download-gene-gene":
        download_gene_gene_network()
    elif args.command == "download-gene-disease":
        retrieve_all_gene_disease_associations(out_file=args.output)
    elif args.command == "download-disease-entries":
        retrieve_disease_entries(mim_numbers=args.mim_numbers, out_file=args.output)

