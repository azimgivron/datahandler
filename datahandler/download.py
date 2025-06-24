"""
CLI to download HumanNet and OMIM flat files, download & parse OMIM genemap2.txt
using genemap2-parser, and construct sparse GxD, GxG, and DxD matrices
(NPZ  labeled CSV).
"""

import argparse
import json
import logging
import os
import subprocess
import urllib.request
import pickle
from pathlib import Path

from datahandler.const import GENE_GENE_SIMILARITY_URL
from datahandler.data_structure import NetworkConfig
import matrix  # assumes matrix.py (with CSV export) is on your PYTHONPATH


def download_from_url(config: NetworkConfig) -> None:
    """Download a file from a URL specified in a NetworkConfig.

    Args:
        config (NetworkConfig): Contains `url` (str) and `filename` (str).

    Raises:
        URLError: If the URL is unreachable.
        IOError: If writing to disk fails.
    """
    urllib.request.urlretrieve(config.url, filename=config.filename)
    logging.info("Downloaded %s", config.filename)


def cmd_download_humannet(args: argparse.Namespace) -> None:
    """Download the HumanNet join file.

    Args:
        args (argparse.Namespace): Command-line args with
            `args.output` specifying the download destination.
    """
    cfg = NetworkConfig(
        url=GENE_GENE_SIMILARITY_URL,
        filename=args.output
    )
    download_from_url(cfg)


def cmd_download_genemap2(args: argparse.Namespace) -> None:
    """Download OMIM genemap2.txt using the OMIM API key.

    Constructs the URL as
    `"https://data.omim.org/downloads/{OMIM_API_KEY}/genemap2.txt"`.

    Args:
        args (argparse.Namespace): Command-line args with
            `args.output` specifying the download destination.

    Raises:
        SystemExit: If the `OMIM_API_KEY` environment variable is not set
            or the download fails.
    """
    api_key = os.getenv("OMIM_API_KEY")
    if not api_key:
        logging.error("OMIM_API_KEY environment variable not set.")
        exit(1)
    url = f"https://data.omim.org/downloads/{api_key}/genemap2.txt"
    try:
        urllib.request.urlretrieve(url, filename=args.output)
        logging.info("Downloaded genemap2.txt to %s", args.output)
    except Exception as e:
        logging.error("Failed to download genemap2.txt: %s", e)
        exit(1)


def cmd_parse_genemap2(args: argparse.Namespace) -> None:
    """Parse genemap2.txt using genemap2-parser and emit JSON.

    This invokes the external `parseGeneMap2` tool to produce `output.pickle`,
    then converts its contents into a simple JSON list of
    `{"mimNumber": int, "geneSymbols": "A, B, ..."}` records.

    Args:
        args (argparse.Namespace): Command-line args with:
            - `args.genemap2_file`: path to `genemap2.txt`
            - `args.output_dir`: directory for the parser’s output.pickle
            - `args.output_json`: path for the resulting JSON file

    Raises:
        subprocess.CalledProcessError: If `parseGeneMap2` exits with an error.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    subprocess.run(
        ["parseGeneMap2", "-i", args.genemap2_file, "-o", args.output_dir],
        check=True
    )

    pickle_path = args.output_dir / "output.pickle"
    with open(pickle_path, "rb") as pf:
        parsed = pickle.load(pf)

    json_records = []
    for rec in parsed:
        mim_val     = rec.get("mim_number")
        syms_val    = rec.get("gene_symbols")
        if not mim_val or not syms_val:
            continue
        try:
            mim_int = int(mim_val)
        except ValueError:
            continue
        # gene_symbols comes back as a comma-separated string
        gene_syms = syms_val if isinstance(syms_val, str) else ", ".join(syms_val)
        json_records.append({
            "mimNumber": mim_int,
            "geneSymbols": gene_syms
        })

    with open(args.output_json, "w", encoding="utf-8") as jf:
        json.dump(json_records, jf, indent=2)
    logging.info("Parsed %d gene–disease records → %s",
                 len(json_records), args.output_json)


def cmd_build_matrices(args: argparse.Namespace) -> None:
    """Build gene–disease, gene–gene, and disease–disease matrices.

    Uses `matrix.construct_gene_disease_matrix`, `matrix.construct_gene_gene_matrix`,
    and `matrix.construct_disease_disease_matrix` to produce NPZ and labeled CSV.

    Args:
        args (argparse.Namespace): Command-line args with:
            - `args.network_file`
            - `args.omim_json`
            - `args.mesh_hierarchy_file`
            - `args.gxd_npz`
            - `args.genes_json`
            - `args.gxg_npz`
            - `args.diseases_json`
            - `args.dxd_npz`
    """
    genes, diseases = matrix.construct_gene_disease_matrix(
        associations_file=args.omim_json,
        out_npz=args.gxd_npz
    )
    with open(args.genes_json, "w", encoding="utf-8") as gh:
        json.dump(genes, gh, indent=2)
    with open(args.diseases_json, "w", encoding="utf-8") as dh:
        json.dump(diseases, dh, indent=2)
    logging.info("Wrote genes → %s and diseases → %s",
                 args.genes_json, args.diseases_json)

    matrix.construct_gene_gene_matrix(
        network_file=args.network_file,
        out_npz=args.gxg_npz,
        genes=genes
    )

    matrix.construct_disease_disease_matrix(
        omim_json=args.omim_json,
        mesh_hierarchy_json=args.mesh_hierarchy_file,
        out_npz=args.dxd_npz,
        diseases=diseases
    )

def main() -> None:
    """Parse CLI arguments and dispatch to appropriate command handler."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(
        description=(
            "Download HumanNet, download & parse OMIM genemap2.txt, "
            "and build sparse matrices (NPZ + labeled CSV)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p1 = sub.add_parser(
        "download-humannet",
        help="Fetch HumanNet join file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p1.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("gene-gene-similarity-network.txt"),
        help="Where to save the HumanNet join file"
    )
    p1.set_defaults(func=cmd_download_humannet)

    p2 = sub.add_parser(
        "download-genemap2",
        help="Download OMIM genemap2.txt using OMIM_API_KEY",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p2.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("genemap2.txt"),
        help="Local path for downloaded genemap2.txt"
    )
    p2.set_defaults(func=cmd_download_genemap2)

    p3 = sub.add_parser(
        "parse-genemap2",
        help="Parse genemap2.txt using genemap2-parser and emit JSON",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p3.add_argument(
        "-i", "--genemap2-file",
        type=Path,
        default=Path("genemap2.txt"),
        help="Path to the local genemap2.txt"
    )
    p3.add_argument(
        "-d", "--output-dir",
        type=Path,
        default=Path("genemap2_parsed"),
        help="Directory where parseGeneMap2 writes output.pickle"
    )
    p3.add_argument(
        "-j", "--output-json",
        type=Path,
        default=Path("all_gene_disease.json"),
        help="Final JSON file of {mimNumber, geneSymbols}"
    )
    p3.set_defaults(func=cmd_parse_genemap2)

    p4 = sub.add_parser(
        "build-matrices",
        help="Construct GxD, GxG, and DxD matrices",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p4.add_argument(
        "-m", "--mesh-hierarchy-file",
        type=Path,
        required=True,
        help="MeSH hierarchy JSON (child→parents map, A & C branches only)"
    )
    p4.add_argument(
        "-n", "--network-file",
        type=Path,
        default=Path("gene-gene-similarity-network.txt"),
        help="HumanNet join file"
    )
    p4.add_argument(
        "-j", "--omim-json",
        type=Path,
        default=Path("all_gene_disease.json"),
        help="Gene–disease JSON (from parse-genemap2)"
    )
    p4.add_argument(
        "--gxd-npz",
        type=Path,
        default=Path("gene_disease.npz"),
        help="Output path for gene–disease CSR matrix (.npz)"
    )
    p4.add_argument(
        "--genes-json",
        type=Path,
        default=Path("genes.json"),
        help="Save ordered gene list as JSON"
    )
    p4.add_argument(
        "--gxg-npz",
        type=Path,
        default=Path("gene_gene.npz"),
        help="Output path for gene–gene CSR matrix (.npz)"
    )
    p4.add_argument(
        "--diseases-json",
        type=Path,
        default=Path("diseases.json"),
        help="Save ordered disease list as JSON"
    )
    p4.add_argument(
        "--dxd-npz",
        type=Path,
        default=Path("disease_disease.npz"),
        help="Output path for disease–disease CSR matrix (.npz)"
    )
    p4.set_defaults(func=cmd_build_matrices)

    args = parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    main()
