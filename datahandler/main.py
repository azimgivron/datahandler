import argparse
import datetime
import logging
import os
import sys
from pathlib import Path

from datahandler.download import (
    build_all_matrices,
    download_and_extract_mesh,
    download_genemap2,
    download_humannet,
    parse_genemap2,
)


def main() -> None:
    """Orchestrate the full datahandler processing pipeline."""
    # — Centralized Logging Configuration to STDOUT —
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s DataHandler: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    root_logger.addHandler(handler)

    parser = argparse.ArgumentParser(
        description="Run full datahandler pipeline into one output directory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Root directory for all intermediate and final files",
    )
    parser.add_argument(
        "-y",
        "--year",
        type=int,
        default=datetime.date.today().year,
        help="Year of MeSH descriptors to download",
    )
    args = parser.parse_args()

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    api_key = os.getenv("OMIM_API_KEY")
    if not api_key:
        logging.error("OMIM_API_KEY environment variable not set.")
        raise RuntimeError("Missing OMIM_API_KEY environment variable")

    logging.info("Starting pipeline; outputs will go to %s", out_dir)

    hn_file = download_humannet(out_dir)
    gm2_file = download_genemap2(out_dir)
    gd_json = parse_genemap2(gm2_file, out_dir)
    _, mesh_hj = download_and_extract_mesh(args.year, out_dir)
    build_all_matrices(hn_file, gd_json, mesh_hj, out_dir)

    logging.info("Pipeline complete. All outputs under %s", out_dir)


if __name__ == "__main__":
    main()
