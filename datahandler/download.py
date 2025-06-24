#!/usr/bin/env python3
"""
End-to-end pipeline to:

 1. Download HumanNet
 2. Download OMIM genemap2.txt
 3. Parse OMIM via genemap2-parser
 4. Download & extract MeSH (A & C branches) hierarchy
 5. Build G×D, G×G, and D×D sparse matrices (NPZ + labeled CSV)

All intermediate and final files are placed under a single output directory.
"""

import collections
import json
import logging
import os
import pickle
import subprocess
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import List, Tuple

import matrix

from datahandler.const import GENE_GENE_SIMILARITY_URL, MESH_URL, OMIM_URL
from datahandler.data_structure import NetworkConfig


def download_humannet(output_dir: Path) -> Path:
    """Download the HumanNet gene–gene similarity join file.

    Args:
        output_dir (Path): Directory under which to save the file.

    Returns:
        Path: Path to the downloaded HumanNet join file.
    """
    out_path = output_dir / "gene-gene-similarity-network.txt"
    cfg = NetworkConfig(url=GENE_GENE_SIMILARITY_URL, filename=str(out_path))
    urllib.request.urlretrieve(cfg.url, filename=cfg.filename)
    logging.info("Downloaded HumanNet → %s", out_path)
    return out_path


def download_genemap2(output_dir: Path) -> Path:
    """Download the OMIM genemap2.txt flat file using the OMIM API key.

    The environment variable `OMIM_API_KEY` must be set.

    Args:
        output_dir (Path): Directory under which to save the file.

    Returns:
        Path: Path to the downloaded genemap2.txt file.

    Raises:
        RuntimeError: If the `OMIM_API_KEY` environment variable is not set.
    """
    api_key = os.getenv("OMIM_API_KEY")
    out_path = output_dir / "genemap2.txt"
    url = f"{OMIM_URL}{api_key}/genemap2.txt"
    urllib.request.urlretrieve(url, filename=str(out_path))
    logging.info("Downloaded OMIM genemap2.txt → %s", out_path)
    return out_path


def parse_genemap2(genemap2_file: Path, output_dir: Path) -> Path:
    """Invoke the genemap2-parser to convert genemap2.txt into JSON.

    Uses the `parseGeneMap2` CLI tool to produce an intermediate pickle, then
    transforms each record into `{"mimNumber": int, "geneSymbols": str}` format.

    Args:
        genemap2_file (Path): Path to the downloaded genemap2.txt.
        output_dir (Path): Directory under which to write parser outputs and JSON.

    Returns:
        Path: Path to the generated JSON file containing gene–disease associations.

    Raises:
        subprocess.CalledProcessError: If the `parseGeneMap2` tool invocation fails.
    """
    intermediate = output_dir / "genemap2_parsed"
    intermediate.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["parseGeneMap2", "-i", str(genemap2_file), "-o", str(intermediate)], check=True
    )

    pickle_path = intermediate / "output.pickle"
    with open(pickle_path, "rb") as pf:
        parsed = pickle.load(pf)

    records: List[dict] = []
    for rec in parsed:
        mim_val = rec.get("mim_number")
        syms_val = rec.get("gene_symbols")
        if mim_val is None or not syms_val:
            continue
        try:
            mim_int = int(mim_val)
        except (TypeError, ValueError):
            continue
        gene_syms = syms_val if isinstance(syms_val, str) else ", ".join(syms_val)
        records.append({"mimNumber": mim_int, "geneSymbols": gene_syms})

    out_json = output_dir / "all_gene_disease.json"
    with open(out_json, "w", encoding="utf-8") as jf:
        json.dump(records, jf, indent=2)
    logging.info("Parsed %d gene–disease records → %s", len(records), out_json)
    return out_json


def download_and_extract_mesh(year: int, output_dir: Path) -> Tuple[Path, Path]:
    """Download the MeSH descriptor zip, extract the XML, and build A/C hierarchy.

    Args:
        year (int): Year of the MeSH release (e.g. 2025).
        output_dir (Path): Directory under which to save and extract files.

    Returns:
        Tuple[Path, Path]:
            - Path to the extracted MeSH XML file.
            - Path to the generated child→parents hierarchy JSON.
    """
    zip_path = output_dir / f"desc{year}.zip"
    xml_path = output_dir / f"desc{year}.xml"
    hierarchy_json = output_dir / "mesh_hierarchy.json"

    zip_url = f"{MESH_URL}desc{year}.zip"
    urllib.request.urlretrieve(zip_url, filename=str(zip_path))
    logging.info("Downloaded MeSH ZIP → %s", zip_path)

    with zipfile.ZipFile(zip_path, "r") as zf:
        xml_files = [n for n in zf.namelist() if n.endswith(".xml")]
        if not xml_files:
            raise RuntimeError(f"No XML found in {zip_path}")
        zf.extract(xml_files[0], path=output_dir)
        extracted_xml = output_dir / xml_files[0]
        extracted_xml.rename(xml_path)
    logging.info("Extracted MeSH XML → %s", xml_path)

    build_mesh_hierarchy_ac(xml_path, hierarchy_json)
    return xml_path, hierarchy_json


def build_mesh_hierarchy_ac(xml_path: Path, json_path: Path) -> None:
    """Filter MeSH XML to the Anatomy (A) and Disease (C) branches.

    Parses the MeSH descriptor XML and writes a JSON mapping each term UI
    to a list of its immediate parents, restricted to the A and C branches.

    Args:
        xml_path (Path): Path to the MeSH descriptor XML file.
        json_path (Path): Path where the filtered hierarchy JSON will be saved.
    """
    root = ET.parse(xml_path).getroot()
    tree_to_ui: dict[str, str] = {}
    for rec in root.findall("DescriptorRecord"):
        ui = rec.findtext("DescriptorUI")
        for tn in rec.findall("TreeNumberList/TreeNumber"):
            t = tn.text or ""
            if t and t[0] in ("A", "C"):
                tree_to_ui[t] = ui

    children: dict[str, List[str]] = collections.defaultdict(list)
    for rec in root.findall("DescriptorRecord"):
        ui = rec.findtext("DescriptorUI")
        for tn in rec.findall("TreeNumberList/TreeNumber"):
            t = tn.text or ""
            if not t or t[0] not in ("A", "C") or "." not in t:
                continue
            parent_t = t.rsplit(".", 1)[0]
            if parent_t[0] not in ("A", "C"):
                continue
            p_ui = tree_to_ui.get(parent_t)
            if p_ui and p_ui != ui:
                children[ui].append(p_ui)

    with open(json_path, "w", encoding="utf-8") as fp:
        json.dump(children, fp, indent=2)
    logging.info("Extracted %d MeSH terms → %s", len(children), json_path)


def build_all_matrices(
    human_net_file: Path, gd_json: Path, mesh_hierarchy_json: Path, output_dir: Path
) -> None:
    """Construct gene–disease, gene–gene, and disease–disease matrices.

    Uses the routines in `matrix.py` to build and save NPZ + labeled CSV.

    Args:
        human_net_file (Path): Path to the HumanNet join file.
        gd_json (Path): Path to the gene–disease JSON file.
        mesh_hierarchy_json (Path): Path to the MeSH A/C hierarchy JSON.
        output_dir (Path): Base directory for saving matrices and index lists.
    """
    # Gene–Disease matrix
    genes, diseases = matrix.construct_gene_disease_matrix(
        associations_file=str(gd_json), out_npz=str(output_dir / "gene_disease.npz")
    )
    with open(output_dir / "genes.json", "w") as f:
        json.dump(genes, f, indent=2)
    with open(output_dir / "diseases.json", "w") as f:
        json.dump(diseases, f, indent=2)

    # Gene–Gene matrix
    matrix.construct_gene_gene_matrix(
        network_file=str(human_net_file),
        out_npz=str(output_dir / "gene_gene.npz"),
        genes=genes,
    )

    # Disease–Disease matrix
    matrix.construct_disease_disease_matrix(
        omim_json=str(gd_json),
        mesh_hierarchy_json=str(mesh_hierarchy_json),
        out_npz=str(output_dir / "disease_disease.npz"),
        diseases=diseases,
    )
