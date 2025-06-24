"""
Module for constructing coherent gene–disease (GxD), gene–gene (GxG), and
disease–disease (DxD) matrices from OMIM, HumanNet, and MeSH hierarchy data,
and exporting both sparse NPZ and labeled CSV representations.
"""

import json
import logging
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, save_npz
from sklearn.metrics.pairwise import cosine_similarity


# --- Data Loaders ---
def load_gene_disease_associations(associations_file: str):
    """
    Load OMIM gene–disease associations JSON.

    Returns:
      records: list of (gene, disease) tuples
      genes: sorted list of unique genes
      diseases: sorted list of unique diseases (MIM IDs)
    """
    with open(associations_file, "r", encoding="utf-8") as f:
        assoc = json.load(f)

    genes_set, diseases_set, records = set(), set(), []
    for entry in assoc:
        d = str(entry.get("mimNumber"))
        diseases_set.add(d)
        for g in entry.get("geneSymbols", "").split(", "):
            if g:
                genes_set.add(g)
                records.append((g, d))

    genes = sorted(genes_set)
    diseases = sorted(diseases_set)
    return records, genes, diseases


# --- Gene–Disease Matrix (GxD) ---
def construct_gene_disease_matrix(associations_file: str, out_npz: str):
    """
    Build and save a binary gene–disease sparse matrix (GxD), and export CSV
    with gene names as rows and disease MIMs as columns.

    Args:
      associations_file: JSON file of OMIM gene–disease records.
      out_npz:   path to write the sparse CSR matrix (.npz).

    Returns:
      genes, diseases: lists used for rows and columns.
    """
    records, genes, diseases = load_gene_disease_associations(associations_file)
    idx_g = {g: i for i, g in enumerate(genes)}
    idx_d = {d: i for i, d in enumerate(diseases)}

    rows, cols, data = [], [], []
    seen = set()
    for g, d in records:
        i, j = idx_g[g], idx_d[d]
        if (i, j) not in seen:
            rows.append(i)
            cols.append(j)
            data.append(1)
            seen.add((i, j))

    mat = coo_matrix(
        (data, (rows, cols)), shape=(len(genes), len(diseases)), dtype=int
    ).tocsr()
    save_npz(out_npz, mat)
    logging.info(
        "G×D NPZ saved: %dx%d nnz=%d → %s", mat.shape[0], mat.shape[1], mat.nnz, out_npz
    )

    # export labeled CSV
    df = pd.DataFrame.sparse.from_spmatrix(mat, index=genes, columns=diseases)
    csv_path = out_npz.replace(".npz", ".csv")
    df.to_csv(csv_path)
    logging.info("G×D CSV saved with labels → %s", csv_path)

    return genes, diseases


# --- Gene–Gene Matrix (GxG) ---
def construct_gene_gene_matrix(network_file: str, out_npz: str, genes):
    """
    Build and save a gene–gene similarity sparse matrix (GxG), and export CSV
    with genes as both row and column labels.

    Args:
      network_file: edge-list file, each line 'gene1 gene2 weight'
      out_npz:      path to write the sparse CSR matrix (.npz)
      genes:        list of genes to include (defines ordering)

    Returns:
      genes: passed through for consistency
    """
    idx_g = {g: i for i, g in enumerate(genes)}
    G = len(genes)

    rows, cols, data = [], [], []
    with open(network_file, "r", encoding="utf-8") as f:
        for raw_line in f:
            parts = raw_line.strip().split()
            if len(parts) < 3:
                continue
            g1, g2, w = (
                parts[0],
                parts[1],
                float(parts[2]) if parts[2] != "NA" else np.nan,
            )
            if g1 in idx_g and g2 in idx_g:
                i, j = idx_g[g1], idx_g[g2]
                rows.append(i)
                cols.append(j)
                data.append(w)

    # add self-similarity =1
    for i in range(G):
        rows.append(i)
        cols.append(i)
        data.append(1.0)

    mat = coo_matrix((data, (rows, cols)), shape=(G, G), dtype=np.float32).tocsr()
    save_npz(out_npz, mat)
    logging.info("G×G NPZ saved: %dx%d nnz=%d → %s", G, G, mat.nnz, out_npz)

    # export labeled CSV
    df = pd.DataFrame.sparse.from_spmatrix(mat, index=genes, columns=genes)
    csv_path = out_npz.replace(".npz", ".csv")
    df.to_csv(csv_path)
    logging.info("G×G CSV saved with labels → %s", csv_path)

    return genes


# --- Disease–Disease Matrix (DxD) via MimMiner pipeline ---
def load_omim_entries(omim_json: str) -> dict:
    """
    Load OMIM entries JSON mapping MIM numbers to dicts with 'TX' and 'CS' fields.
    """
    with open(omim_json, "r", encoding="utf-8") as f:
        return json.load(f)


def load_mesh_hierarchy(mesh_hierarchy_json: str) -> dict:
    """
    Load MeSH hierarchy: child term → list of parent terms.
    """
    with open(mesh_hierarchy_json, "r", encoding="utf-8") as f:
        return json.load(f)


def build_children_map(hierarchy: dict) -> dict:
    """
    Invert MeSH hierarchy: parent term → list of direct hyponyms.
    """
    children = defaultdict(list)
    for child, parents in hierarchy.items():
        for p in parents:
            children[p].append(child)
    return children


def extract_mesh_terms(text: str, vocab: set) -> list:
    """
    Extract MeSH concepts by simple string matching (case-insensitive).
    """
    up = text.upper()
    return [term for term in vocab if term.upper() in up]


def construct_disease_disease_matrix(
    omim_json: str,
    mesh_hierarchy_json: str,
    out_npz: str,
    diseases,
) -> list:
    """
    Build and save a disease–disease cosine similarity sparse matrix (DxD), and
    export CSV with disease MIMs as both row and column labels.

    Args:
      omim_json:            OMIM entries JSON with 'TX' & 'CS' text
      mesh_hierarchy_json:  MeSH hierarchy JSON (child→parents)
      out_npz:              path to write the sparse CSR matrix (.npz)
      diseases:             list of MIM IDs defining ordering

    Returns:
      diseases: passed through for consistency
    """
    records = load_omim_entries(omim_json)
    hierarchy = load_mesh_hierarchy(mesh_hierarchy_json)
    children = build_children_map(hierarchy)
    vocab = set(hierarchy.keys())
    D = len(diseases)

    # Steps 1-2: raw and hierarchy-refined counts
    term_doc_freq = defaultdict(int)
    refined_counts = {}
    for mim in diseases:
        entry = records.get(str(mim), {})
        text = (entry.get("TX", "") + " " + entry.get("CS", "")).upper()
        raw = defaultdict(int)
        for t in extract_mesh_terms(text, vocab):
            raw[t] += 1

        ref = raw.copy()
        for parent, hypos in children.items():
            nh = len(hypos)
            if nh:
                s = sum(raw.get(h, 0) for h in hypos)
                if s:
                    ref[parent] += s / nh

        refined_counts[mim] = ref
        for term, cnt in ref.items():
            if cnt > 0:
                term_doc_freq[term] += 1

    # Step 3: IDF
    idf = {t: np.log(D / df) for t, df in term_doc_freq.items()}
    terms = sorted(idf)

    # Step 4: construct TF×IDF matrix
    rows, cols, data = [], [], []
    for i, mim in enumerate(diseases):
        vec = np.array(
            [refined_counts[mim].get(t, 0) * idf[t] for t in terms], dtype=float
        )
        norm = np.linalg.norm(vec)
        if norm:
            vec /= norm
        nz = np.nonzero(vec)[0]
        rows.extend([i] * len(nz))
        cols.extend(nz.tolist())
        data.extend(vec[nz].tolist())

    X = coo_matrix((data, (rows, cols)), shape=(D, len(terms)), dtype=np.float32)

    # Step 5: cosine similarity → full DxD
    S = cosine_similarity(csr_matrix(X))
    M = coo_matrix(S).tocsr()
    save_npz(out_npz, M)
    logging.info("D×D NPZ saved: %dx%d nnz=%d → %s", D, D, M.nnz, out_npz)

    # export labeled CSV
    df = pd.DataFrame.sparse.from_spmatrix(M, index=diseases, columns=diseases)
    csv_path = out_npz.replace(".npz", ".csv")
    df.to_csv(csv_path)
    logging.info("D×D CSV saved with labels → %s", csv_path)

    return diseases
