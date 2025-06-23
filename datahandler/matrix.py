"""
Module for constructing coherent gene–disease (GxD), gene–gene (GxG), and disease–disease (DxD) matrices
from OMIM, HumanNet, and MeSH hierarchy data.

Note on TFxIDF:
We do not use sklearn.feature_extraction.text.TfidfVectorizer directly because:
  1. Controlled Vocabulary: we extract curated MeSH concepts (multi-word), not raw tokens.
  2. Hierarchy Propagation: counts are refined via division by number of hyponyms per Eq.(1),
     which TfidfVectorizer cannot support.
  3. Custom IDF: IDF is computed on refined MeSH counts, not raw token counts.
  4. Matrix Orientation: we require GxD and DxD sparse matrices; TfidfVectorizer produces DxT.

This pipeline reproduces the MimMiner disease–disease similarity and coherent GxD, GxG
matrices as described in the target publication.

See for more details:
    *   Driel, Marc & Bruggeman, Jorn & Vriend, Gert & Brunner, Han. (2006).
        A text-mining analysis of the human phenome <http://www.cmbi.ru.nl/MimMiner.
        European journal of human genetics : EJHG. 14. 535-42. 10.1038/sj.ejhg.5201585.
"""
import logging
import json
import numpy as np
from collections import defaultdict
from scipy.sparse import coo_matrix, save_npz, csr_matrix
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
    with open(associations_file) as f:
        assoc = json.load(f)
    genes_set, diseases_set, records = set(), set(), []
    for entry in assoc:
        d = str(entry.get('mimNumber'))
        diseases_set.add(d)
        for g in entry.get('geneSymbols', '').split(', '):
            if g:
                genes_set.add(g)
                records.append((g, d))
    genes = sorted(genes_set)
    diseases = sorted(diseases_set)
    return records, genes, diseases

# --- Gene–Disease Matrix (GxD) ---
def construct_gene_disease_matrix(associations_file: str, out_npz: str):
    """
    Build and save a binary gene–disease sparse matrix (GxD).

    The resulting coo_matrix has shape (G, D) where G = #genes, D = #diseases.
    Returns the genes and diseases lists for consistent ordering.

    Note:
      We use a sparse COO to store only nonzero associations, rather than dense arrays.
    """
    records, genes, diseases = load_gene_disease_associations(associations_file)
    idx_g = {g: i for i, g in enumerate(genes)}
    idx_d = {d: i for i, d in enumerate(diseases)}
    rows, cols, data, seen = [], [], [], set()
    for g, d in records:
        i, j = idx_g[g], idx_d[d]
        if (i, j) not in seen:
            rows.append(i); cols.append(j); data.append(1)
            seen.add((i, j))
    mat = coo_matrix((data, (rows, cols)), shape=(len(genes), len(diseases)), dtype=int)
    save_npz(out_npz, mat)
    logging.info("Gene–disease matrix saved: %dx%d, nnz=%d -> %s",
                 mat.shape[0], mat.shape[1], mat.nnz, out_npz)
    return genes, diseases

# --- Gene–Gene Matrix (GxG) ---
def construct_gene_gene_matrix(network_file: str, out_npz: str, genes):
    """
    Build and save a gene–gene similarity sparse matrix (GxG) filtered to given genes.

    Input edge list: lines of 'gene1 gene2 weight'.
    Ensures diagonal self-similarity = 1.
    """
    idx_g = {g: i for i, g in enumerate(genes)}
    G = len(genes)
    rows, cols, data = [], [], []
    with open(network_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3: continue
            g1, g2, w = parts[0], parts[1], float(parts[2])
            if g1 in idx_g and g2 in idx_g:
                rows.append(idx_g[g1]); cols.append(idx_g[g2]); data.append(w)
    # self-similarity
    for i in range(G):
        rows.append(i); cols.append(i); data.append(1.0)
    mat = coo_matrix((data, (rows, cols)), shape=(G, G))
    save_npz(out_npz, mat)
    logging.info("Gene–gene matrix saved: %dx%d, nnz=%d -> %s",
                 G, G, mat.nnz, out_npz)
    return genes

# --- Disease–Disease Matrix (DxD) via MimMiner pipeline ---

def load_omim_entries(omim_json: str) -> dict:
    """
    Load OMIM entries JSON mapping MIM numbers to dicts with 'TX' and 'CS' fields.
    """
    with open(omim_json) as f:
        return json.load(f)


def load_mesh_hierarchy(mesh_hierarchy_json: str) -> dict:
    """
    Load MeSH hierarchy: child term → list of parent terms.
    """
    with open(mesh_hierarchy_json) as f:
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
    text_up = text.upper()
    return [t for t in vocab if t.upper() in text_up]


def construct_disease_disease_matrix(
    omim_json: str,
    mesh_hierarchy_json: str,
    out_npz: str,
    diseases
) -> list:
    """
    Build and save a disease–disease cosine similarity sparse matrix (DxD).

    Pipeline:
      1. Raw counts of MeSH in TX+CS.
      2. Hierarchy-aware refinement per Eq.(1): divide hyponym counts by n_h and add to parent.
      3. Compute IDF on refined counts (Eq.(2)).
      4. Form TFxIDF vectors, L2-normalize each (Eq.(3)).
      5. Compute full cosine similarity (Eq.(4)).
    Requires diseases list for consistent ordering.
    """
    records   = load_omim_entries(omim_json)
    hierarchy = load_mesh_hierarchy(mesh_hierarchy_json)
    children  = build_children_map(hierarchy)
    vocab     = set(hierarchy.keys())
    D = len(diseases)

    # Steps 1-2: raw + refined counts
    term_doc_freq = defaultdict(int)
    refined = {}
    for mim in diseases:
        entry = records.get(str(mim), {})
        text  = (entry.get('TX','') + ' ' + entry.get('CS','')).upper()
        raw   = defaultdict(int)
        for t in extract_mesh_terms(text, vocab): raw[t] += 1
        ref   = raw.copy()
        for p, hypos in children.items():
            nh = len(hypos)
            if nh:
                s = sum(raw.get(h,0) for h in hypos)
                if s: ref[p] += s/nh
        refined[mim] = ref
        for t, c in ref.items():
            if c>0: term_doc_freq[t]+=1

    # Step 3: IDF
    idf       = {t: np.log(D/df) for t, df in term_doc_freq.items()}
    terms     = sorted(idf)

    # Step 4: TFxIDF and L2 norm → sparse term matrix
    rows, cols, data = [], [], []
    for i, mim in enumerate(diseases):
        vec = np.array([refined[mim].get(t,0)*idf[t] for t in terms])
        norm= np.linalg.norm(vec)
        if norm: vec/=norm
        for j in np.nonzero(vec)[0]:
            rows.append(i); cols.append(j); data.append(vec[j])
    X = coo_matrix((data,(rows,cols)),shape=(D,len(terms)))

    # Step 5: cosine similarity → DxD
    S = cosine_similarity(csr_matrix(X))
    M = coo_matrix(S)
    save_npz(out_npz, M)
    logging.info("Disease–disease matrix saved: %dx%d, nnz=%d -> %s",
                 D, D, M.nnz, out_npz)
    return diseases
