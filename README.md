# DataHandler

**An end-to-end pipeline** for downloading and processing HumanNet, OMIM and MeSH data into sparse gene–disease, gene–gene and disease–disease matrices (NPZ + labeled CSV).

---

## Table of Contents

1. [Overview](#overview)  
2. [Requirements](#requirements)  
3. [Installation](#installation)  
4. [Usage](#usage)  
5. [Generated Files](#generated-files)  
   1. [Downloaded Files](#downloaded-files)  
   2. [Intermediate Files](#intermediate-files)  
   3. [JSON Outputs](#json-outputs)  
   4. [Sparse Matrices](#sparse-matrices)  

---

## Overview

This toolchain performs the following steps, depositing **all** outputs under a single directory:

1. **Download HumanNet** gene–gene similarity join file  
2. **Download OMIM** `genemap2.txt` via API key  
3. **Parse OMIM** using the `genemap2-parser` into a simple JSON of gene–disease associations  
4. **Download MeSH** descriptor archive, extract the XML, and filter to the Anatomy (A) & Diseases (C) branches  
5. **Construct sparse matrices**:  
   - Gene–Disease (G×D)  
   - Gene–Gene (G×G)  
   - Disease–Disease (D×D)  
   
Each matrix is saved in both SciPy’s `.npz` format and a labeled `.csv` for easy inspection.

---

## Requirements

- Python ≥ 3.11  
- A valid OMIM API key in `OMIM_API_KEY`  
- Internet access to download data  
- [`genemap2-parser`](https://github.com/OMIM-org/genemap2-parser) on your `$PATH`  
- Dependencies listed in `pyproject.toml` / `requirements.txt`:  
  - `requests`, `numpy`, `pandas`, `scipy`, `scikit-learn`  

---

## Installation

```bash
# Clone and install
git clone https://github.com/azimgivron/datahandler.git
cd datahandler
pip install .
````

---

## Usage

Run the full pipeline with:

```bash
python main.py --output-dir output --year 2025
```

* `--output-dir` (default: `output/`): root for all files
* `--year` (default: 2025): MeSH descriptor version

---

## Generated Files

### Downloaded Files

* **`gene-gene-similarity-network.txt`**
  Raw HumanNet v1 “join” file (edge-list with per-source LLS and integrated IntNet score).

* **`genemap2.txt`**
  Raw OMIM gene–phenotype flat file from `data.omim.org`.

* **`desc{YEAR}.zip`**
  Official MeSH Descriptor archive (contains `desc{YEAR}.xml`).

* **`desc{YEAR}.xml`**
  Extracted MeSH descriptor XML, contains every MeSH term and its TreeNumbers.

### Intermediate Files

* **`genemap2_parsed/output.pickle`**
  Raw output of the `parseGeneMap2` tool (Python pickle of parsed records).

### JSON Outputs

* **`all_gene_disease.json`**
  List of `{ "mimNumber": <int>, "geneSymbols": "<A, B, …>" }` records.

* **`mesh_hierarchy.json`**
  Child→parents mapping of MeSH UIs, restricted to the Anatomy (A) & Diseases (C) branches.

* **`genes.json`**
  Ordered list of gene identifiers (Entrez symbols) used as rows for G×D & G×G.

* **`diseases.json`**
  Ordered list of OMIM MIM numbers used as columns for G×D and rows/columns for D×D.

### Sparse Matrices

For each matrix, a `.npz` (SciPy CSR) and a labeled `.csv` are produced:

* **`gene_disease.npz`** / **`gene_disease.csv`**
  Binary G×D matrix (genes × diseases) of OMIM associations.

* **`gene_gene.npz`** / **`gene_gene.csv`**
  Weighted G×G matrix (genes × genes) from HumanNet IntNet scores (self-diagonal = 1).

* **`disease_disease.npz`** / **`disease_disease.csv`**
  Cosine‐similarity D×D matrix of diseases, computed via a MeSH‐based TF×IDF pipeline.

---

**All files** end up under your specified `output/` directory, ready for downstream network analyses or machine learning workflows.
