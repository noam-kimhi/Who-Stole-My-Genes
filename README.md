# 🧬 Who Stole My Genes?

**Who Stole My Genes?** is a computational biology project that hunts for *horizontal gene transfer* (HGT) candidates across bacterial genomes — without the computational weight of whole-genome phylogenetics.  
Instead of asking _"Which species are related?"_, we ask:
- Which proteins are suspiciously similar across phylogenetically distant species?
- Can graph topology alone expose genes that crossed species boundaries?

This project was created by [**Or Forshmit**](https://github.com/OrF8), [**Noam Kimhi**](https://github.com/noam-kimhi), [**Roee Tekoah**](https://github.com/roeetekoah), [**Dor Stein**](https://github.com/dorstein0909), and [**Noam Korkos**](https://github.com/NoamKorkos)  
as part of the course [**76558 – Algorithms in Computational Biology**](https://shnaton.huji.ac.il/index.php/NewSyl/76558/1/2026/) at the Hebrew University of Jerusalem ([**HUJI**](https://en.huji.ac.il/)).

Full paper available [here](LaTeX/Who%20Stole%20My%20Genes%20-%20Detecting%20HGT%20Candidates%20Using%20Graph-Based%20Analysis.pdf).

<p align="center">
  <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/OrF8/CBIO-Hackathon-HGT/main/badges/python-percentage.json&style=round-square&logo=python" alt="python-percentage">
  <img src="https://img.shields.io/badge/BioPython-1.85-009688" alt="biopython">
  <img src="https://img.shields.io/badge/NetworkX-3.2.1-0A66C2" alt="networkx">
  <img src="https://img.shields.io/badge/MMSeqs2-external-FF6C37" alt="mmseqs2">
  <img src="https://img.shields.io/badge/NCBI_RefSeq-data-4CAF50" alt="ncbi-refseq">
</p>

---

## 🔗 Table of Contents

- [📍 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [🔬 Methodology](#-methodology)
  - [Method 1 — Orthologous Clustering](#method-1--orthologous-clustering)
  - [Method 2 — Alignment-Free k-mer Pipeline](#method-2--alignment-free-k-mer-pipeline)
- [📁 Project Structure](#-project-structure)
- [🚀 Getting Started](#-getting-started)
  - [☑️ Prerequisites](#%EF%B8%8F-prerequisites)
  - [⚙️ Installation](#%EF%B8%8F-installation)
  - [🤖 Usage](#-usage)
- [📊 Results](#-results)
- [🛠️ Technologies](#%EF%B8%8F-technologies)
- [👥 Contributors](#-contributors)
- [🎓 Course Context](#-course-context)
- [📄 License](#-license)

---

## 📍 Overview

Horizontal Gene Transfer (HGT) is one of evolution's most disruptive forces — bacteria can acquire entirely new capabilities by directly absorbing genes from unrelated organisms.  
Detecting these events computationally is hard: high similarity across distant species can arise from HGT, but also from convergent evolution, strong functional conservation, or assembly artifacts.

This project explores a **lightweight, graph-based approach** to HGT candidate detection.  
Rather than running computationally expensive whole-genome alignments or phylogenetic reconstructions, we model cross-species protein similarity as a **sparse graph** and leverage its global topology to flag anomalous edges — connections that are too strong given the taxonomic distance between their endpoints.

Two independent methodologies were developed and compared:  
- **Method 1** operates at the *protein level*, grouping proteins into homologous clusters and analyzing each cluster's internal similarity graph.  
- **Method 2** operates at the *proteome level*, using alignment-free *k*-mer Jaccard similarity to build a large cross-species graph and extract statistical outliers.

---

## ✨ Key Features

- **Two complementary HGT detection pipelines** — alignment-based and alignment-free
- **MMSeqs2 clustering** for scalable homolog grouping (Method 1)
- **k-mer inverted index** for fast candidate edge generation without alignment (Method 2)
- **Taxonomic distance integration** via the NCBI Taxonomy database
- **Graph-theoretic scoring** — suspicious edge detection, node HGT scores, betweenness, z-scores
- **3D interactive visualizations** of protein similarity graphs (Plotly)
- **Neighbor-Joining phylogenetic trees** for top HGT candidates (Method 1)
- **Simulation framework** to stress-test pipeline sensitivity against ancient/ameliorated HGTs (Method 2)

---

## 🔬 Methodology

### Method 1 — Orthologous Clustering

This method focuses on individual proteins, clustering them into homologous groups and searching for cross-species anomalies within each group.

```
Protein FASTAs (15 bacterial proteomes)
        │
        ▼
  Step 1: MMSeqs2 Clustering
  (identity ≥ 50%, coverage ≥ 80%)
        │
        ▼
  Step 2: Pairwise Alignment (BLOSUM62, local)
  Edge created if: identity ≥ 0.5 AND coverage ≥ 0.8 AND aligned_len ≥ 50
  SimilarityScore(u,v) = Identity × Coverage
        │
        ▼
  Step 3: NCBI Taxonomy Annotation
  Taxonomic distance d(u,v) = lowest shared rank index
        │
        ▼
  Step 4: HGT Scoring
  Suspicious edge: d_tax ≥ 2 AND SimilarityScore ≥ 0.6
  HGT score(i) = Σ SimilarityScore over suspicious edges
  Candidate: score ≥ P70 AND suspicious_degree ≥ 2
        │
        ▼
  Output: Top HGT candidates + 3D graph + NJ phylogenetic tree
```

### Method 2 — Alignment-Free k-mer Pipeline

This method scales to 48 species (19 families) without pairwise alignment, relying instead on *k*-mer composition statistics and graph-level anomaly detection.

```
Protein FASTAs (48 species, 19 families from RefSeq)
        │
        ▼
  Step 1: k-mer Extraction
  Build inverted index; compute shared k-mer count & Jaccard per candidate pair
        │
        ▼
  Step 2: Candidate Edge Generation
  Filter by minimum shared k-mers and top-M per protein
        │
        ▼
  Step 3: Graph Pruning
  Percentile Jaccard threshold + top-X edges per node
        │
        ▼
  Step 4: Graph Construction & Connected Components
        │
        ▼
  Step 5: Per-species-pair Robust Statistics
  Compute z-scores for each edge relative to its species-pair background
        │
        ▼
  Step 6: Node & Component Feature Extraction
  Betweenness centrality, clustering coefficient, component concentration
        │
        ▼
  Step 7: HGT Candidate Ranking
  Composite HGT-likeness score → ranked candidate list
        │
        ▼
  Output: hgt_candidates.tsv, all_scores.tsv, edge/node/component features
```

---

## 📁 Project Structure

```sh
└── CBIO-Hackathon-HGT/
    ├── Method 1/                        # Alignment-based pipeline
    │   ├── create_graph.py              # Entry point (CLI)
    │   ├── help_noam.py                 # Utility script
    │   ├── hgt_graph/                   # Core library
    │   │   ├── cli.py                   # Argument parsing & orchestration
    │   │   ├── constants.py             # Thresholds, paths, keys
    │   │   ├── graph/                   # Graph construction & scoring
    │   │   ├── io/                      # Protein I/O utilities
    │   │   ├── similarity/              # Pairwise alignment logic
    │   │   ├── taxonomy/                # NCBI taxonomy integration
    │   │   └── viz/                     # Plotly 3D & phylogenetic tree export
    │   ├── MMSeqs2 Files/               # MMSeqs2 clustering outputs
    │   ├── data/                        # Sample cluster CSV files
    │   └── results/                     # Generated graphs & phylogenetic trees
    │
    ├── Method 2/                        # Alignment-free pipeline
    │   ├── graph_hgt_pipeline.py        # Public API module
    │   ├── graph_pruning.py             # Edge pruning utilities
    │   ├── simulaiton.py                # Main simulation runner (note: filename has a typo)
    │   ├── ancient_hgt_simulation.py    # Stress-test: ancient/ameliorated HGTs
    │   ├── tax_distances.txt            # Precomputed taxonomic distances
    │   ├── config/
    │   │   └── species.txt              # 60 target bacterial species
    │   ├── src/
    │   │   ├── graph_construction/      # k-mer extraction, FASTA parsing, pruning
    │   │   └── hgt_pipeline/            # Pipeline stages, scoring, ranking
    │   ├── tools/                       # Standalone helper scripts
    │   ├── data/                        # RefSeq downloads
    │   ├── golden/                      # Ground-truth reference data
    │   ├── artifacts/                   # Intermediate pipeline artifacts
    │   └── tests/                       # Regression tests
    │
    ├── LaTeX/                           # Full academic paper (PDF + source)
    └── requirements.txt                 # Python dependencies
```

---

## 🚀 Getting Started

### ☑️ Prerequisites

- Python ≥ 3.10
- pip
- [MMSeqs2](https://github.com/soedinglab/MMseqs2) (required for Method 1 clustering step)

### ⚙️ Installation

1. Clone the repository:
```sh
git clone https://github.com/OrF8/CBIO-Hackathon-HGT
```

2. Navigate to the project directory:
```sh
cd CBIO-Hackathon-HGT
```

3. Install Python dependencies:
```sh
pip install -r requirements.txt
```

### 🤖 Usage

#### Method 1 — Alignment-Based Pipeline

Run the full pipeline on a cluster CSV file:
```sh
cd "Method 1"
python create_graph.py \
    --data data/cluster_999_size12.csv \
    --taxonomy_cache data/taxonomy_data/taxonomy_cache_cluster_999_size12.json \
    --score_percentile 70.0
```

This produces an interactive 3D Plotly graph (`results/`) and a Neighbor-Joining phylogenetic tree for the top HGT candidate.

#### Method 2 — Alignment-Free Pipeline

**Step 1 — Build candidate edges from protein FASTAs:**
```sh
cd "Method 2"
python src/graph_construction/kmer_candidates_from_faa.py \
    --manifest data/out_refseq/manifest.tsv \
    --downloads_dir data/out_refseq/downloads \
    --out candidates.tsv \
    --k 5 --min_len 50 --max_postings 2000 --min_shared 3 --top_m 50
```

**Step 2 — Prune the candidate graph:**
```sh
python graph_pruning.py --in candidates.tsv --out pruned_edges.tsv
```

**Step 3 — Run the HGT scoring pipeline:**
```sh
python graph_hgt_pipeline.py \
    --in_edges pruned_edges.tsv \
    --out_dir results/ \
    --weight_for_z jaccard \
    --z0 3.0
```

Output files in `results/`:  
| File | Description |
|---|---|
| `hgt_candidates.tsv` | Top-200 ranked HGT candidate proteins |
| `all_scores.tsv` | HGT-likeness scores for all proteins |
| `edge_features.tsv` | Per-edge z-scores and Jaccard statistics |
| `protein_features.tsv` | Per-node betweenness, clustering, component features |
| `component_features.tsv` | Per-component concentration and z-score statistics |

---

## 📊 Results

Both methods were applied to diverse bacterial datasets spanning gram-positive, gram-negative, and extremophilic organisms.

**Method 1** successfully flagged known HGT-associated genes (e.g., *eptC*, *ymdF*) and produced phylogenetic trees that visually confirm their anomalous cross-species similarity:

| EptC — Similarity Graph | EptC — Phylogenetic Tree |
|---|---|
| ![EptC graph](LaTeX/figures/EptC_graph.png) | ![EptC tree](LaTeX/figures/EptC_tree.png) |

**Method 2** demonstrated strong statistical signal in high-*z* edges and identified clusters with elevated species-pair boundary crossings. Component-level analysis revealed concentrated HGT-like components distinguishable from background noise:

| Component A                             | Component B                             | Component C                             |
|-----------------------------------------|-----------------------------------------|-----------------------------------------|
| ![compA](LaTeX/figures/component_A.png) | ![compB](LaTeX/figures/component_B.png) | ![compC](LaTeX/figures/component_C.png) |

The simulation in `ancient_hgt_simulation.py` further characterizes Method 2's sensitivity to ancient, ameliorated HGT signals versus conserved hub proteins.

---

## 🛠️ Technologies

<p align="center">
  <img src="https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white" alt="python">
  <img src="https://img.shields.io/badge/BioPython-1.85-009688" alt="biopython">
  <img src="https://img.shields.io/badge/NetworkX-3.2.1-0A66C2" alt="networkx">
  <img src="https://img.shields.io/badge/NumPy-2.0.2-013243?logo=numpy" alt="numpy">
  <img src="https://img.shields.io/badge/Pandas-2.3.3-150458?logo=pandas" alt="pandas">
  <img src="https://img.shields.io/badge/Matplotlib-3.9.4-11557C" alt="matplotlib">
  <img src="https://img.shields.io/badge/Plotly-6.5.1-3F4F75?logo=plotly" alt="plotly">
  <img src="https://img.shields.io/badge/scikit--learn-1.6.1-F7931E?logo=scikit-learn" alt="sklearn">
  <img src="https://img.shields.io/badge/SciPy-1.13.1-8CAAE6?logo=scipy" alt="scipy">
  <img src="https://img.shields.io/badge/MMSeqs2-external-FF6C37" alt="mmseqs2">
  <img src="https://img.shields.io/badge/NCBI_RefSeq-data-4CAF50" alt="ncbi">
</p>

| Tool / Library | Role |
|---|---|
| **MMSeqs2** | Ultra-fast protein sequence clustering (Method 1) |
| **BioPython** | Pairwise alignment (BLOSUM62), NCBI Entrez taxonomy queries |
| **NetworkX** | Graph construction, connected components, betweenness centrality |
| **Plotly** | Interactive 3D protein similarity graph visualization |
| **Matplotlib** | Phylogenetic tree rendering, diagnostic plots |
| **scikit-learn** | ROC-AUC evaluation in simulation benchmarks |
| **NumPy / Pandas** | Numerical computation and tabular data handling |
| **NCBI RefSeq** | Source of bacterial proteome FASTA files |
| **NCBI Taxonomy** | Taxonomic lineage annotation and distance computation |

---

## 👥 Contributors

<table>
  <tr>
    <td align="center">
      <a href="https://github.com/OrF8"><b>Or Forshmit</b></a>
    </td>
    <td align="center">
      <a href="https://github.com/noam-kimhi"><b>Noam Kimhi</b></a>
    </td>
    <td align="center">
      <a href="https://github.com/roeetekoah"><b>Roee Tekoah</b></a>
    </td>
    <td align="center">
      <a href="https://github.com/dorstein0909"><b>Dor Stein</b></a>
    </td>
    <td align="center">
      <a href="https://github.com/NoamKorkos"><b>Noam Korkos</b></a>
    </td>
  </tr>
</table>

<details closed>
<summary>Contributor Graph</summary>
<br>
<p align="left">
  <a href="https://github.com/OrF8/CBIO-Hackathon-HGT/graphs/contributors">
    <img src="https://contrib.rocks/image?repo=OrF8/CBIO-Hackathon-HGT">
  </a>
</p>
</details>

---

## 🎓 Course Context

This project was developed as part of the **Hackathon** component of:

> **76558 – Algorithms in Computational Biology**  
> [Hebrew University of Jerusalem (HUJI)](https://en.huji.ac.il/)  
> [Course Syllabus](https://shnaton.huji.ac.il/index.php/NewSyl/76558/1/2026/)

### ⚠️ Contribution Policy

This is an **academic course project**, not a community-driven open-source project.  
We are **not seeking pull requests or code contributions**.  
However, we welcome:
- Feedback on our analysis and methods
- Questions or discussion about the graph-based HGT detection approach
- Suggestions for further biological validation

---

## 📄 License

This project is licensed under the GNU General Public License v3.0.

Copyright (c) 2026 Or Forshmit, Noam Kimhi, Roee Tekoah, Dor Stein, Noam Korkos.

See the [LICENSE](LICENSE) file for details.
