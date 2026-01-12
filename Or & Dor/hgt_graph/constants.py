import os
from Bio.Align import substitution_matrices
from typing import List, Tuple

RESULTS_DIR: str = os.path.join('.', 'results')
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)
OUT_PATH_FMT: str = os.path.join(RESULTS_DIR, 'protein_similarity_graph_{}.html')
DATA_DIR: str = os.path.join('.', 'data')
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)
TAX_DATA_DIR: str = os.path.join(DATA_DIR, 'taxonomy_data')
if not os.path.exists(TAX_DATA_DIR):
    os.makedirs(TAX_DATA_DIR)
TAX_SAVE_PATH_FMT: str = os.path.join(TAX_DATA_DIR, 'taxonomy_cache_{}.json')
BLOSUM62 = substitution_matrices.load("BLOSUM62")
PROTEIN_ID_KEY: str = 'protein_id'
PROTEIN_NAME_KEY: str = 'protein_name'
ORGANISM_KEY: str = 'organism'
SEQUENCE_KEY: str = 'sequence'
SEQ_LENGTH_KEY: str = 'seq_length'
EDGE_WEIGHT_KEY: str = 'weight'
IDENTITY_KEY: str = 'identity'
COVERAGE_KEY: str = 'coverage'
ALIGNED_LENGTH_KEY: str = 'aligned_len'
TAX_DIST_KEY: str = 'taxonomic_distance'
ENTREZ_EMAIL: str = 'or.forshmit@mail.huji.ac.il'
ENTREZ_TOOL: str = 'hgt_graph'
SUS_EDGE_COLOR: str = 'rgba(220,50,50,0.95)'
REG_EDGE_COLOR: str = 'rgba(50,50,50,0.7)'
SUS_NODE_BORDER_COLOR: str = 'rgba(255,0,0,1.0)'
REG_NODE_BORDER_COLOR: str = 'rgba(30,30,30,0.6)'
GAP_OPEN_PENALTY: int = -10
GAP_EXTEND_PENALTY: int = -1
DEFAULT_MIN_ID: float = 0.4
DEFAULT_MIN_COVERAGE: float = 0.7
MIN_ALIGNED_LENGTH: int = 50
TAX_RANKS: List[str] = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
TAX_RANKS_CLOSE_TO_BROAD = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
HGT_TAX_DISTANCE_MIN: int = 2
HGT_MIN_WEIGHT: float = 0.6
TOP_HGT_N: int = 1
TOP_N_EDGES: int = 5
MIN_SUS_EDGES: int = 2
COLOR_BOUNDS: Tuple[int, int] = (30, 220)
PLOTLY_MIN_W: float = 1.0
PLOTLY_MAX_W: float = 10.0
NCBI_SLEEP: float = 0.34