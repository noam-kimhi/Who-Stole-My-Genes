import pandas as pd
import argparse
from pathlib import Path
from typing import Set, Tuple
from .constants import OUT_PATH_FMT, TAX_SAVE_PATH_FMT, HGT_TAX_DISTANCE_MIN, HGT_MIN_WEIGHT, TOP_HGT_N, ORGANISM_KEY
from .io.proteins import create_node_attributes, extract_seqs
from .graph.build import build_nx_graph
from .similarity.edges import build_edges_from_alignments
from .taxonomy.cache import get_or_build_taxonomy_cache
from .taxonomy.annotate import annotate_graph_with_taxonomy, annotate_edges_with_tax_distance
from .graph.scoring import compute_hgt_scores, percentile, top_hgt_candidates, top_hgt_edges
from .viz.plotly3d import export_plotly_3d

def main() -> None:
    parser = argparse.ArgumentParser(description="Create protein similarity graph and identify HGT candidates.")
    parser.add_argument('--data', type=str, help='Path to input CSV data file.', required=True)
    parser.add_argument('--taxonomy_cache', type=str, help='Path to taxonomy cache JSON file.')
    parser.add_argument("--score_percentile", type=float, default=70.0,
                        help="Percentile for HGT score threshold (default 70.0).")

    args = parser.parse_args()

    data_path: Path = Path(args.data)
    path_suffix: str = data_path.stem

    df: pd.DataFrame = pd.read_csv(data_path)

    out_path: str = OUT_PATH_FMT.format(path_suffix)

    node_attrs = create_node_attributes(df)
    seqs = extract_seqs(df)

    edges = build_edges_from_alignments(seqs)

    G = build_nx_graph(node_attrs, edges)

    taxonomy_cache = get_or_build_taxonomy_cache(
        args.taxonomy_cache,
        node_attrs=node_attrs,
        save_path=TAX_SAVE_PATH_FMT.format(path_suffix),
    )

    annotate_graph_with_taxonomy(G, taxonomy_cache)
    annotate_edges_with_tax_distance(G)

    scores, suspicious_degree, suspicious_edges = compute_hgt_scores(
        G,
        min_tax_distance=HGT_TAX_DISTANCE_MIN,
        min_weight=HGT_MIN_WEIGHT
    )

    nonzero_scores = [s for s in scores.values() if s > 0]
    if not nonzero_scores:
        print("No suspicious edges found under current thresholds.")
        top = []
        highlight_nodes = set()
    else:
        p: float = args.score_percentile
        score_threshold = percentile(nonzero_scores, p)
        top = top_hgt_candidates(G, scores, suspicious_degree, score_threshold, TOP_HGT_N)
        highlight_nodes = {node for node, _ in top}

        print(f"Score threshold (P{p} of non-zero): {score_threshold:.4f}")

    suspicious_edges_display: Set[Tuple[str, str]] = top_hgt_edges(suspicious_edges, scores)

    print("\nTop HGT candidates:")
    for node, score in top:
        org = G.nodes[node].get(ORGANISM_KEY, "")
        print(f"{node}\tHGT_score={score:.4f}\t{org}")

    export_plotly_3d(G, out_path, suspicious_edges_display, scores, highlight_nodes)


if __name__ == "__main__":
    main()
