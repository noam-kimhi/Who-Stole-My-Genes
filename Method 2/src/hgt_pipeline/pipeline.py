#!/usr/bin/env python3
"""End-to-end graph HGT pipeline orchestration."""

from __future__ import annotations

import argparse
import time
from contextlib import contextmanager
from pathlib import Path

from .stages.component_features import (
    compute_component_features,
    write_component_features,
)
from .stages.edge_io import dedupe_edges, read_edges
from .stages.graph_ops import (
    attach_edge_z_to_graph,
    build_graph,
    compute_components,
)
from .stages.node_features import compute_node_features, write_node_features
from .stages.pair_stats import (
    compute_edge_features,
    compute_pair_robust_stats,
    print_pair_stats_sanity_table,
    write_edge_features,
)
from .stages.ranking import (
    score_hgt_likeness,
    write_all_scores,
    write_hgt_candidates,
)


@contextmanager
def stage(name: str):
    """Context manager to time a pipeline stage."""
    t0 = time.perf_counter()
    print(f"[STAGE] {name} ...")
    try:
        yield
    finally:
        dt = time.perf_counter() - t0
        print(f"[DONE ] {name} in {dt:.3f}s")


def log_counts(edges=None, graph=None, pair_stats=None, prefix: str = "[INFO]"):
    """Print quick counts for sanity."""
    parts = []
    if edges is not None:
        parts.append(f"edges={len(edges)}")
    if graph is not None:
        parts.append(f"nodes={graph.number_of_nodes()} edges={graph.number_of_edges()}")
    if pair_stats is not None:
        parts.append(f"species_pairs={len(pair_stats)}")
    if parts:
        print(f"{prefix} " + " | ".join(parts))


def run_pipeline(
    in_path: Path,
    out_dir: Path,
    weight_for_z: str = "jaccard",
    min_pair_edges_for_z: int = 30,
    z0: float = 3.0,
    compute_betweenness: bool = True,
) -> None:
    """Run the full graph-analysis pipeline."""
    out_dir.mkdir(parents=True, exist_ok=True)

    with stage("1) Read edges"):
        edges = dedupe_edges(read_edges(in_path))
        log_counts(edges=edges)

    with stage("2) Build graph"):
        graph = build_graph(edges, use_networkx=True)
        log_counts(graph=graph)

    with stage("3) Robust per-species-pair stats"):
        pair_stats = compute_pair_robust_stats(
            edges,
            weight=weight_for_z,
            min_pair_edges=min_pair_edges_for_z,
        )
        log_counts(pair_stats=pair_stats)
        print_pair_stats_sanity_table(
            pair_stats,
            min_pair_edges_for_z=min_pair_edges_for_z,
            top_n=15,
        )

    with stage("4) Edge features + write edge_features.tsv"):
        edge_feats = compute_edge_features(
            edges,
            pair_stats=pair_stats,
            weight=weight_for_z,
            min_pair_edges_for_z=min_pair_edges_for_z,
        )
        write_edge_features(out_dir / "edge_features.tsv", edge_feats)

    with stage("5) Attach z to graph edges"):
        attach_edge_z_to_graph(graph, edge_feats)

    total_edges = graph.number_of_edges()
    with_z = sum(1 for _, _, d in graph.edges(data=True) if "z_robust" in d)
    print(f"[SANITY] edges with z: {with_z}/{total_edges}")

    with stage("6) Connected components"):
        component_id_of, components = compute_components(graph)
        print(f"[INFO] components={len(components)}")

    with stage("7) Component features + write component_features.tsv"):
        comp_feats = compute_component_features(graph, components, z0=z0)
        write_component_features(out_dir / "component_features.tsv", comp_feats)

    with stage("8) Node features + write protein_features.tsv"):
        node_feats = compute_node_features(
            graph,
            component_id_of=component_id_of,
            compute_betweenness=compute_betweenness,
            betweenness_mode="high_z",
            z_threshold_for_high_subgraph=z0,
        )
        write_node_features(out_dir / "protein_features.tsv", node_feats)

    with stage("9) Rank proteins + write candidate tables"):
        ranked = score_hgt_likeness(
            node_feats,
            comp_feats,
            min_component_size=10,
            use_betweenness=compute_betweenness,
        )
        scores = {u: s for u, s in ranked}
        write_hgt_candidates(out_dir / "hgt_candidates.tsv", node_feats, scores, top_n=200)
        write_all_scores(out_dir / "all_scores.tsv", node_feats, scores)

    if ranked:
        print("[SANITY] top 10 HGT-like proteins:")
        for u, s in ranked[:10]:
            print(f"  {s:10.3f}  {u}")

        vals = sorted(s for _, s in ranked if s > 0)
        if vals:
            print(
                f"[SANITY] score stats: n={len(vals)} "
                f"min={vals[0]:.6g} p50={vals[len(vals)//2]:.6g} "
                f"p95={vals[int(0.95 * len(vals))]:.6g} max={vals[-1]:.6g}"
            )


def main() -> None:
    ap = argparse.ArgumentParser(description="Graph construction + HGT-like analysis.")
    ap.add_argument("--in_edges", required=True, type=Path, help="Input CSV/TSV of undirected pruned edges")
    ap.add_argument("--out_dir", required=True, type=Path, help="Output directory for TSV deliverables")
    ap.add_argument("--weight_for_z", choices=["jaccard", "shared"], default="jaccard")
    ap.add_argument("--min_pair_edges_for_z", type=int, default=30)
    ap.add_argument("--z0", type=float, default=3.0, help="High-z threshold for component concentration + subgraph")
    ap.add_argument("--no_betweenness", action="store_true", help="Skip betweenness for faster runs")
    args = ap.parse_args()

    run_pipeline(
        in_path=args.in_edges,
        out_dir=args.out_dir,
        weight_for_z=args.weight_for_z,
        min_pair_edges_for_z=args.min_pair_edges_for_z,
        z0=args.z0,
        compute_betweenness=not args.no_betweenness,
    )


if __name__ == "__main__":
    main()
