#!/usr/bin/env python3
"""
explain_component.py

Text-only explanation of a specific connected component for report writing.

Inputs:
  --edges              out/edge_features.tsv (recommended; must include z_robust if you want z-based edge ranking)
  --protein_features   out/protein_features.tsv
  --component_features out/component_features.tsv
  --hgt_candidates     out/hgt_candidates.tsv (optional but recommended; provides final scores)

Usage examples:
  python explain_component.py --component_id 5 --edges out/edge_features.tsv --protein_features out/protein_features.tsv --component_features out/component_features.tsv --hgt_candidates out/hgt_candidates.tsv

  python explain_component.py --component_id 5 --top_nodes 20 --top_edges 30 --edges out/edge_features.tsv --protein_features out/protein_features.tsv --component_features out/component_features.tsv --hgt_candidates out/hgt_candidates.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_float(x: str) -> Optional[float]:
    x = (x or "").strip()
    if x == "":
        return None
    return float(x)


@dataclass
class NodeRow:
    u: str
    species: str
    component_id: int
    deg_xsp: int
    n_species: int
    max_z: Optional[float]
    top5_mean_z: Optional[float]
    max_species_fraction: float
    entropy_norm: float
    betweenness: Optional[float]
    clustering: Optional[float]


@dataclass
class CompRow:
    component_id: int
    size: int
    num_edges: int
    n_species: int
    entropy_norm: float
    high_z_frac: Optional[float]
    max_z: Optional[float]


@dataclass
class EdgeRow:
    u: str
    v: str
    species_u: str
    species_v: str
    shared: int
    jac: float
    z: Optional[float]


def load_protein_features(path: Path) -> Dict[str, NodeRow]:
    out: Dict[str, NodeRow] = {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            u = row["u"]
            out[u] = NodeRow(
                u=u,
                species=row["species"],
                component_id=int(row["component_id"]),
                deg_xsp=int(row["deg_xsp"]),
                n_species=int(row["n_species"]),
                max_z=parse_float(row.get("max_z", "")),
                top5_mean_z=parse_float(row.get("top5_mean_z", "")),
                max_species_fraction=float(row["max_species_fraction"]),
                entropy_norm=float(row["species_entropy_norm"]),
                betweenness=parse_float(row.get("betweenness", "")),
                clustering=parse_float(row.get("clustering_coeff", "")),
            )
    return out


def load_component_features(path: Path) -> Dict[int, CompRow]:
    out: Dict[int, CompRow] = {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            cid = int(row["component_id"])
            out[cid] = CompRow(
                component_id=cid,
                size=int(row["size"]),
                num_edges=int(row["num_edges"]),
                n_species=int(row["n_species"]),
                entropy_norm=float(row["species_entropy_norm"]),
                high_z_frac=parse_float(row.get("high_z_frac", "")),
                max_z=parse_float(row.get("max_z", "")),
            )
    return out


def load_scores(hgt_candidates_path: Optional[Path]) -> Dict[str, float]:
    """
    Load final scores from hgt_candidates.tsv if provided.
    Returns u -> score
    """
    if hgt_candidates_path is None:
        return {}
    out: Dict[str, float] = {}
    with open(hgt_candidates_path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if "u" not in (r.fieldnames or []) or "score" not in (r.fieldnames or []):
            raise ValueError(f"{hgt_candidates_path} must contain columns: u, score")
        for row in r:
            out[row["u"]] = float(row["score"])
    return out


def load_component_subgraph_edges(
    edges_path: Path,
    nodes_in_component: set[str],
) -> List[EdgeRow]:
    """
    Load only edges whose endpoints are both in nodes_in_component.
    Expects edge_features.tsv columns:
      u v species_u species_v shared_kmers jaccard z_robust
    """
    edges: List[EdgeRow] = []
    with open(edges_path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fields = set(r.fieldnames or [])
        required = {"u", "v", "shared_kmers", "jaccard", "species_u", "species_v"}
        if not required.issubset(fields):
            raise ValueError(f"{edges_path} missing required columns. Found: {r.fieldnames}")
        has_z = "z_robust" in fields

        for row in r:
            u = row["u"]
            v = row["v"]
            if u not in nodes_in_component or v not in nodes_in_component:
                continue
            edges.append(
                EdgeRow(
                    u=u,
                    v=v,
                    species_u=row.get("species_u", ""),
                    species_v=row.get("species_v", ""),
                    shared=int(row["shared_kmers"]),
                    jac=float(row["jaccard"]),
                    z=parse_float(row["z_robust"]) if has_z else None,
                )
            )
    return edges


def density(num_nodes: int, num_edges: int) -> float:
    if num_nodes <= 1:
        return 0.0
    return (2.0 * num_edges) / (num_nodes * (num_nodes - 1.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--component_id", type=int, required=True)
    ap.add_argument("--edges", type=Path, required=True, help="out/edge_features.tsv (tab-separated)")
    ap.add_argument("--protein_features", type=Path, required=True)
    ap.add_argument("--component_features", type=Path, required=True)
    ap.add_argument("--hgt_candidates", type=Path, default=None, help="out/hgt_candidates.tsv (optional)")
    ap.add_argument("--top_nodes", type=int, default=15)
    ap.add_argument("--top_edges", type=int, default=20)
    ap.add_argument("--top_species", type=int, default=15)
    args = ap.parse_args()

    nodes = load_protein_features(args.protein_features)
    comps = load_component_features(args.component_features)
    scores = load_scores(args.hgt_candidates)

    cid = args.component_id
    comp = comps.get(cid)
    if comp is None:
        raise SystemExit(f"Component {cid} not found in {args.component_features}")

    # collect nodes in the component
    comp_nodes = [u for u, n in nodes.items() if n.component_id == cid]
    comp_set = set(comp_nodes)

    # load edges within component
    comp_edges = load_component_subgraph_edges(args.edges, comp_set)

    # sanity
    dens = density(comp.size, comp.num_edges)

    print(f"[COMPONENT] id={cid}")
    print(f"  size={comp.size}  edges={comp.num_edges}  density={dens:.4f}")
    print(f"  n_species={comp.n_species}  H_norm={comp.entropy_norm:.3f}  high_z_frac={comp.high_z_frac}  max_z={comp.max_z}")
    print()

    # species distribution
    sp_counts = Counter(nodes[u].species for u in comp_nodes)
    print("[SPECIES] top species counts:")
    for sp, c in sp_counts.most_common(args.top_species):
        print(f"  {sp:<30} {c}")
    print()

    # node rankings
    def get_score(u: str) -> float:
        return scores.get(u, 0.0)

    # top by final score (if available)
    if scores:
        top_by_score = sorted(comp_nodes, key=get_score, reverse=True)[:args.top_nodes]
        print(f"[NODES] top by FINAL score (from hgt_candidates.tsv), k={args.top_nodes}:")
        for u in top_by_score:
            n = nodes[u]
            print(
                f"  score={get_score(u):9.6f}  "
                f"deg={n.deg_xsp:2d}  ns={n.n_species:2d}  "
                f"top5z={'' if n.top5_mean_z is None else f'{n.top5_mean_z:.3f}':>8}  "
                f"bw={'' if n.betweenness is None else f'{n.betweenness:.6g}':>10}  "
                f"clust={'' if n.clustering is None else f'{n.clustering:.3f}':>6}  "
                f"maxSpFrac={n.max_species_fraction:.3f}  "
                f"{u}"
            )
        print()

    # top by top5_mean_z
    top_by_z = sorted(comp_nodes, key=lambda u: (nodes[u].top5_mean_z or float("-inf")), reverse=True)[:args.top_nodes]
    print(f"[NODES] top by top5_mean_z, k={args.top_nodes}:")
    for u in top_by_z:
        n = nodes[u]
        print(
            f"  top5z={'' if n.top5_mean_z is None else f'{n.top5_mean_z:.3f}':>8}  "
            f"deg={n.deg_xsp:2d}  ns={n.n_species:2d}  "
            f"bw={'' if n.betweenness is None else f'{n.betweenness:.6g}':>10}  "
            f"clust={'' if n.clustering is None else f'{n.clustering:.3f}':>6}  "
            f"{u}"
        )
    print()

    # top by betweenness
    top_by_bw = sorted(comp_nodes, key=lambda u: (nodes[u].betweenness or 0.0), reverse=True)[:args.top_nodes]
    print(f"[NODES] top by betweenness, k={args.top_nodes}:")
    for u in top_by_bw:
        n = nodes[u]
        print(
            f"  bw={'' if n.betweenness is None else f'{n.betweenness:.6g}':>10}  "
            f"clust={'' if n.clustering is None else f'{n.clustering:.3f}':>6}  "
            f"deg={n.deg_xsp:2d}  ns={n.n_species:2d}  "
            f"top5z={'' if n.top5_mean_z is None else f'{n.top5_mean_z:.3f}':>8}  "
            f"{u}"
        )
    print()

    # top by low clustering (bridge proxy)
    def bridge_proxy(u: str) -> float:
        c = nodes[u].clustering
        if c is None:
            return 0.0
        return 1.0 - c

    top_by_bridge = sorted(comp_nodes, key=bridge_proxy, reverse=True)[:args.top_nodes]
    print(f"[NODES] top by (1 - clustering) bridge proxy, k={args.top_nodes}:")
    for u in top_by_bridge:
        n = nodes[u]
        bp = bridge_proxy(u)
        print(
            f"  bridge={bp:6.3f}  "
            f"clust={'' if n.clustering is None else f'{n.clustering:.3f}':>6}  "
            f"deg={n.deg_xsp:2d}  ns={n.n_species:2d}  "
            f"top5z={'' if n.top5_mean_z is None else f'{n.top5_mean_z:.3f}':>8}  "
            f"{u}"
        )
    print()

    # edge rankings
    if not comp_edges:
        print("[EDGES] no edges found for this component (unexpected).")
        return

    has_any_z = any(e.z is not None for e in comp_edges)
    if has_any_z:
        comp_edges_sorted = sorted(comp_edges, key=lambda e: (e.z if e.z is not None else float("-inf")), reverse=True)
        print(f"[EDGES] top by z_robust, k={args.top_edges}:")
        for e in comp_edges_sorted[:args.top_edges]:
            print(
                f"  z={e.z:.6g}  jac={e.jac:.6f}  shared={e.shared:4d}  "
                f"{e.species_u} | {e.u}   --   {e.species_v} | {e.v}"
            )
        print()

    comp_edges_sorted_j = sorted(comp_edges, key=lambda e: e.jac, reverse=True)
    print(f"[EDGES] top by jaccard, k={args.top_edges}:")
    for e in comp_edges_sorted_j[:args.top_edges]:
        print(
            f"  jac={e.jac:.6f}  shared={e.shared:4d}  z={e.z}  "
            f"{e.species_u} | {e.u}   --   {e.species_v} | {e.v}"
        )
    print()


if __name__ == "__main__":
    main()
