#!/usr/bin/env python3
"""
explain_top_candidates.py

Explain top HGT-like proteins in a human-readable, text-only report.

Modes:
  - Step 2 (recommended): provide --hgt_candidates out/hgt_candidates.tsv
      -> uses the pipeline's final ranking + score (component-relative, log-space)
  - Step 1 fallback: if --hgt_candidates not provided
      -> ranks by a simple standalone score from protein_features.tsv

Inputs:
  --edges              Prefer out/edge_features.tsv so neighbor z_robust is available
  --protein_features   out/protein_features.tsv
  --component_features out/component_features.tsv
  --hgt_candidates     out/hgt_candidates.tsv (optional but recommended)

Output:
  Prints per-candidate:
    - pipeline score (if provided) or standalone score
    - node features
    - component summary
    - neighbor species concentration
    - top neighbors by (z, jaccard)
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict, Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def sniff_delimiter(path: Path) -> str:
    with open(path, "r", encoding="utf-8", newline="") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "\t" in line and "," not in line:
                return "\t"
            if "," in line and "\t" not in line:
                return ","
            return ","
    return ","


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
class AdjEdge:
    v: str
    sp_v: str
    shared: int
    jac: float
    z: Optional[float]


def parse_float(x: str) -> Optional[float]:
    x = (x or "").strip()
    if x == "":
        return None
    return float(x)


def load_protein_features(path: Path) -> Dict[str, NodeRow]:
    m: Dict[str, NodeRow] = {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            u = row["u"]
            m[u] = NodeRow(
                u=u,
                species=row["species"],
                component_id=int(row["component_id"]),
                deg_xsp=int(row["deg_xsp"]),
                n_species=int(row["n_species"]),
                max_z=parse_float(row["max_z"]),
                top5_mean_z=parse_float(row["top5_mean_z"]),
                max_species_fraction=float(row["max_species_fraction"]),
                entropy_norm=float(row["species_entropy_norm"]),
                betweenness=parse_float(row.get("betweenness", "")),
                clustering=parse_float(row.get("clustering_coeff", "")),
            )
    return m


def load_component_features(path: Path) -> Dict[int, CompRow]:
    m: Dict[int, CompRow] = {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            cid = int(row["component_id"])
            m[cid] = CompRow(
                component_id=cid,
                size=int(row["size"]),
                num_edges=int(row["num_edges"]),
                n_species=int(row["n_species"]),
                entropy_norm=float(row["species_entropy_norm"]),
                high_z_frac=parse_float(row.get("high_z_frac", "")),
                max_z=parse_float(row.get("max_z", "")),
            )
    return m


def load_adjacency(edges_path: Path) -> Dict[str, List[AdjEdge]]:
    """
    Build adjacency list from edges file.

    If you pass out/edge_features.tsv, it includes z_robust and the script will show it.
    """
    delim = sniff_delimiter(edges_path)
    adj: Dict[str, List[AdjEdge]] = defaultdict(list)

    with open(edges_path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter=delim)
        has_z = "z_robust" in (r.fieldnames or [])

        for row in r:
            u = row["u"]
            v = row["v"]
            sp_u = row.get("species_u", "")
            sp_v = row.get("species_v", "")
            shared = int(row["shared_kmers"])
            jac = float(row["jaccard"])
            z = parse_float(row["z_robust"]) if has_z else None

            adj[u].append(AdjEdge(v=v, sp_v=sp_v, shared=shared, jac=jac, z=z))
            adj[v].append(AdjEdge(v=u, sp_v=sp_u, shared=shared, jac=jac, z=z))

    return adj


def standalone_score(n: NodeRow) -> float:
    """
    Step-1 style heuristic score (legacy):
      score = top5_mean_z * log(1+deg) * max_species_fraction
    """
    import math
    z = n.top5_mean_z or 0.0
    if z <= 0:
        return 0.0
    return z * math.log1p(n.deg_xsp) * n.max_species_fraction


def load_hgt_candidates(path: Path, top_n: int) -> List[Tuple[str, float]]:
    """
    Read top-N (u, score) from hgt_candidates.tsv in file order.
    Expected columns include: u, score
    """
    out: List[Tuple[str, float]] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if "u" not in (r.fieldnames or []) or "score" not in (r.fieldnames or []):
            raise ValueError(f"{path} must include columns 'u' and 'score'. Found: {r.fieldnames}")

        for row in r:
            u = row["u"]
            s = float(row["score"])
            out.append((u, s))
            if len(out) >= top_n:
                break
    return out


def explain_one(
    u: str,
    node: NodeRow,
    comp: Optional[CompRow],
    adj: Dict[str, List[AdjEdge]],
    top_k_neighbors: int,
) -> str:
    lines: List[str] = []
    lines.append(f"Protein: {u}")
    lines.append(f"  species={node.species}  component={node.component_id}  deg={node.deg_xsp}  n_species={node.n_species}")
    lines.append(
        f"  top5_mean_z={node.top5_mean_z}  max_z={node.max_z}  "
        f"max_species_fraction={node.max_species_fraction:.3f}  entropy_norm={node.entropy_norm:.3f}"
    )
    lines.append(f"  betweenness={node.betweenness}  clustering={node.clustering}")

    if comp is None:
        lines.append("  component_summary: (missing)")
    else:
        lines.append(
            f"  component_summary: size={comp.size} edges={comp.num_edges} n_species={comp.n_species} "
            f"H_norm={comp.entropy_norm:.3f} high_z_frac={comp.high_z_frac} max_z={comp.max_z}"
        )

    neigh = adj.get(u, [])
    if not neigh:
        lines.append("  neighbors: (none found in edges file)")
        return "\n".join(lines)

    sp_counts = Counter(e.sp_v for e in neigh)
    top_sp, top_c = sp_counts.most_common(1)[0]
    lines.append(f"  neighbor_species_top: {top_sp}  ({top_c}/{len(neigh)} = {top_c/len(neigh):.3f})")

    # Sort neighbors: by z desc (if present), then by jaccard desc
    if any(e.z is not None for e in neigh):
        neigh_sorted = sorted(
            neigh,
            key=lambda e: ((e.z if e.z is not None else float("-inf")), e.jac),
            reverse=True,
        )
    else:
        neigh_sorted = sorted(neigh, key=lambda e: e.jac, reverse=True)

    lines.append(f"  top_neighbors (k={top_k_neighbors}):")
    for e in neigh_sorted[:top_k_neighbors]:
        lines.append(f"    v={e.v}  sp={e.sp_v}  jac={e.jac:.6f}  shared={e.shared}  z={e.z}")

    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", required=True, type=Path, help="Edges file (prefer out/edge_features.tsv)")
    ap.add_argument("--protein_features", required=True, type=Path, help="out/protein_features.tsv")
    ap.add_argument("--component_features", required=True, type=Path, help="out/component_features.tsv")
    ap.add_argument("--hgt_candidates", type=Path, default=None, help="out/hgt_candidates.tsv (recommended)")
    ap.add_argument("--top_n", type=int, default=20)
    ap.add_argument("--top_k_neighbors", type=int, default=10)
    args = ap.parse_args()

    nodes = load_protein_features(args.protein_features)
    comps = load_component_features(args.component_features)
    adj = load_adjacency(args.edges)

    if args.hgt_candidates is not None:
        ranked_u = load_hgt_candidates(args.hgt_candidates, args.top_n)
        print(f"[REPORT] Top {args.top_n} candidates from {args.hgt_candidates.name} (pipeline Step-3 ranking)")
        print(f"[NOTE] Neighbor z values shown only if --edges points to edge_features.tsv.")
        print()
        for i, (u, score) in enumerate(ranked_u, 1):
            node = nodes.get(u)
            if node is None:
                print(f"=== #{i} score={score:.6f} ===")
                print(f"Protein: {u}")
                print("  (missing from protein_features.tsv)")
                print()
                continue
            comp = comps.get(node.component_id)
            print(f"=== #{i} score={score:.6f} ===")
            print(explain_one(u, node, comp, adj, args.top_k_neighbors))
            print()
    else:
        # Legacy fallback
        ranked = sorted(nodes.values(), key=standalone_score, reverse=True)[:args.top_n]
        print(f"[REPORT] Top {args.top_n} candidates (standalone Step-1 score) using {args.protein_features.name}")
        print(f"[NOTE] For pipeline Step-2 ranking, pass --hgt_candidates out/hgt_candidates.tsv.")
        print()
        for i, n in enumerate(ranked, 1):
            comp = comps.get(n.component_id)
            print(f"=== #{i} score={standalone_score(n):.6f} ===")
            print(explain_one(n.u, n, comp, adj, args.top_k_neighbors))
            print()


if __name__ == "__main__":
    main()
