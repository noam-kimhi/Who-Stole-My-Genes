#!/usr/bin/env python3
"""
top_anomaly_edges.py

Extract top anomalous edges by z_robust across the whole graph or within selected components.
Writes a TSV table and an optional plot.

Inputs:
  --edges out/edge_features.tsv (must include z_robust)
  --protein_features out/protein_features.tsv (optional; used for filtering by component)
  --components e.g. "5,32" (optional; if absent, uses all edges)

Outputs:
  out_dir/top_anomaly_edges.tsv
  out_dir/top_anomaly_edges.png
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt


def parse_float(x: str) -> Optional[float]:
    x = (x or "").strip()
    if x == "":
        return None
    return float(x)


def load_node_component_map(protein_features: Path) -> Dict[str, int]:
    m: Dict[str, int] = {}
    with open(protein_features, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            m[row["u"]] = int(row["component_id"])
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", type=Path, required=True)
    ap.add_argument("--protein_features", type=Path, default=None)
    ap.add_argument("--components", type=str, default=None, help="comma-separated component IDs, optional")
    ap.add_argument("--top_n", type=int, default=20)
    ap.add_argument("--out_dir", type=Path, default=Path("figs"))
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    comp_filter = None
    node_to_comp = None
    if args.components is not None:
        if args.protein_features is None:
            raise SystemExit("If --components is set, you must also pass --protein_features.")
        comp_filter = set(int(x.strip()) for x in args.components.split(",") if x.strip())
        node_to_comp = load_node_component_map(args.protein_features)

    rows: List[Tuple[float, dict]] = []
    with open(args.edges, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if "z_robust" not in (r.fieldnames or []):
            raise SystemExit("edges must be edge_features.tsv with z_robust.")

        for row in r:
            z = parse_float(row.get("z_robust", ""))
            if z is None:
                continue

            if comp_filter is not None and node_to_comp is not None:
                u = row["u"]
                v = row["v"]
                cu = node_to_comp.get(u)
                cv = node_to_comp.get(v)
                if cu is None or cv is None or cu != cv or cu not in comp_filter:
                    continue

            rows.append((z, row))

    rows.sort(key=lambda t: t[0], reverse=True)
    top = rows[: args.top_n]

    # Write TSV
    out_tsv = args.out_dir / "top_anomaly_edges.tsv"
    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["rank", "z_robust", "jaccard", "shared_kmers", "species_u", "u", "species_v", "v"])
        for i, (z, row) in enumerate(top, 1):
            w.writerow([
                i,
                f"{z:.6f}",
                row["jaccard"],
                row["shared_kmers"],
                row.get("species_u", ""),
                row["u"],
                row.get("species_v", ""),
                row["v"],
            ])

    print(f"[OK] wrote {out_tsv}")

    # Plot
    zs = [z for z, _ in top]
    labels = [f"{i}" for i in range(1, len(zs) + 1)]

    plt.figure(figsize=(10, 6))
    plt.barh(labels[::-1], zs[::-1])  # rank 1 at top
    plt.xlabel("z_robust (surprise vs species-pair baseline)")
    title = "Top anomalous edges by z_robust"
    if comp_filter is not None:
        title += f" (components {sorted(comp_filter)})"
    plt.title(title)
    plt.tight_layout()
    out_png = args.out_dir / "top_anomaly_edges.png"
    plt.savefig(out_png, dpi=220)
    plt.close()
    print(f"[OK] wrote {out_png}")


if __name__ == "__main__":
    main()
