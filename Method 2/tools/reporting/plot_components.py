#!/usr/bin/env python3
"""
plot_components.py

Unified component plotting script.

Key options:
  - --node_size_mode score|constant
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


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


@dataclass
class EdgeRow:
    u: str
    v: str
    jac: float
    shared: int
    z: Optional[float]


def load_nodes(protein_features: Path) -> Dict[str, NodeRow]:
    out: Dict[str, NodeRow] = {}
    with open(protein_features, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            u = row["u"]
            out[u] = NodeRow(u=u, species=row["species"], component_id=int(row["component_id"]))
    return out


def load_scores(hgt_candidates: Optional[Path]) -> Dict[str, float]:
    if hgt_candidates is None:
        return {}
    out: Dict[str, float] = {}
    with open(hgt_candidates, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            out[row["u"]] = float(row["score"])
    return out


def load_edges_for_component(edges_path: Path, nodes_in_comp: set[str]) -> List[EdgeRow]:
    edges: List[EdgeRow] = []
    with open(edges_path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fields = set(r.fieldnames or [])
        if "z_robust" not in fields:
            raise ValueError("edges file must include z_robust (use edge_features.tsv).")
        for row in r:
            u = row["u"]
            v = row["v"]
            if u not in nodes_in_comp or v not in nodes_in_comp:
                continue
            edges.append(
                EdgeRow(
                    u=u,
                    v=v,
                    jac=float(row["jaccard"]),
                    shared=int(row["shared_kmers"]),
                    z=parse_float(row.get("z_robust", "")),
                )
            )
    return edges


def percentile_rank(xs: List[float], x: float) -> float:
    if not xs:
        return 0.0
    le = 0
    for v in xs:
        if v <= x:
            le += 1
    return le / len(xs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", type=Path, required=True)
    ap.add_argument("--protein_features", type=Path, required=True)
    ap.add_argument("--hgt_candidates", type=Path, default=None)
    ap.add_argument("--component_ids", type=str, required=True, help="comma-separated, e.g. 5,8,32")
    ap.add_argument("--out_dir", type=Path, default=Path("figs"))

    ap.add_argument("--node_size_mode", type=str, default="score", choices=["score", "constant"])
    ap.add_argument("--constant_node_size", type=float, default=300.0)
    ap.add_argument("--score_size_min", type=float, default=45.0)
    ap.add_argument("--score_size_max", type=float, default=640.0)
    ap.add_argument("--score_size_gamma", type=float, default=1.4)
    ap.add_argument("--outlined_size_mult_constant", type=float, default=1.9)
    ap.add_argument("--outlined_size_mult_score", type=float, default=1.0)
    ap.add_argument("--top_score_outline", type=int, default=6)
    ap.add_argument("--z_min_highlight", type=float, default=5.0)
    ap.add_argument("--max_highlight_edges", type=int, default=60)
    ap.add_argument("--color_mode", type=str, default="highlight", choices=["highlight", "frequency"])
    ap.add_argument("--max_colored_species", type=int, default=12)

    ap.add_argument("--no_fade_non_outlined", action="store_true")
    ap.add_argument("--fade_alpha", type=float, default=0.30)
    ap.add_argument("--fade_scale", type=float, default=0.60)
    ap.add_argument("--context_edge_alpha", type=float, default=0.06)

    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--figsize", type=str, default="11.69,8.27")
    ap.add_argument("--spring_k_scale", type=float, default=2.2,
                    help="Larger values spread nodes more in spring layout.")
    ap.add_argument("--spring_iterations", type=int, default=250)
    ap.add_argument("--legend_fontsize", type=int, default=10)
    ap.add_argument("--caption_fontsize", type=int, default=11)
    ap.add_argument("--dpi", type=int, default=900)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_nodes(args.protein_features)
    scores = load_scores(args.hgt_candidates)
    have_scores = bool(scores)

    comp_ids = [int(x.strip()) for x in args.component_ids.split(",") if x.strip()]
    comp_to_nodes: Dict[int, List[str]] = defaultdict(list)
    for u, n in nodes.items():
        comp_to_nodes[n.component_id].append(u)

    fig_w, fig_h = (float(x) for x in args.figsize.split(","))

    for cid in comp_ids:
        comp_nodes = comp_to_nodes.get(cid, [])
        if not comp_nodes:
            print(f"[WARN] component {cid} not found.")
            continue
        comp_set = set(comp_nodes)
        edges = load_edges_for_component(args.edges, comp_set)

        g = nx.Graph()
        species_by_node: Dict[str, str] = {}
        score_by_node: Dict[str, float] = {}
        for u in comp_nodes:
            sp = nodes[u].species
            sc = scores.get(u, 0.0)
            g.add_node(u, species=sp, score=sc)
            species_by_node[u] = sp
            score_by_node[u] = sc
        for e in edges:
            g.add_edge(e.u, e.v, z=e.z, jac=e.jac, shared=e.shared)

        comp_nodes_now = list(g.nodes())
        sp_list = [species_by_node[u] for u in comp_nodes_now]
        sp_counts = Counter(sp_list)

        edges_with_z: List[Tuple[str, str, float]] = []
        for u, v, d in g.edges(data=True):
            z = d.get("z")
            if z is not None:
                edges_with_z.append((u, v, float(z)))
        edges_with_z.sort(key=lambda t: t[2], reverse=True)

        highlight_edges = [(u, v, z) for (u, v, z) in edges_with_z if z >= args.z_min_highlight]
        highlight_edges = highlight_edges[:args.max_highlight_edges]

        if args.color_mode == "highlight":
            sp_in_hl = Counter()
            for u, v, _ in highlight_edges:
                sp_in_hl[species_by_node[u]] += 1
                sp_in_hl[species_by_node[v]] += 1
            if have_scores:
                topk = sorted(comp_nodes_now, key=lambda u: score_by_node.get(u, 0.0), reverse=True)[: args.top_score_outline]
                for u in topk:
                    sp_in_hl[species_by_node[u]] += 1
            top_species = [sp for sp, _ in sp_in_hl.most_common(args.max_colored_species)]
        else:
            top_species = [sp for sp, _ in sp_counts.most_common(args.max_colored_species)]

        palette = [
            "#1f77b4", "#2ca02c", "#9467bd", "#17becf", "#bcbd22", "#8c564b",
            "#7f7f7f", "#e377c2", "#1b9e77", "#66a61e", "#7570b3", "#a6761d",
        ]
        sp_to_color: Dict[str, str] = {}
        for i, sp in enumerate(top_species):
            sp_to_color[sp] = palette[i % len(palette)]
        other_color = "#D0D0D0"

        node_size_map: Dict[str, float] = {}
        if args.node_size_mode == "score" and have_scores:
            comp_scores = [score_by_node.get(u, 0.0) for u in comp_nodes_now]
            for u in comp_nodes_now:
                p = percentile_rank(comp_scores, score_by_node.get(u, 0.0))
                node_size_map[u] = args.score_size_min + (args.score_size_max - args.score_size_min) * (p ** args.score_size_gamma)
        else:
            for u in comp_nodes_now:
                node_size_map[u] = args.constant_node_size

        outlined: set[str] = set()
        if have_scores:
            topk = sorted(comp_nodes_now, key=lambda u: score_by_node.get(u, 0.0), reverse=True)[: args.top_score_outline]
            outlined = set(topk)

        base_edges = list(g.edges())
        hl = [(u, v) for (u, v, _) in highlight_edges]
        hl_widths = [min(3.2, 0.9 + 0.12 * z) for (_, _, z) in highlight_edges]

        hl_nodes: set[str] = set()
        for u, v, _ in highlight_edges:
            hl_nodes.add(u)
            hl_nodes.add(v)

        n_layout = g.number_of_nodes()
        k = args.spring_k_scale / (max(1, n_layout) ** 0.5)
        pos = nx.spring_layout(g, seed=args.seed, k=k, iterations=args.spring_iterations)
        plt.figure(figsize=(fig_w, fig_h))
        nx.draw_networkx_edges(g, pos, edgelist=base_edges, width=0.6, alpha=0.10, edge_color="#808080")
        if hl:
            nx.draw_networkx_edges(g, pos, edgelist=hl, width=hl_widths, alpha=0.9, edge_color="black")

        def node_col(u: str) -> str:
            return sp_to_color.get(species_by_node[u], other_color)

        non_outline_nodes = [u for u in g.nodes() if u not in outlined]
        outline_nodes = [u for u in g.nodes() if u in outlined]

        if args.no_fade_non_outlined:
            non_outline_hi = non_outline_nodes
            non_outline_lo = []
        else:
            non_outline_hi = [u for u in non_outline_nodes if u in hl_nodes]
            non_outline_lo = [u for u in non_outline_nodes if u not in hl_nodes]

        faded_size_mult = args.fade_scale if args.node_size_mode == "constant" else 1.0

        if non_outline_lo:
            nx.draw_networkx_nodes(
                g,
                pos,
                nodelist=non_outline_lo,
                node_color=[node_col(u) for u in non_outline_lo],
                node_size=[node_size_map[u] * faded_size_mult for u in non_outline_lo],
                alpha=args.fade_alpha,
                linewidths=0.0,
            )
        if non_outline_hi:
            nx.draw_networkx_nodes(
                g,
                pos,
                nodelist=non_outline_hi,
                node_color=[node_col(u) for u in non_outline_hi],
                node_size=[node_size_map[u] for u in non_outline_hi],
                alpha=1.0,
                linewidths=0.0,
            )
        if outline_nodes:
            outline_mult = args.outlined_size_mult_score if args.node_size_mode == "score" else args.outlined_size_mult_constant
            nx.draw_networkx_nodes(
                g,
                pos,
                nodelist=outline_nodes,
                node_color=[node_col(u) for u in outline_nodes],
                node_size=[node_size_map[u] * outline_mult for u in outline_nodes],
                alpha=1.0,
                linewidths=3.2,
                edgecolors="#D00000",
            )

        n = g.number_of_nodes()
        m = g.number_of_edges()
        dens = (2.0 * m) / (n * (n - 1.0)) if n > 1 else 0.0
        component_titles = {
            5: "Broad Multi-Taxa Exchange Module",
            8: "Chained Cross-Species Connectivity",
            32: "Localized High-Confidence Bridges",
        }
        title_suffix = component_titles.get(cid, "Cross-Species Transfer Signal Map")
        plt.title(f"Component {cid}: {title_suffix}", fontsize=16)
        plt.axis("off")

        handles = []
        for sp in top_species:
            handles.append(Patch(facecolor=sp_to_color[sp], edgecolor="none", label=sp))
        handles.append(Patch(facecolor=other_color, edgecolor="none", label="Other species"))
        handles.append(Line2D([0], [0], color="#808080", lw=2, alpha=0.25, label="All edges"))
        handles.append(Line2D([0], [0], color="black", lw=2.5, alpha=0.9, label=f"High-surprise edges (z >= {args.z_min_highlight})"))
        handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                              markeredgecolor="#D00000", markeredgewidth=2.5, markersize=12,
                              label=f"Outlined proteins: top-{args.top_score_outline} by score"))
        plt.gcf().legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.085),
            ncol=4,
            frameon=True,
            fontsize=args.legend_fontsize,
        )

        caption = (
            f"Caption: Node color = species; node size mode = {args.node_size_mode}. "
            "Red outlines mark top-scoring proteins. "
            f"Black edges are high-surprise links (z >= {args.z_min_highlight}); gray edges provide context. "
            f"Stats: nodes={n}, edges={m}, density={dens:.3f}."
        )
        plt.figtext(0.5, 0.02, caption, ha="center", va="bottom", fontsize=args.caption_fontsize, wrap=True)

        out_path = args.out_dir / f"component_{cid}.png"
        plt.tight_layout(rect=[0, 0.18, 1, 0.94])
        plt.savefig(out_path, dpi=args.dpi)
        plt.close()
        print(f"[OK] wrote {out_path}")


if __name__ == "__main__":
    main()
