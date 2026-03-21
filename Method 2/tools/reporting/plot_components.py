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
import math
import textwrap
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import colors as mcolors
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


def distinct_species_palette() -> List[str]:
    # High-contrast categorical colors (colorblind-friendlier than default cycle).
    return [
        "#0072B2", "#56B4E9", "#009E73", "#1B9E77", "#66A61E", "#A6D854",
        "#F0E442", "#B3DE69", "#332288", "#393B79", "#7570B3", "#984EA3",
        "#6A3D9A", "#8C6D31", "#A6761D", "#B15928", "#D62728", "#E41A1C",
        "#17BECF", "#5AB4AC",
    ]


def rgb_dist(c1: str, c2: str) -> float:
    r1, g1, b1 = mcolors.to_rgb(c1)
    r2, g2, b2 = mcolors.to_rgb(c2)
    return math.sqrt((r1 - r2) ** 2 + (g1 - g2) ** 2 + (b1 - b2) ** 2)


def assign_species_colors(
    top_species: List[str],
    palette: List[str],
    species_adj: Dict[str, set[str]],
    species_counts: Counter,
) -> Dict[str, str]:
    # Greedy max-separation coloring on the species-interaction graph.
    ordered_species = sorted(
        top_species,
        key=lambda sp: (len(species_adj.get(sp, set())), species_counts.get(sp, 0)),
        reverse=True,
    )
    assigned: Dict[str, str] = {}
    color_use: Counter = Counter()
    used_once: set[str] = set()
    for sp in ordered_species:
        if len(used_once) < len(palette):
            candidates = [c for c in palette if c not in used_once]
        else:
            candidates = list(palette)
        best_color = candidates[0]
        best_score = (-1.0, -1.0, float("-inf"))
        for c in candidates:
            neigh = [assigned[n] for n in species_adj.get(sp, set()) if n in assigned]
            if neigh:
                min_neigh = min(rgb_dist(c, nc) for nc in neigh)
            else:
                min_neigh = 10.0
            if assigned:
                min_global = min(rgb_dist(c, ac) for ac in assigned.values())
            else:
                min_global = 10.0
            usage_bonus = -float(color_use[c])
            score = (min_neigh, min_global, usage_bonus)
            if score > best_score:
                best_score = score
                best_color = c
        assigned[sp] = best_color
        color_use[best_color] += 1
        used_once.add(best_color)
    return assigned


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
    ap.add_argument("--z_min_highlight", type=float, default=3.0)
    ap.add_argument("--max_highlight_edges", type=int, default=60)
    ap.add_argument("--color_mode", type=str, default="highlight", choices=["highlight", "frequency"])
    ap.add_argument("--max_colored_species", type=int, default=12)

    ap.add_argument("--no_fade_non_outlined", action="store_true")
    ap.add_argument("--fade_alpha", type=float, default=0.30)
    ap.add_argument("--fade_scale", type=float, default=0.60)
    ap.add_argument("--context_edge_alpha", type=float, default=0.06)

    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--figsize", type=str, default="11.69,9.8")
    ap.add_argument("--spring_k_scale", type=float, default=2.2,
                    help="Larger values spread nodes more in spring layout.")
    ap.add_argument("--layout_y_scale", type=float, default=1.28,
                    help="Multiply y coordinates after spring layout to stretch vertical separation.")
    ap.add_argument("--layout_x_scale", type=float, default=0.94,
                    help="Multiply x coordinates after spring layout to slightly compress horizontal spread.")
    ap.add_argument("--spring_iterations", type=int, default=250)
    ap.add_argument("--title_fontsize", type=int, default=20)
    ap.add_argument("--legend_fontsize", type=int, default=12)
    ap.add_argument("--caption_fontsize", type=int, default=12)
    ap.add_argument("--dpi", type=int, default=900)
    ap.add_argument(
        "--highlight_policy_label",
        type=str,
        default="",
        help="Optional caption/legend label for highlighted-edge policy (e.g., 'Top 5%% by z_robust').",
    )
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
            # Ensure every species incident to bold edges is explicitly color-coded.
            top_species = [sp for sp, _ in sp_in_hl.most_common()]
            if not top_species:
                top_species = [sp for sp, _ in sp_counts.most_common(args.max_colored_species)]
        else:
            top_species = [sp for sp, _ in sp_counts.most_common(args.max_colored_species)]

        palette = distinct_species_palette()
        top_set = set(top_species)
        species_adj: Dict[str, set[str]] = defaultdict(set)
        for u, v in g.edges():
            su = species_by_node[u]
            sv = species_by_node[v]
            if su == sv:
                continue
            if su in top_set and sv in top_set:
                species_adj[su].add(sv)
                species_adj[sv].add(su)
        sp_to_color = assign_species_colors(top_species, palette, species_adj, sp_counts)
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

        base_edges = list(g.edges())
        hl = [(u, v) for (u, v, _) in highlight_edges]
        hl_widths = [min(3.2, 0.9 + 0.12 * z) for (_, _, z) in highlight_edges]

        hl_nodes: set[str] = set()
        for u, v, _ in highlight_edges:
            hl_nodes.add(u)
            hl_nodes.add(v)

        n = g.number_of_nodes()
        m = g.number_of_edges()
        dens = (2.0 * m) / (n * (n - 1.0)) if n > 1 else 0.0
        if args.node_size_mode == "score":
            size_text = "Node size reflects the assigned HGT candidate score (larger means higher score)."
        else:
            size_text = "Node size is constant in this rendering mode."
        if args.highlight_policy_label.strip():
            highlight_policy = args.highlight_policy_label.strip()
        elif args.z_min_highlight <= -900:
            highlight_policy = f"Top-{args.max_highlight_edges} edges by z_robust"
        else:
            highlight_policy = f"Hard threshold: z >= {args.z_min_highlight:g}"
        policy_lc = highlight_policy.lower()
        if "top 5%" in policy_lc or "top5%" in policy_lc:
            legend_policy = "top5%"
        elif policy_lc.startswith("top-") and "edges" in policy_lc:
            legend_policy = f"top{args.max_highlight_edges}"
        elif "z >=" in policy_lc:
            legend_policy = "z>=" + highlight_policy.split("z >=", 1)[1].strip()
        elif "hard threshold:" in policy_lc and "z >=" in policy_lc:
            legend_policy = "z>=" + highlight_policy.split("z >=", 1)[1].strip()
        else:
            legend_policy = "policy"

        caption = (
            "Each node is a protein. Each weighted edge encodes similarity between two proteins. "
            "Node color indicates the source species of that protein. "
            f"{size_text} "
            f"Highlight policy: {highlight_policy}. "
            "Black edges are highlighted links; gray edges provide context. "
            f"Stats: nodes={n}, edges={m}, density={dens:.3f}."
        )
        wrapped_caption = textwrap.fill(caption, width=125)

        legend_item_count = len(top_species) + 4
        legend_cols = 3
        legend_rows = max(1, math.ceil(legend_item_count / legend_cols))
        extra_legend_height = max(0.0, (legend_rows - 2) * 0.7)
        caption_lines = wrapped_caption.count("\n") + 1
        extra_caption_height = max(0.55, (caption_lines - 1) * 0.22)
        fig_h_local = fig_h + extra_legend_height + extra_caption_height

        n_layout = g.number_of_nodes()
        k = args.spring_k_scale / (max(1, n_layout) ** 0.5)
        pos = nx.spring_layout(g, seed=args.seed, k=k, iterations=args.spring_iterations)
        if args.layout_y_scale != 1.0 or args.layout_x_scale != 1.0:
            pos = {
                u: (xy[0] * args.layout_x_scale, xy[1] * args.layout_y_scale)
                for u, xy in pos.items()
            }
        plt.figure(figsize=(fig_w, fig_h_local))
        nx.draw_networkx_edges(g, pos, edgelist=base_edges, width=0.6, alpha=0.10, edge_color="#808080")
        if hl:
            nx.draw_networkx_edges(g, pos, edgelist=hl, width=hl_widths, alpha=0.9, edge_color="black")

        def node_col(u: str) -> str:
            return sp_to_color.get(species_by_node[u], other_color)

        all_nodes = list(g.nodes())

        if args.no_fade_non_outlined:
            non_outline_hi = all_nodes
            non_outline_lo = []
        else:
            non_outline_hi = [u for u in all_nodes if u in hl_nodes]
            non_outline_lo = [u for u in all_nodes if u not in hl_nodes]

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
                linewidths=0.45,
                edgecolors="#111111",
            )

        component_titles = {
            5: "Broad Multi-Taxa Exchange Module",
            8: "Chained Cross-Species Connectivity",
            32: "Localized High-Confidence Bridges",
        }
        title_suffix = component_titles.get(cid, "Cross-Species Transfer Signal Map")
        plt.title(f"Component {cid}: {title_suffix}", fontsize=args.title_fontsize)
        plt.axis("off")

        handles = []
        for sp in top_species:
            handles.append(Patch(facecolor=sp_to_color[sp], edgecolor="none", label=sp))
        handles.append(Patch(facecolor=other_color, edgecolor="none", label="Other species"))
        handles.append(Line2D([0], [0], color="#808080", lw=2, alpha=0.25, label="All edges"))
        handles.append(Line2D([0], [0], color="black", lw=2.5, alpha=0.9, label=f"Highlighted edges ({legend_policy})"))
        plt.gcf().legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.085),
            ncol=legend_cols,
            frameon=True,
            fontsize=args.legend_fontsize,
        )

        plt.figtext(
            0.06,
            0.02,
            wrapped_caption,
            ha="left",
            va="bottom",
            fontsize=args.caption_fontsize,
            multialignment="left",
        )

        out_path = args.out_dir / f"component_{cid}.png"
        plt.tight_layout(rect=(0.04, 0.18, 0.96, 0.94))
        plt.savefig(out_path, dpi=args.dpi)
        plt.close()
        print(f"[OK] wrote {out_path}")


if __name__ == "__main__":
    main()
