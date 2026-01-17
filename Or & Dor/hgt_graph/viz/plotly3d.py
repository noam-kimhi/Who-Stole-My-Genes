import math
import os
import networkx as nx
import plotly.graph_objs as go
from typing import Dict, List, Set, Tuple, Mapping, Any, Union
from ..constants import (EDGE_WEIGHT_KEY, IDENTITY_KEY, COVERAGE_KEY, ALIGNED_LENGTH_KEY, TAX_DIST_KEY,
                         ORGANISM_KEY, PROTEIN_NAME_KEY, SEQ_LENGTH_KEY, PLOTLY_MIN_W, PLOTLY_MAX_W,
                         SUS_EDGE_COLOR, REG_EDGE_COLOR, SUS_NODE_BORDER_COLOR, REG_NODE_BORDER_COLOR,
                         ORG_NAME_LBL, PathLike, NORMAL_NODES, HIGHLIGHT_NODES, X_VALS_KEY, Y_VALS_KEY,
                         Z_VALS_KEY, HOVER_TEXT_KEY, TEXT_SIZE_KEY, TEXT_COLOR_KEY, LABEL_KEY)
from .colors import stable_color_from_string, rescale_weight


def build_edge_traces(G: nx.Graph, pos: Dict[str, Tuple[float, float, float]],
                      suspicious_edges: Set[Tuple[str, str]], show_edge_hover: bool = True) -> List[go.Scatter3d]:
    """
    Build edge traces for 3D Plotly visualization.
    :param G: The input NetworkX graph.
    :param pos: The 3D positions of nodes.
    :param suspicious_edges: A set of edges considered suspicious for HGT.
    :param show_edge_hover: Whether to show hover info on edges.
    :return: A list of Plotly Scatter3d traces for edges.
    """
    traces = []

    node_id_to_org_name = {nid: attrs.get(ORGANISM_KEY, '') for nid, attrs in G.nodes(data=True)}

    weights = [float(edata.get(EDGE_WEIGHT_KEY, 0.0)) for _, _, edata in G.edges(data=True)]
    if not weights:
        return []  # no edge traces
    min_w, max_w = min(weights), max(weights)

    for u, v, edata in G.edges(data=True):
        x0, y0, z0 = pos[u]
        x1, y1, z1 = pos[v]

        w = float(edata.get(EDGE_WEIGHT_KEY, 0.0))
        width = rescale_weight(w, min_w, max_w, PLOTLY_MIN_W, PLOTLY_MAX_W)
        identity = edata.get(IDENTITY_KEY, 0.0)
        coverage = edata.get(COVERAGE_KEY, 0.0)
        aligned_len = edata.get(ALIGNED_LENGTH_KEY, 0)
        tax_d = edata.get(TAX_DIST_KEY, -1)

        key = (u, v) if u < v else (v, u)
        is_suspicious = key in suspicious_edges

        color = SUS_EDGE_COLOR if is_suspicious else REG_EDGE_COLOR

        if is_suspicious:
            width = min(width * 1.35, 14.0)  # Enhance width for suspicious edges

        hover_text = None
        if show_edge_hover:
            hover_text = (
                f'{u} ({node_id_to_org_name.get(u)}) <> {v} ({node_id_to_org_name.get(v)})<br>'
                f'weight={w:.4f}<br>'
                f'identity={identity:.4f}<br>'
                f'coverage={coverage:.4f}<br>'
                f'tax_distance={tax_d}<br>'
                f'aligned_len={aligned_len}'
            )

        traces.append(
            go.Scatter3d(
                x=[x0, x1],
                y=[y0, y1],
                z=[z0, z1],
                mode='lines',
                line=dict(
                    width=width,
                    color=color
                ),
                hoverinfo='text' if show_edge_hover else 'none',
                hovertext=hover_text,
                showlegend=False
            )
        )

    return traces


def make_node_traces(G: nx.Graph, pos: Dict[str, Tuple[float, float, float]], highlight_nodes: Set[str],
                     hgt_scores: Dict[str, float], node_label: str) -> List[go.Scatter3d]:
    """
    Create node traces for 3D Plotly visualization.
    :param G: The input NetworkX graph.
    :param pos: The 3D positions of nodes.
    :param highlight_nodes: A set of node IDs to highlight.
    :param hgt_scores: A dictionary of HGT suspicion scores for nodes.
    :param node_label: Choice of node label: "name", "id", or "none".
    :return: A list of Plotly Scatter3d traces for nodes.
    """
    def node_display(attrs: Mapping[str, Any]) -> str:
        org: str = str(attrs.get(ORGANISM_KEY, ''))
        protein_name: str = str(attrs.get(PROTEIN_NAME_KEY, ''))
        return org if (node_label == ORG_NAME_LBL and org) else protein_name

    # containers for the two groups
    groups: Dict[str, Dict[str, List[Union[float, str]]]] = {
        NORMAL_NODES: {
            X_VALS_KEY: [], Y_VALS_KEY: [], Z_VALS_KEY: [], HOVER_TEXT_KEY: [],
            TEXT_SIZE_KEY: [], TEXT_COLOR_KEY: [], LABEL_KEY: []
        },
        HIGHLIGHT_NODES: {
            X_VALS_KEY: [], Y_VALS_KEY: [], Z_VALS_KEY: [], HOVER_TEXT_KEY: [],
            TEXT_SIZE_KEY: [], TEXT_COLOR_KEY: [], LABEL_KEY: []
        },
    }

    for nid, attrs in G.nodes(data=True):
        x, y, z = pos[nid]
        organism = str(attrs.get(ORGANISM_KEY, 'unknown'))
        name = str(attrs.get(PROTEIN_NAME_KEY, ''))
        seq_len = int(attrs.get(SEQ_LENGTH_KEY, 0))
        phylum = attrs.get('phylum')
        genus = attrs.get('genus')
        score = float(hgt_scores.get(nid, 0.0))

        disp = node_display(attrs)

        # Determine node size based on sequence length
        base_size = 8 + 6 * math.log10(seq_len + 1) if seq_len > 0 else 10
        # Increase size for highlighted nodes
        size = base_size * 1.25 if nid in highlight_nodes else base_size

        hover = (
            f'<b>{disp}</b><br>'
            f'<b>ID:</b> {nid}<br>'
            f'<b>Name:</b> {name}<br>'
            f'<b>Organism:</b> {organism}<br>'
            f'<b>Genus:</b> {genus}<br>'
            f'<b>Phylum:</b> {phylum}<br>'
            f'<b>Length:</b> {seq_len}<br>'
            f'<b>HGT score:</b> {score:.4f}'
        )

        grp = HIGHLIGHT_NODES if nid in highlight_nodes else NORMAL_NODES
        groups[grp][X_VALS_KEY].append(x)
        groups[grp][Y_VALS_KEY].append(y)
        groups[grp][Z_VALS_KEY].append(z)
        groups[grp][HOVER_TEXT_KEY].append(hover)
        groups[grp][TEXT_SIZE_KEY].append(size)
        groups[grp][TEXT_COLOR_KEY].append(stable_color_from_string(organism))
        groups[grp][LABEL_KEY].append(disp if node_label != 'none' else '')

    # normal nodes
    traces = [go.Scatter3d(
        x=groups[NORMAL_NODES][X_VALS_KEY],
        y=groups[NORMAL_NODES][Y_VALS_KEY],
        z=groups[NORMAL_NODES][Z_VALS_KEY],
        mode='markers+text' if node_label != 'none' else 'markers',
        text=groups[NORMAL_NODES][LABEL_KEY] if node_label != 'none' else None,
        textposition='top center',
        hoverinfo='text',
        hovertext=groups[NORMAL_NODES][HOVER_TEXT_KEY],
        marker=dict(
            size=groups[NORMAL_NODES][TEXT_SIZE_KEY],
            color=groups[NORMAL_NODES][TEXT_COLOR_KEY],
            opacity=0.85,
            line=dict(width=1, color=REG_NODE_BORDER_COLOR),
        ),
        name='nodes'
    )]

    # highlighted nodes
    if groups[HIGHLIGHT_NODES][X_VALS_KEY]:
        traces.append(
            go.Scatter3d(
                x=groups[HIGHLIGHT_NODES][X_VALS_KEY],
                y=groups[HIGHLIGHT_NODES][Y_VALS_KEY],
                z=groups[HIGHLIGHT_NODES][Z_VALS_KEY],
                mode='markers+text' if node_label != 'none' else 'markers',
                text=groups[HIGHLIGHT_NODES][LABEL_KEY] if node_label != 'none' else None,
                textposition='top center',
                hoverinfo='text',
                hovertext=groups[HIGHLIGHT_NODES][HOVER_TEXT_KEY],
                marker=dict(
                    size=groups[HIGHLIGHT_NODES][TEXT_SIZE_KEY],
                    color=groups[HIGHLIGHT_NODES][TEXT_COLOR_KEY],
                    opacity=0.95,
                    line=dict(width=5, color=SUS_NODE_BORDER_COLOR),
                ),
                name='HGT candidates'
            )
        )

    return traces


def export_plotly_3d(G: nx.Graph, out_html: PathLike, sus_edges: Set[Tuple[str, str]], hgt_scores: Dict[str, float],
                     highlight_nodes: Set[str], node_label: str = ORG_NAME_LBL, show_edge_hover: bool = True) -> None:
    """
    Export the given NetworkX graph to a rotatable 3D Plotly HTML visualization.
    :param G: The input NetworkX graph.
    :param out_html: The output HTML file path.
    :param sus_edges: A set of edges considered suspicious for HGT.
    :param hgt_scores: A dictionary of HGT suspicion scores for nodes.
    :param highlight_nodes: A set of node IDs to highlight.
    :param node_label: Choice of node label: "name" or "id".
    :param show_edge_hover: Whether to show hover info on edges.
    """
    pos: Dict[str, Tuple[float, float, float]] = nx.spring_layout(G, dim=3, seed=42, weight=EDGE_WEIGHT_KEY)

    edge_traces = build_edge_traces(G, pos, sus_edges, show_edge_hover)

    node_traces = make_node_traces(G, pos, highlight_nodes, hgt_scores, node_label)

    fig = go.Figure(data=edge_traces + node_traces)
    fig.update_layout(
        title="Protein Similarity Graph (3D)",
        showlegend=False,
        scene=dict(
            xaxis=dict(showbackground=False, visible=False),
            yaxis=dict(showbackground=False, visible=False),
            zaxis=dict(showbackground=False, visible=False),
        ),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"Saved 3D interactive graph to: {os.path.abspath(out_html)}")
