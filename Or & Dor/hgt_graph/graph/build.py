import networkx as nx
from typing import Dict, List, Tuple, Union
from ..constants import IDENTITY_KEY


def build_nx_graph(node_attrs: Dict[str, Dict[str, Union[str, int]]],
                   edges: List[Tuple[str, str, Dict[str, float]]]) -> nx.Graph:
    """
    Build a NetworkX graph from node attributes and edges.
    :param node_attrs: A dictionary of node attributes.
    :param edges: A list of edges with weights.
    :return: A NetworkX graph.
    """
    G = nx.Graph()
    for node_id, attrs in node_attrs.items():
        G.add_node(node_id, **attrs)
    for u, v, eattrs in edges:
        G.add_edge(u, v, **eattrs)
    return G


def build_identity_map_from_edges(edges: List[Tuple[str, str, Dict[str, float]]]) -> Dict[Tuple[str, str], float]:
    """
    Build a pairwise identity map from graph edges.
    :param edges: The edges of the graph with attributes.
    :return: A dictionary mapping (node1, node2) to their identity score.
    """
    ident_map: Dict[Tuple[str, str], float] = {}
    for u, v, attrs in edges:
        if IDENTITY_KEY in attrs:
            ident_map[(u, v)] = float(attrs[IDENTITY_KEY])
    return ident_map
