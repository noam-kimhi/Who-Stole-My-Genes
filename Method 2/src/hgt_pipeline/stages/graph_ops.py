from typing import Dict, List, Tuple

from .edge_io import Edge
from .pair_stats import EdgeFeatures

try:
    import networkx as nx
except ImportError:  # pragma: no cover
    nx = None


def build_graph(edges: List[Edge], use_networkx: bool = True):
    """
    Construct the undirected graph object.

    If use_networkx=True, returns a networkx.Graph with:
      node attrs: species
      edge attrs: shared_kmers, jaccard
    """
    if use_networkx and nx is None:
        raise RuntimeError("networkx not installed; set use_networkx=False or install it.")

    if use_networkx:
        G = nx.Graph()
        for e in edges:
            if e.u not in G:
                G.add_node(e.u, species=e.species_u)
            if e.v not in G:
                G.add_node(e.v, species=e.species_v)
            G.add_edge(e.u, e.v, shared_kmers=e.shared_kmers, jaccard=e.jaccard)
        return G

    raise NotImplementedError


def attach_edge_z_to_graph(G, edge_feats: List[EdgeFeatures]) -> None:
    """
    Attach robust z-score and pair-baseline info to graph edges.
    """
    missing = 0
    for ef in edge_feats:
        u, v = ef.u, ef.v
        if not G.has_edge(u, v):
            missing += 1
            continue

        G[u][v]["z_robust"] = ef.z_robust
        G[u][v]["pair_n"] = ef.pair_n
        G[u][v]["pair_median"] = ef.pair_median
        G[u][v]["pair_mad"] = ef.pair_mad

    if missing > 0:
        raise KeyError(
            f"attach_edge_z_to_graph: {missing} edges from edge_feats "
            f"were not found in the graph"
        )


def compute_components(G) -> Tuple[Dict[str, int], List[List[str]]]:
    """
    Compute connected components of an undirected graph.
    """
    if nx is None:
        raise RuntimeError("networkx is required for compute_components()")

    comps = [list(c) for c in nx.connected_components(G)]
    comps.sort(key=len, reverse=True)

    component_id_of: Dict[str, int] = {}
    for cid, nodes in enumerate(comps):
        for u in nodes:
            component_id_of[u] = cid

    return component_id_of, comps
