import math
import networkx as nx
from typing import Dict, Set, Tuple, List, Iterable
from ..constants import HGT_TAX_DISTANCE_MIN, EDGE_WEIGHT_KEY, TAX_DIST_KEY, HGT_MIN_WEIGHT, TOP_HGT_N, MIN_SUS_EDGES, TOP_N_EDGES

def compute_hgt_scores(G: nx.Graph, min_tax_distance: int = HGT_TAX_DISTANCE_MIN,
                       min_weight: float = HGT_MIN_WEIGHT) -> Tuple[Dict[str, float],
                                                                    Dict[str, int],
                                                                    Set[Tuple[str, str]]]:
    """
    Compute HGT suspicion scores for each node based on edges that meet criteria.
    :param G: The input NetworkX graph.
    :param min_tax_distance: The minimum taxonomic distance to consider an edge suspicious.
    :param min_weight: The minimum edge weight to consider an edge suspicious.
    :return: A tuple containing:
             - A dictionary of node scores.
             - A dictionary of suspicious edge counts per node.
             - A set of suspicious edges.
    """
    scores = {n: 0.0 for n in G.nodes()}
    sus_degree = {n: 0 for n in G.nodes()}
    suspicious_edges: Set[Tuple[str, str]] = set()

    for u, v, edata in G.edges(data=True):
        w = float(edata.get(EDGE_WEIGHT_KEY, 0.0))
        d = int(edata.get(TAX_DIST_KEY, -1))

        if w >= min_weight and d >= min_tax_distance:
            key = (u, v) if u < v else (v, u)
            suspicious_edges.add(key)

            scores[u] += w
            scores[v] += w
            sus_degree[u] += 1
            sus_degree[v] += 1

    return scores, sus_degree, suspicious_edges


def percentile(values: List[float], p: float = 95.0) -> float:
    """
    Compute the p-th percentile of a list of values.
    :param values: A list of numeric values.
    :param p: The desired percentile (0-100).
    :return: The p-th percentile value.
    """
    if not values:
        return 0.0
    xs = sorted(values)
    k = (len(xs) - 1) * (p / 100.0)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return xs[int(k)]
    return xs[f] + (k - f) * (xs[c] - xs[f])


def top_hgt_candidates(nodes: Iterable[str], scores: Dict[str, float], sus_degree: Dict[str, int],
                       score_thresh: float, top_n: int = TOP_HGT_N) -> List[Tuple[str, float]]:
    """
    Get the top N HGT candidate nodes based on scores and suspicious degree.
    :param nodes: An iterable of node IDs.
    :param scores: A dictionary of node scores.
    :param sus_degree: A dictionary of suspicious edge counts per node.
    :param score_thresh: The minimum score threshold to consider a node.
    :param top_n: The number of top candidates to return.
    :return: A list of tuples containing node IDs and their scores.
    """
    eligible = [
        (node, scores[node])
        for node in nodes
        if scores[node] >= score_thresh and sus_degree[node] >= MIN_SUS_EDGES
    ]

    eligible_sorted = sorted(eligible, key=lambda x: x[1], reverse=True)
    return eligible_sorted[:top_n]


def top_hgt_edges(suspicious_edges: Set[Tuple[str, str]], scores: Dict[str, float],
                  top_n: int = TOP_N_EDGES) -> Set[Tuple[str, str]]:
    """
    Get the top N suspicious edges based on the sum of their nodes' scores.
    :param suspicious_edges: A set of suspicious edges (tuples of node IDs).
    :param scores: A dictionary of node scores.
    :param top_n: The number of top edges to return.
    :return: A set of tuples containing edges and their combined scores.
    """
    def edge_score(edge: Tuple[str, str]) -> float:
        u, v = edge
        return scores[u] + scores[v]

    sorted_edges = sorted(suspicious_edges, key=edge_score, reverse=True)
    top_edges = {edge for edge in sorted_edges[:top_n]}
    return top_edges
