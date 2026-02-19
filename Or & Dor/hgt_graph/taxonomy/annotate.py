import networkx as nx
from typing import Dict, Mapping, Any, List
from ..constants import ORGANISM_KEY, TAX_DIST_KEY, TAX_RANKS, TAX_RANKS_CLOSE_TO_BROAD, SPEC_NAME_KEY
from .normalize import normalize_species_name


def annotate_graph_with_taxonomy(G: nx.Graph, taxonomy_cache: Dict[str, Mapping[str, Any]]) -> None:
    """
    Annotate the graph nodes with taxonomic lineage information from the taxonomy cache.
    :param G: The input NetworkX graph.
    :param taxonomy_cache: A dictionary mapping normalized species names to their taxonomic lineage.
    """
    for node_id, attrs in G.nodes(data=True):
        org = str(attrs.get(ORGANISM_KEY, ''))
        sp = normalize_species_name(org)
        lineage = taxonomy_cache.get(sp, {})
        G.nodes[node_id][SPEC_NAME_KEY] = sp

        for rank in TAX_RANKS:
            if rank in lineage:
                G.nodes[node_id][rank] = lineage[rank]
            else:
                G.nodes[node_id][rank] = None  # Set missing ranks to None for consistency


def taxonomic_distance(u_attrs: Mapping[str, Any], v_attrs: Mapping[str, Any]) -> int:
    """
    Compute taxonomic distance between two nodes based on their taxonomic attributes.
    :param u_attrs: The attributes of node u.
    :param v_attrs: The attributes of node v.
    :return: The taxonomic distance as an integer.
    """
    ranks: List[str] = TAX_RANKS_CLOSE_TO_BROAD
    for d, rank in enumerate(ranks):
        u_val = u_attrs.get(rank)
        v_val = v_attrs.get(rank)
        if u_val is not None and v_val is not None and u_val == v_val:
            return d
    return len(ranks)


def has_any_taxonomy(attrs: Mapping[str, Any]) -> bool:
    """
    Check if the given attributes contain any taxonomic information.
    :param attrs: The node attributes to check.
    :return: True if any taxonomic rank is present, False otherwise.
    """
    return any(attrs.get(rank) is not None for rank in TAX_RANKS)


def annotate_edges_with_tax_distance(G: nx.Graph) -> None:
    """
    Annotate graph edges with taxonomic distance between connected nodes.
    :param G: The input NetworkX graph.
    """
    for u, v in G.edges():
        # if either node lacks taxonomy, mark unknown
        if not has_any_taxonomy(G.nodes[u]) or not has_any_taxonomy(G.nodes[v]):
            G.edges[u, v][TAX_DIST_KEY] = -1
            continue

        d: int = taxonomic_distance(G.nodes[u], G.nodes[v])
        G.edges[u, v][TAX_DIST_KEY] = d
