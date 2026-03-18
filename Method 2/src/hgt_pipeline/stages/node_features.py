import csv
import math
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

try:
    import networkx as nx
except ImportError:  # pragma: no cover
    nx = None


@dataclass
class NodeFeatures:
    """Computed node-level features."""

    u: str
    species: str
    component_id: int
    deg_xsp: int
    n_species: int
    top1_jac: float
    top5_mean_jac: float
    top1_shared: int
    top5_mean_shared: float
    max_z: Optional[float]
    top5_mean_z: Optional[float]
    max_species_fraction: float
    species_entropy: float
    species_entropy_norm: float
    n_eff_species: float
    participation_coeff: float
    clustering_coeff: Optional[float]
    betweenness: Optional[float]


def topk_mean(xs: Sequence[float], k: int) -> float:
    if not xs:
        return 0.0
    if k <= 0:
        raise ValueError("k must be positive")
    ys = sorted(xs, reverse=True)
    ys = ys[:k]
    return sum(ys) / len(ys)


def entropy_from_counts(counts: Dict[str, int]) -> float:
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c <= 0:
            continue
        p = c / total
        h -= p * math.log(p)
    return h


def participation_coefficient(counts: Dict[str, int], deg: int) -> float:
    if deg <= 0:
        return 0.0
    s = 0.0
    for c in counts.values():
        frac = c / deg
        s += frac * frac
    return 1.0 - s


def compute_node_features(
    G,
    component_id_of: Dict[str, int],
    z_attr: str = "z_robust",
    compute_clustering: bool = True,
    compute_betweenness: bool = True,
    betweenness_mode: str = "high_z",
    z_threshold_for_high_subgraph: float = 3.0,
) -> List[NodeFeatures]:
    if nx is None:
        raise RuntimeError("networkx is required for compute_node_features()")
    if betweenness_mode not in {"full", "per_component", "high_z"}:
        raise ValueError("betweenness_mode must be one of: full, per_component, high_z")

    clustering: Optional[Dict[str, float]] = None
    if compute_clustering:
        clustering = nx.clustering(G)

    betweenness: Optional[Dict[str, float]] = None
    if compute_betweenness:
        if betweenness_mode == "full":
            betweenness = nx.betweenness_centrality(G, normalized=True)
        elif betweenness_mode == "per_component":
            betweenness = {}
            for nodes in nx.connected_components(G):
                if len(nodes) <= 2:
                    for u in nodes:
                        betweenness[u] = 0.0
                    continue
                H = G.subgraph(nodes)
                b = nx.betweenness_centrality(H, normalized=True)
                betweenness.update(b)
        else:
            print(f"[BETWEENNESS] mode=high_z z_attr={z_attr} z_threshold={z_threshold_for_high_subgraph}")
            H = nx.Graph()
            kept_edges = 0
            t0 = time.perf_counter()
            last_log = t0
            scanned = 0
            total_edges = G.number_of_edges()
            for u, v, data in G.edges(data=True):
                scanned += 1
                z = data.get(z_attr)
                if z is None or z < z_threshold_for_high_subgraph:
                    now = time.perf_counter()
                    if now - last_log >= 2.0:
                        pct = 100.0 * scanned / max(1, total_edges)
                        print(f"[BETWEENNESS] scanning edges: {scanned}/{total_edges} ({pct:.1f}%) kept_edges={kept_edges}")
                        last_log = now
                    continue

                if u not in H:
                    H.add_node(u, **G.nodes[u])
                if v not in H:
                    H.add_node(v, **G.nodes[v])
                H.add_edge(u, v)
                kept_edges += 1

                now = time.perf_counter()
                if now - last_log >= 2.0:
                    pct = 100.0 * scanned / max(1, total_edges)
                    print(f"[BETWEENNESS] scanning edges: {scanned}/{total_edges} ({pct:.1f}%) kept_edges={kept_edges}")
                    last_log = now

            dt_scan = time.perf_counter() - t0
            print(f"[BETWEENNESS] high-z subgraph built in {dt_scan:.2f}s: H_nodes={H.number_of_nodes()} H_edges={H.number_of_edges()}")

            if kept_edges == 0:
                betweenness = {u: 0.0 for u in G.nodes()}
            else:
                t1 = time.perf_counter()
                comps = [len(c) for c in nx.connected_components(H)]
                comps.sort(reverse=True)
                largest = comps[0] if comps else 0
                ncomp = len(comps)
                med = comps[len(comps) // 2] if comps else 0
                print(f"[BETWEENNESS] H components: n={ncomp} largest={largest} median={med}")

                APPROX_IF_LARGEST_GE = 800
                APPROX_K = 200
                if largest >= APPROX_IF_LARGEST_GE:
                    print(f"[BETWEENNESS] H largest component={largest} -> using APPROX betweenness (k={APPROX_K})")
                    b = nx.betweenness_centrality(H, k=APPROX_K, seed=42, normalized=True)
                else:
                    print(f"[BETWEENNESS] computing EXACT betweenness on H ...")
                    b = nx.betweenness_centrality(H, normalized=True)

                dt_bw = time.perf_counter() - t1
                print(f"[BETWEENNESS] betweenness done in {dt_bw:.2f}s")
                betweenness = {u: b.get(u, 0.0) for u in G.nodes()}

    feats: List[NodeFeatures] = []
    for u in G.nodes():
        sp_u = G.nodes[u].get("species", "")
        neigh_species = Counter()
        jac_list: List[float] = []
        shared_list: List[float] = []
        z_list: List[float] = []

        for v, edata in G[u].items():
            sp_v = G.nodes[v].get("species", "")
            neigh_species[sp_v] += 1
            j = edata.get("jaccard")
            if j is not None:
                jac_list.append(float(j))
            sh = edata.get("shared_kmers")
            if sh is not None:
                shared_list.append(float(sh))
            z = edata.get(z_attr)
            if z is not None:
                z_list.append(float(z))

        deg = len(G[u])
        n_species = len(neigh_species)
        top1_jac = max(jac_list) if jac_list else 0.0
        top5_mean_jac = topk_mean(jac_list, 5)
        top1_shared = int(max(shared_list)) if shared_list else 0
        top5_mean_shared = topk_mean(shared_list, 5)
        max_z = max(z_list) if z_list else None
        top5_mean_z = topk_mean(z_list, 5) if z_list else None

        if deg > 0:
            max_species_fraction = max(neigh_species.values()) / deg
        else:
            max_species_fraction = 0.0

        H = entropy_from_counts(neigh_species)
        if n_species >= 2:
            H_norm = H / math.log(n_species)
        else:
            H_norm = 0.0

        n_eff = math.exp(H)
        part = participation_coefficient(neigh_species, deg)
        cc = clustering.get(u) if clustering is not None else None
        bw = betweenness.get(u) if betweenness is not None else None

        feats.append(
            NodeFeatures(
                u=u,
                species=sp_u,
                component_id=component_id_of.get(u, -1),
                deg_xsp=deg,
                n_species=n_species,
                top1_jac=top1_jac,
                top5_mean_jac=top5_mean_jac,
                top1_shared=top1_shared,
                top5_mean_shared=top5_mean_shared,
                max_z=max_z,
                top5_mean_z=top5_mean_z,
                max_species_fraction=max_species_fraction,
                species_entropy=H,
                species_entropy_norm=H_norm,
                n_eff_species=n_eff,
                participation_coeff=part,
                clustering_coeff=cc,
                betweenness=bw,
            )
        )

    return feats


def write_node_features(path: Path, feats: List[NodeFeatures]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "u",
                "species",
                "component_id",
                "deg_xsp",
                "n_species",
                "top1_jac",
                "top5_mean_jac",
                "top1_shared",
                "top5_mean_shared",
                "max_z",
                "top5_mean_z",
                "max_species_fraction",
                "species_entropy",
                "species_entropy_norm",
                "n_eff_species",
                "participation_coeff",
                "clustering_coeff",
                "betweenness",
            ]
        )
        for x in feats:
            w.writerow(
                [
                    x.u,
                    x.species,
                    x.component_id,
                    x.deg_xsp,
                    x.n_species,
                    f"{x.top1_jac:.6f}",
                    f"{x.top5_mean_jac:.6f}",
                    x.top1_shared,
                    f"{x.top5_mean_shared:.6f}",
                    "" if x.max_z is None else f"{x.max_z:.6f}",
                    "" if x.top5_mean_z is None else f"{x.top5_mean_z:.6f}",
                    f"{x.max_species_fraction:.6f}",
                    f"{x.species_entropy:.6f}",
                    f"{x.species_entropy_norm:.6f}",
                    f"{x.n_eff_species:.6f}",
                    f"{x.participation_coeff:.6f}",
                    "" if x.clustering_coeff is None else f"{x.clustering_coeff:.6f}",
                    "" if x.betweenness is None else f"{x.betweenness:.6f}",
                ]
            )
