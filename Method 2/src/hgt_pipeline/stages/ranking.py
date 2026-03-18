import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .component_features import ComponentFeatures
from .node_features import NodeFeatures
from .pair_stats import mad, median


def _minmax_scale(values: List[float]) -> List[float]:
    """Min-max scale to [0,1]. If constant, returns all zeros."""
    if not values:
        return []
    vmin = min(values)
    vmax = max(values)
    if vmax <= vmin:
        return [0.0] * len(values)
    return [(v - vmin) / (vmax - vmin) for v in values]


def score_hgt_likeness(
    node_feats: List[NodeFeatures],
    comp_feats: List[ComponentFeatures],
    min_component_size: int = 10,
    use_betweenness: bool = True,
    H_min: float = 0.75,
    f_min: Optional[float] = None,
    density_eta: float = 2.0,
    cluster_eta: float = 1.0,
    gamma_highz: float = 1.0,
    lambda_betweenness: float = 1.0,
    z_cap: float = 10.0,
    scale_floor_frac: float = 1e-3,
    context_floor: float = 1e-6,
    min_deg: int = 4,
    min_n_species: int = 2,
    require_positive_z: bool = True,
    entropy_node_floor: float = 0.15,
    require_any_anomaly: bool = True,
    z_min: float = 2.0,
    require_multispecies_or_bridge: bool = True,
    min_species_for_multispecies: int = 3,
    bw_min: float = 1e-12,
) -> List[Tuple[str, float]]:
    comp_by_id: Dict[int, ComponentFeatures] = {c.component_id: c for c in comp_feats}

    by_c: Dict[int, List[NodeFeatures]] = defaultdict(list)
    for nf in node_feats:
        by_c[nf.component_id].append(nf)

    def robust_center_scale(xs: List[float]) -> Tuple[float, float]:
        m = median(xs)
        d = mad(xs, m)
        s = 1.4826 * d
        floor = scale_floor_frac * (abs(m) + 1.0)
        return m, max(s, floor)

    comp_baselines: Dict[int, Dict[str, Tuple[float, float]]] = {}
    for cid, lst in by_c.items():
        z_vals = [(nf.top5_mean_z or 0.0) for nf in lst]
        f_vals = [nf.max_species_fraction for nf in lst]
        d_vals = [float(nf.deg_xsp) for nf in lst]
        if use_betweenness:
            b_vals = [(nf.betweenness or 0.0) for nf in lst]
        else:
            b_vals = [0.0 for _ in lst]

        comp_baselines[cid] = {
            "z": robust_center_scale(z_vals),
            "f": robust_center_scale(f_vals),
            "d": robust_center_scale(d_vals),
            "b": robust_center_scale(b_vals),
        }

    def clamp_pos(z: float) -> float:
        if z <= 0.0:
            return 0.0
        return z_cap if z > z_cap else z

    ranked: List[Tuple[str, float]] = []

    w_z = 1.0
    w_f = 2.0
    w_d = 0.5
    w_b = 1.0

    for nf in node_feats:
        cid = nf.component_id
        comp = comp_by_id.get(cid)
        if comp is None:
            ranked.append((nf.u, 0.0))
            continue

        comp_size = comp.size
        if comp_size < min_component_size:
            ranked.append((nf.u, 0.0))
            continue

        Hc = comp.species_entropy_norm
        if Hc < H_min:
            ranked.append((nf.u, 0.0))
            continue

        if nf.deg_xsp < min_deg:
            ranked.append((nf.u, 0.0))
            continue

        if nf.n_species < min_n_species:
            ranked.append((nf.u, 0.0))
            continue

        if f_min is not None and nf.max_species_fraction < f_min:
            ranked.append((nf.u, 0.0))
            continue

        z_u_raw = nf.top5_mean_z
        if require_positive_z:
            if z_u_raw is None or z_u_raw <= 0.0:
                ranked.append((nf.u, 0.0))
                continue

        if z_u_raw is None or z_u_raw < z_min:
            ranked.append((nf.u, 0.0))
            continue

        if require_multispecies_or_bridge:
            bw = (nf.betweenness or 0.0) if use_betweenness else 0.0
            if nf.n_species < min_species_for_multispecies and bw <= bw_min:
                ranked.append((nf.u, 0.0))
                continue

        if comp_size >= 2:
            dens = (2.0 * comp.num_edges) / (comp_size * (comp_size - 1.0))
        else:
            dens = 0.0
        dens = min(1.0, max(0.0, dens))
        dens_pen = max(1e-6, (1.0 - dens) ** density_eta)

        highz = comp.high_z_frac if comp.high_z_frac is not None else 0.0
        context = max(highz ** gamma_highz, context_floor)

        size_gate = math.log1p(comp_size)
        deg_gate = math.log1p(nf.deg_xsp)
        if deg_gate <= 0.0:
            ranked.append((nf.u, 0.0))
            continue

        cc = nf.clustering_coeff if nf.clustering_coeff is not None else 1.0
        cc = min(1.0, max(0.0, cc))
        bridge_local = max(1e-6, (1.0 - cc) ** cluster_eta)

        ent_u = nf.species_entropy_norm
        if ent_u < entropy_node_floor:
            ent_pen = 0.1 + 0.9 * (ent_u / max(entropy_node_floor, 1e-9))
        else:
            ent_pen = 1.0

        z_u = nf.top5_mean_z or 0.0
        f_u = nf.max_species_fraction
        d_u = float(nf.deg_xsp)
        b_u = (nf.betweenness or 0.0) if use_betweenness else 0.0

        (mz, sz) = comp_baselines[cid]["z"]
        (mf, sf) = comp_baselines[cid]["f"]
        (md, sd) = comp_baselines[cid]["d"]
        (mb, sb) = comp_baselines[cid]["b"]

        Zz = clamp_pos((z_u - mz) / sz)
        Zf = clamp_pos((f_u - mf) / sf)
        Zd = clamp_pos((d_u - md) / sd)
        Zb = clamp_pos((b_u - mb) / sb) if use_betweenness else 0.0

        if require_any_anomaly and (Zz == 0.0 and Zf == 0.0 and Zb == 0.0):
            ranked.append((nf.u, 0.0))
            continue

        logS = 0.0
        logS += math.log(size_gate)
        logS += math.log(context)
        logS += math.log(deg_gate)
        logS += math.log(dens_pen)
        logS += math.log(bridge_local)
        logS += math.log(ent_pen)
        logS += w_z * math.log1p(Zz)
        logS += w_f * math.log1p(Zf)
        logS += w_d * math.log1p(0.25 * Zd)
        if use_betweenness:
            logS += w_b * math.log1p(lambda_betweenness * Zb)

        ranked.append((nf.u, math.exp(logS)))

    ranked.sort(key=lambda t: t[1], reverse=True)
    return ranked


def write_hgt_candidates(
    path: Path,
    node_feats: List[NodeFeatures],
    scores: Dict[str, float],
    top_n: int = 200,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(node_feats, key=lambda nf: scores.get(nf.u, 0.0), reverse=True)[:top_n]

    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "u",
                "species",
                "component_id",
                "score",
                "deg_xsp",
                "n_species",
                "max_z",
                "top5_mean_z",
                "max_species_fraction",
                "species_entropy_norm",
                "betweenness",
                "clustering_coeff",
            ]
        )

        for nf in rows:
            score = scores.get(nf.u, 0.0)
            w.writerow(
                [
                    nf.u,
                    nf.species,
                    nf.component_id,
                    f"{score:.6f}",
                    nf.deg_xsp,
                    nf.n_species,
                    "" if nf.max_z is None else f"{nf.max_z:.6f}",
                    "" if nf.top5_mean_z is None else f"{nf.top5_mean_z:.6f}",
                    f"{nf.max_species_fraction:.6f}",
                    f"{nf.species_entropy_norm:.6f}",
                    "" if nf.betweenness is None else f"{nf.betweenness:.6f}",
                    "" if nf.clustering_coeff is None else f"{nf.clustering_coeff:.6f}",
                ]
            )


def write_all_scores(
    path: Path,
    node_feats: List[NodeFeatures],
    scores: Dict[str, float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["u", "species", "component_id", "score"])
        for nf in node_feats:
            s = scores.get(nf.u, 0.0)
            w.writerow([nf.u, nf.species, nf.component_id, f"{s:.6f}"])
