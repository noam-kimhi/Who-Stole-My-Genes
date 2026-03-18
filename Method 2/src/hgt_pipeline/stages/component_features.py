import csv
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from .node_features import entropy_from_counts, topk_mean


@dataclass
class ComponentFeatures:
    """Computed component-level features."""

    component_id: int
    size: int
    num_edges: int
    n_species: int
    species_entropy: float
    species_entropy_norm: float
    high_z_frac: Optional[float]
    max_z: Optional[float]
    top10_mean_z: Optional[float]


def compute_component_features(
    G,
    components: List[List[str]],
    z_attr: str = "z_robust",
    z0: float = 3.0,
) -> List[ComponentFeatures]:
    feats: List[ComponentFeatures] = []

    for cid, nodes in enumerate(components):
        node_set = set(nodes)
        size = len(nodes)

        species_counts = Counter(G.nodes[u].get("species", "") for u in nodes)
        n_species = sum(1 for _, c in species_counts.items() if c > 0)

        H = entropy_from_counts(species_counts)
        if n_species >= 2:
            H_norm = H / math.log(n_species)
        else:
            H_norm = 0.0

        num_edges = 0
        z_vals: List[float] = []
        high_z_count = 0

        for u in nodes:
            for v, edata in G[u].items():
                if v not in node_set:
                    continue
                if u >= v:
                    continue
                num_edges += 1
                z = edata.get(z_attr)
                if z is None:
                    continue
                zf = float(z)
                z_vals.append(zf)
                if zf >= z0:
                    high_z_count += 1

        if z_vals:
            high_z_frac = high_z_count / num_edges if num_edges > 0 else 0.0
            max_z = max(z_vals)
            top10_mean_z = topk_mean(z_vals, 10)
        else:
            high_z_frac = None
            max_z = None
            top10_mean_z = None

        feats.append(
            ComponentFeatures(
                component_id=cid,
                size=size,
                num_edges=num_edges,
                n_species=n_species,
                species_entropy=H,
                species_entropy_norm=H_norm,
                high_z_frac=high_z_frac,
                max_z=max_z,
                top10_mean_z=top10_mean_z,
            )
        )

    return feats


def write_component_features(path: Path, feats: List[ComponentFeatures]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "component_id",
                "size",
                "num_edges",
                "n_species",
                "species_entropy",
                "species_entropy_norm",
                "high_z_frac",
                "max_z",
                "top10_mean_z",
            ]
        )
        for c in feats:
            w.writerow(
                [
                    c.component_id,
                    c.size,
                    c.num_edges,
                    c.n_species,
                    f"{c.species_entropy:.6f}",
                    f"{c.species_entropy_norm:.6f}",
                    "" if c.high_z_frac is None else f"{c.high_z_frac:.6f}",
                    "" if c.max_z is None else f"{c.max_z:.6f}",
                    "" if c.top10_mean_z is None else f"{c.top10_mean_z:.6f}",
                ]
            )
