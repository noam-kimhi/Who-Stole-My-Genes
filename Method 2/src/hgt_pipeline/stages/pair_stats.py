import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .edge_io import Edge


@dataclass
class PairRobustStats:
    """Robust baseline stats for a species-pair."""

    n: int
    median: float
    mad: float
    scale: float  # 1.4826*MAD + eps


@dataclass
class EdgeFeatures:
    """Computed edge-level features."""

    u: str
    v: str
    species_u: str
    species_v: str
    shared_kmers: int
    jaccard: float
    z_robust: Optional[float]  # None if insufficient data
    pair_median: Optional[float]
    pair_mad: Optional[float]
    pair_n: int


def median(xs: Sequence[float]) -> float:
    """Compute the median of a non-empty sequence of numbers."""
    if not xs:
        raise ValueError("median() arg is an empty sequence")

    xs_sorted = sorted(xs)
    n = len(xs_sorted)
    mid = n // 2

    if n % 2 == 1:
        return xs_sorted[mid]
    return 0.5 * (xs_sorted[mid - 1] + xs_sorted[mid])


def mad(xs: Sequence[float], med: float) -> float:
    """Compute the Median Absolute Deviation (MAD) from a given median."""
    if not xs:
        raise ValueError("mad() arg is an empty sequence")

    abs_devs = [abs(x - med) for x in xs]
    return median(abs_devs)


def compute_pair_robust_stats(
    edges: List[Edge],
    weight: str = "jaccard",
    min_pair_edges: int = 30,
    eps: float = 1e-9,
    shared_transform: str = "log1p",
) -> Dict[Tuple[str, str], PairRobustStats]:
    """
    For each unordered species pair (A,B), compute robust baseline stats over edge weights.
    """
    by_pair: Dict[Tuple[str, str], List[float]] = defaultdict(list)

    if weight not in {"jaccard", "shared"}:
        raise ValueError("weight must be 'jaccard' or 'shared'")
    if shared_transform not in {"raw", "log1p"}:
        raise ValueError("shared_transform must be 'raw' or 'log1p'")

    if weight == "jaccard":
        for e in edges:
            by_pair[e.species_pair()].append(e.jaccard)
    else:
        if shared_transform == "log1p":
            for e in edges:
                by_pair[e.species_pair()].append(math.log1p(e.shared_kmers))
        else:
            for e in edges:
                by_pair[e.species_pair()].append(float(e.shared_kmers))

    stats: Dict[Tuple[str, str], PairRobustStats] = {}
    for key, ws in by_pair.items():
        m = median(ws)
        d = mad(ws, m)
        s = 1.4826 * d + eps
        stats[key] = PairRobustStats(n=len(ws), median=m, mad=d, scale=s)

    return stats


def print_pair_stats_sanity_table(
    pair_stats: Dict[Tuple[str, str], PairRobustStats],
    min_pair_edges_for_z: int = 30,
    top_n: int = 15,
) -> None:
    """Print a quick sanity table to visually verify species-pair baselines."""
    rows = []
    for (a, b), st in pair_stats.items():
        ok = "YES" if st.n >= min_pair_edges_for_z else "NO"
        rows.append((st.n, a, b, st.median, st.mad, st.scale, ok))

    rows.sort(reverse=True, key=lambda t: t[0])
    rows = rows[:top_n]

    print("\n[PAIR-STATS] Robust baselines per species-pair (top by edge count)")
    print(f"{'n':>7}  {'species_A':<24}  {'species_B':<24}  {'median':>10}  {'MAD':>10}  {'scale':>10}  {'ok_for_z':>8}")
    print("-" * (7 + 2 + 24 + 2 + 24 + 2 + 10 + 2 + 10 + 2 + 10 + 2 + 8))
    for n, a, b, med, d, s, ok in rows:
        print(f"{n:7d}  {a:<24.24}  {b:<24.24}  {med:10.6f}  {d:10.6f}  {s:10.6f}  {ok:>8}")
    print()


def compute_edge_features(
    edges: List[Edge],
    pair_stats: Dict[Tuple[str, str], PairRobustStats],
    weight: str = "jaccard",
    min_pair_edges_for_z: int = 30,
) -> List[EdgeFeatures]:
    """Compute edge-level features, including robust z-score per species-pair."""
    out: List[EdgeFeatures] = []

    if weight not in {"jaccard", "shared"}:
        raise ValueError("weight must be 'jaccard' or 'shared'")

    for e in edges:
        key = e.species_pair()
        st = pair_stats.get(key)

        if weight == "jaccard":
            w = e.jaccard
        else:
            w = float(e.shared_kmers)

        if st is None:
            z = None
            pair_median = None
            pair_mad = None
            pair_n = 0
        else:
            pair_median = st.median
            pair_mad = st.mad
            pair_n = st.n

            if st.n < min_pair_edges_for_z:
                z = None
            else:
                z = (w - st.median) / st.scale

        out.append(
            EdgeFeatures(
                u=e.u,
                v=e.v,
                species_u=e.species_u,
                species_v=e.species_v,
                shared_kmers=e.shared_kmers,
                jaccard=e.jaccard,
                z_robust=z,
                pair_median=pair_median,
                pair_mad=pair_mad,
                pair_n=pair_n,
            )
        )

    return out


def write_edge_features(path: Path, feats: List[EdgeFeatures]) -> None:
    """Write edge_features.tsv with robust z and per-pair baseline fields."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "u",
                "v",
                "species_u",
                "species_v",
                "shared_kmers",
                "jaccard",
                "z_robust",
                "pair_median",
                "pair_mad",
                "pair_n",
            ]
        )
        for e in feats:
            w.writerow(
                [
                    e.u,
                    e.v,
                    e.species_u,
                    e.species_v,
                    e.shared_kmers,
                    f"{e.jaccard:.6f}",
                    "" if e.z_robust is None else f"{e.z_robust:.6f}",
                    "" if e.pair_median is None else f"{e.pair_median:.6f}",
                    "" if e.pair_mad is None else f"{e.pair_mad:.6f}",
                    e.pair_n,
                ]
            )
