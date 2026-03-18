import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class Edge:
    """Undirected edge record."""

    u: str
    v: str
    shared_kmers: int
    jaccard: float
    species_u: str
    species_v: str

    def species_pair(self) -> Tuple[str, str]:
        """Return unordered species-pair key (A,B) with A <= B."""
        a, b = self.species_u, self.species_v
        return (a, b) if a <= b else (b, a)

    def node_pair(self) -> Tuple[str, str]:
        """Return unordered node-pair key (u,v) with u <= v."""
        a, b = self.u, self.v
        return (a, b) if a <= b else (b, a)


def sniff_delimiter(path: Path) -> str:
    """
    Try to detect delimiter (',' or '\\t') from the first non-empty line.
    """
    with open(path, "r", encoding="utf-8", newline="") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "\t" in line and "," not in line:
                return "\t"
            if "," in line and "\t" not in line:
                return ","
            # If both exist, fallback to CSV logic:
            return ","
    return ","


def read_edges(path: Path) -> List[Edge]:
    """
    Read edge list from CSV/TSV with header:
      u, v, shared_kmers, jaccard, species_u, species_v

    Returns a list of Edge records.

    NOTE: This assumes the file represents undirected edges already.
    If duplicates exist (u,v) and (v,u), you should dedupe (see dedupe_edges()).
    """
    delim = sniff_delimiter(path)
    edges: List[Edge] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter=delim)
        required = {"u", "v", "shared_kmers", "jaccard", "species_u", "species_v"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}. Found: {r.fieldnames}")

        for row in r:
            edges.append(
                Edge(
                    u=row["u"],
                    v=row["v"],
                    shared_kmers=int(row["shared_kmers"]),
                    jaccard=float(row["jaccard"]),
                    species_u=row["species_u"],
                    species_v=row["species_v"],
                )
            )
    return edges


def dedupe_edges(edges: List[Edge]) -> List[Edge]:
    """
    Dedupe duplicates for an undirected graph.

    Policy:
      For each unordered protein pair {u, v},
      keep exactly one edge - the one with the maximum Jaccard similarity.
    """
    best: Dict[Tuple[str, str], Edge] = {}
    for e in edges:
        key = e.node_pair()
        prev = best.get(key)
        if prev is None:
            best[key] = e
        else:
            if e.jaccard > prev.jaccard:
                best[key] = e
    return list(best.values())

