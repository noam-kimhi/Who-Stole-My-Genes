from Bio.Align import PairwiseAligner
from typing import Dict, List, Tuple
from ..constants import (DEFAULT_MIN_ID, DEFAULT_MIN_COVERAGE, MIN_ALIGNED_LENGTH, EDGE_WEIGHT_KEY,
                         IDENTITY_KEY, COVERAGE_KEY, ALIGNED_LENGTH_KEY)
from .local_alignment import local_alignment_metrics
from .aligner import make_aligner


def build_edges_from_alignments(seqs: Dict[str, str],
                                min_identity: float = DEFAULT_MIN_ID,
                                min_coverage: float = DEFAULT_MIN_COVERAGE,
                                min_aligned_length: int = MIN_ALIGNED_LENGTH) -> List[Tuple[str, str, Dict[str, float]]]:
    """
    Build edges based on local alignment metrics between sequences.
    :param seqs: A dictionary of sequences with protein IDs as keys.
    :param min_identity: The minimum identity threshold.
    :param min_coverage: The minimum coverage threshold.
    :param min_aligned_length: The minimum aligned length threshold.
    :return: A list of edges with weights.
    """
    ids = sorted(seqs.keys())
    edges: List[Tuple[str, str, Dict[str, float]]] = []

    aligner: PairwiseAligner = make_aligner()

    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            u, v = ids[i], ids[j]
            s1, s2 = seqs[u], seqs[v]

            identity, coverage, aligned_len = local_alignment_metrics(aligner, s1, s2)

            if identity >= min_identity and coverage >= min_coverage and aligned_len >= min_aligned_length:
                similarity_score: float = identity * coverage
                edges.append((
                    u, v,
                    {
                        EDGE_WEIGHT_KEY: similarity_score,
                        IDENTITY_KEY: identity,
                        COVERAGE_KEY: coverage,
                        ALIGNED_LENGTH_KEY: aligned_len,
                    }
                ))

    return edges
