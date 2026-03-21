#!/usr/bin/env python3
"""
kmer_candidates_from_faa.py

Build k-mer inverted index (Option A) and generate candidate pairs from downloaded
RefSeq assembly protein FASTAs (*.faa.gz).

Inputs:
  - manifest.tsv (from refseq_fetch_proteins.py)
  - downloads_dir containing *_protein.faa.gz files

Outputs:
  - candidates.tsv: sparse candidate edges with shared-kmer count and Jaccard

Example run:
  python kmer_candidates_from_faa.py \
    --manifest data/out_refseq/manifest.tsv \
    --downloads_dir data/out_refseq/downloads \
    --out candidates.tsv \
    --k 5 --min_len 50 --max_postings 2000 --min_shared 3 --top_m 50
"""

import argparse
import csv
import time
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple, Union

try:
    from graph_construction.k_mer_encoding import kmers_encoded_set
    from graph_construction.fasta_parsing import (
        iter_fasta_gz,
        protein_id_from_header,
        load_manifest_species_map,
    )
except ImportError:
    from k_mer_encoding import kmers_encoded_set
    from fasta_parsing import iter_fasta_gz, protein_id_from_header, load_manifest_species_map

# ---------------- Build protein list, K, and inverted index ----------------
def load_proteins_from_downloads(
    downloads_dir: Path,
    dirname_to_species: Dict[str, str],
    min_len: int,
) -> Tuple[List[str], List[int], List[str], List[str]]:
    """
    Returns:
      prot_uid: list[str]  unique protein id (dirname|protein_accession)
      species_id: list[int]
      species_name: list[str]  (parallel list)
      seqs: list[str]
    """
    # assign small int to each species
    species2id: Dict[str, int] = {}
    species_name: List[str] = []

    prot_uid: List[str] = []
    species_id: List[int] = []
    seqs: List[str] = []

    faa_files = sorted(downloads_dir.glob("*_protein.faa.gz"))
    if not faa_files:
        raise RuntimeError(f"No *_protein.faa.gz files found in {downloads_dir}")

    for faa in faa_files:
        fname = faa.name
        if not fname.endswith("_protein.faa.gz"):
            continue

        dirname = fname[:-len("_protein.faa.gz")]
        sp = dirname_to_species.get(dirname)
        if sp is None:
            # not in manifest; skip
            continue

        if sp not in species2id:
            species2id[sp] = len(species2id)
            species_name.append(sp)

        spid = species2id[sp]

        print(f"[LOAD] {fname}  ->  {sp}")
        n_in_file = 0
        n_kept = 0
        for header, seq in iter_fasta_gz(faa):
            n_in_file += 1
            if len(seq) < min_len:
                continue
            pid = protein_id_from_header(header)
            uid = f"{dirname}|{pid}"
            prot_uid.append(uid)
            species_id.append(spid)
            seqs.append(seq)
            n_kept += 1

        print(f"       proteins: {n_in_file} records, kept: {n_kept} (min_len={min_len})")

    if not prot_uid:
        raise RuntimeError("Loaded 0 proteins. Check manifest.tsv mapping and downloads_dir.")
    return prot_uid, species_id, species_name, seqs


def build_kmer_index(
    prot_uid: List[str],
    species_id: List[int],
    seqs: List[str],
    k: int,
    max_postings: int,
) -> tuple[list[list[int]], dict[int, list[int]]]:
    index: Dict[int, List[int]] = defaultdict(list)
    K: List[List[int]] = []

    total = len(seqs)
    for i, seq in enumerate(seqs):
        if (i + 1) % 5000 == 0 or i == 0 or (i + 1) == total:
            print(f"[KMERS] {i+1}/{total}")

        ks_set = kmers_encoded_set(seq, k)
        ks = sorted(ks_set)          # list, smaller than set long-term
        K.append(ks)

        for code in ks:
            index[code].append(i)

    before = len(index)
    # Posting list statistics (diagnostic)
    if index:
        lens = sorted(len(lst) for lst in index.values())
        n = len(lens)
        print(
            f"[INDEX] posting sizes: "
            f"max={lens[-1]}, "
            f"p99={lens[int(0.99 * (n - 1))]}, "
            f"p95={lens[int(0.95 * (n - 1))]}, "
            f"median={lens[n // 2]}"
        )

    if max_postings > 0:
        index = {code: lst for code, lst in index.items() if len(lst) <= max_postings}
        kept = set(index.keys())
        for i in range(len(K)):
            # filter list in-place by rebuilding (fast enough)
            K[i] = [x for x in K[i] if x in kept]

    after = len(index)
    print(f"[INDEX] kmers before prune: {before:,}, after prune: {after:,} (max_postings={max_postings})")

    return K, index


# ---------------- Candidate generation ----------------

def jaccard(intersection: int, a: int, b: int) -> float:
    denom = a + b - intersection
    return 0.0 if denom <= 0 else intersection / denom

def generate_candidates(
    prot_uid: List[str],
    species_id: List[int],
    species_name: List[str],
    K: Union[List[Set[int]], List[List[int]]],
    index: Dict[int, List[int]],
    min_shared: int,
    top_m: int,
    out_path: Path,
    cross_species_only: bool = True,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    total = len(prot_uid)
    k_sizes = [len(ks) for ks in K]

    idx_get = index.get
    sid = species_id

    # metrics
    proteins_with_edges = 0
    total_items_after_filter = 0
    max_items_after_filter = 0

    total_posting_visits = 0
    codes_touched = 0

    min_inter = 10**9
    max_inter = 0
    sum_inter = 0

    proteins_with_lists = 0
    proteins_hitting_cap = 0

    # store best per unordered pair: (a,b)->(shared, jaccard)
    best_edge: Dict[Tuple[int, int], Tuple[int, float]] = {}

    for p in range(total):
        if (p + 1) % 2000 == 0 or p == 0 or (p + 1) == total:
            print(f"[CAND] {p+1}/{total}")

        sp_p = sid[p]
        shared = Counter()

        for code in K[p]:
            post = idx_get(code)
            if not post:
                continue

            codes_touched += 1
            total_posting_visits += len(post)

            for q in post:
                if q == p:
                    continue
                if cross_species_only and sid[q] == sp_p:
                    continue
                shared[q] += 1

        if not shared:
            continue

        items = [(q, c) for q, c in shared.items() if c >= min_shared]
        if not items:
            continue

        items.sort(key=lambda t: t[1], reverse=True)
        items = items[:top_m]

        proteins_with_lists += 1
        if len(items) == top_m:
            proteins_hitting_cap += 1

        proteins_with_edges += 1
        total_items_after_filter += len(items)
        if len(items) > max_items_after_filter:
            max_items_after_filter = len(items)

        a_sz = k_sizes[p]
        for q, inter in items:
            b_sz = k_sizes[q]
            jac = jaccard(inter, a_sz, b_sz)

            a = p
            b = q
            if a > b:
                a, b = b, a

            prev = best_edge.get((a, b))
            if prev is None or inter > prev[0] or (inter == prev[0] and jac > prev[1]):
                best_edge[(a, b)] = (inter, jac)

            # shared stats
            if inter < min_inter:
                min_inter = inter
            if inter > max_inter:
                max_inter = inter
            sum_inter += inter

    # write once
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["u", "v", "shared_kmers", "jaccard", "species_u", "species_v"])

        for (a, b), (inter, jac) in best_edge.items():
            u = prot_uid[a]
            v = prot_uid[b]
            sp_u = species_name[sid[a]]
            sp_v = species_name[sid[b]]
            w.writerow([u, v, inter, f"{jac:.6f}", sp_u, sp_v])

    edges_written = len(best_edge)
    print(f"[OK] wrote candidates: {out_path}  edges_undirected={edges_written:,}")
    print(f"[CAND-STATS] proteins with >=1 candidate list: {proteins_with_edges:,}/{total:,}")
    if proteins_with_lists:
        frac = proteins_hitting_cap / proteins_with_lists
        print(f"[CAND-STATS] proteins hitting top_m cap: {proteins_hitting_cap:,}/{proteins_with_lists:,} ({frac:.1%})")
    if proteins_with_edges:
        print(f"[CAND-STATS] avg kept per protein (among those with lists): "
              f"{total_items_after_filter / proteins_with_edges:.2f} (top_m={top_m}), max={max_items_after_filter}")
    if edges_written:
        print(f"[CAND-STATS] shared_kmers: min={min_inter}, mean={sum_inter / max(1, (proteins_with_edges * max(1, int(total_items_after_filter / proteins_with_edges)))):.2f}, max={max_inter}")
    print(f"[CAND-WORK] codes_touched={codes_touched:,}, total_posting_visits={total_posting_visits:,}, "
          f"avg_post_len={total_posting_visits / max(codes_touched,1):.2f}")




# ---------------- CLI ----------------
def main():
    ap = argparse.ArgumentParser(description="From *.faa.gz to inverted index and k-mer candidates (Option A).")
    ap.add_argument("--manifest", required=True, type=Path, help="Path to manifest.tsv")
    ap.add_argument("--downloads_dir", required=True, type=Path, help="Directory with *_protein.faa.gz")
    ap.add_argument("--out", required=True, type=Path, help="Output candidates.tsv")
    ap.add_argument("--k", type=int, default=5, help="k-mer length (default: 5)")
    ap.add_argument("--min_len", type=int, default=50, help="min protein length (default: 50)")
    ap.add_argument("--max_postings", type=int, default=2000,
                    help="drop kmers whose posting list > this (default: 2000). Use 0 to disable.")
    ap.add_argument("--min_shared", type=int, default=3, help="min shared kmers for candidate (default: 3)")
    ap.add_argument("--top_m", type=int, default=50, help="keep top-M candidates per protein (default: 50)")
    ap.add_argument("--cross_species_only", action="store_true", help="keep only pairs across different species")
    args = ap.parse_args()

    dirname_to_species = load_manifest_species_map(args.manifest)

    print(f"[CFG] k={args.k} max_postings={args.max_postings} min_shared={args.min_shared} top_m={args.top_m} cross_species_only={args.cross_species_only}")

    # ---------------- LOAD ----------------
    t0 = time.perf_counter()
    prot_uid, species_id, species_name, seqs = load_proteins_from_downloads(
        downloads_dir=args.downloads_dir,
        dirname_to_species=dirname_to_species,
        min_len=args.min_len,
    )
    t1 = time.perf_counter()

    total_aa = sum(len(s) for s in seqs)
    avg_len = total_aa / len(seqs)
    print(f"[INFO] avg protein length = {avg_len:.1f} aa")
    print(f"[INFO] expected max windows/protein for k={args.k}: ~{max(avg_len - args.k + 1, 0):.1f}")

    print(f"[INFO] loaded proteins: {len(prot_uid)} across {len(species_name)} species")
    print(f"[INFO] input size: total amino acids = {total_aa:,}")
    print(f"[TIME] load FASTA + filter: {t1 - t0:.3f}s")

    sp_counts = Counter(species_id)
    print("[INFO] proteins per species:")
    for spid, cnt in sorted(sp_counts.items(), key=lambda t: t[1], reverse=True):
        print(f"  - {species_name[spid]}: {cnt:,}")

    # ---------------- KMERS + INDEX ----------------
    t2 = time.perf_counter()
    K, index = build_kmer_index(
        prot_uid=prot_uid,
        species_id=species_id,
        seqs=seqs,
        k=args.k,
        max_postings=args.max_postings,
    )
    t3 = time.perf_counter()

    # These summaries are helpful; cost is small compared to building the index.
    total_unique_kmers_summed = sum(len(ks) for ks in K)
    print(f"[INFO] kmers summary: sum(|K[p]|) = {total_unique_kmers_summed:,}, avg/protein = {total_unique_kmers_summed / max(len(K),1):.1f}")
    print(f"[INFO] index summary: distinct kmers = {len(index):,}")
    print(f"[TIME] k-mers + inverted index (incl prune): {t3 - t2:.3f}s")

    # ---------------- CANDIDATES ----------------
    t4 = time.perf_counter()
    generate_candidates(
        prot_uid=prot_uid,
        species_id=species_id,
        species_name=species_name,
        K=K,
        index=index,
        min_shared=args.min_shared,
        top_m=args.top_m,
        out_path=args.out,
        cross_species_only=args.cross_species_only
    )
    t5 = time.perf_counter()

    print(f"[TIME] candidate generation + write: {t5 - t4:.3f}s")
    print(f"[TIME] total pipeline: {t5 - t0:.3f}s")


if __name__ == "__main__":
    main()
