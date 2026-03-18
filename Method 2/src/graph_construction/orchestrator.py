#!/usr/bin/env python3
"""Orchestrate graph-construction recipes (FASTA -> candidates -> pruned edges)."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from graph_construction.fasta_parsing import load_manifest_species_map
from graph_construction.graph_pruning import (
    keep_q_percentile_edges,
    keep_top_X_edges_per_node,
)
from graph_construction.kmer_candidates_from_faa import (
    build_kmer_index,
    generate_candidates,
    load_proteins_from_downloads,
)


def build_candidates(
    manifest: Path,
    downloads_dir: Path,
    out_candidates: Path,
    k: int,
    min_len: int,
    max_postings: int,
    min_shared: int,
    top_m: int,
    cross_species_only: bool,
) -> None:
    dirname_to_species = load_manifest_species_map(manifest)
    prot_uid, species_id, species_name, seqs = load_proteins_from_downloads(
        downloads_dir=downloads_dir,
        dirname_to_species=dirname_to_species,
        min_len=min_len,
    )
    kmers_by_protein, kmer_index = build_kmer_index(
        prot_uid=prot_uid,
        species_id=species_id,
        seqs=seqs,
        k=k,
        max_postings=max_postings,
    )
    generate_candidates(
        prot_uid=prot_uid,
        species_id=species_id,
        species_name=species_name,
        K=kmers_by_protein,
        index=kmer_index,
        min_shared=min_shared,
        top_m=top_m,
        out_path=out_candidates,
        cross_species_only=cross_species_only,
    )


def prune_candidates(
    in_candidates: Path,
    out_edges: Path,
    q: float,
    top_x: int,
) -> None:
    df = pd.read_csv(in_candidates, sep="\t")
    filtered = keep_q_percentile_edges(df, q=q)
    pruned = keep_top_X_edges_per_node(filtered, X=top_x)
    out_edges.parent.mkdir(parents=True, exist_ok=True)
    pruned.to_csv(out_edges, sep="\t", index=False)
    print(f"[OK] wrote pruned edges: {out_edges}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Graph-construction orchestration recipes.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_build = sub.add_parser("build-candidates", help="Build candidates.tsv from manifest + FASTA downloads.")
    p_build.add_argument("--manifest", type=Path, required=True)
    p_build.add_argument("--downloads_dir", type=Path, required=True)
    p_build.add_argument("--out_candidates", type=Path, required=True)
    p_build.add_argument("--k", type=int, default=6)
    p_build.add_argument("--min_len", type=int, default=50)
    p_build.add_argument("--max_postings", type=int, default=100)
    p_build.add_argument("--min_shared", type=int, default=6)
    p_build.add_argument("--top_m", type=int, default=10)
    p_build.add_argument("--allow_within_species", action="store_true")

    p_prune = sub.add_parser("prune-candidates", help="Prune candidates.tsv into edges.tsv.")
    p_prune.add_argument("--in_candidates", type=Path, required=True)
    p_prune.add_argument("--out_edges", type=Path, required=True)
    p_prune.add_argument("--q", type=float, default=0.9)
    p_prune.add_argument("--top_x", type=int, default=20)

    p_construct = sub.add_parser("construct-edges", help="Run build-candidates then prune-candidates.")
    p_construct.add_argument("--manifest", type=Path, required=True)
    p_construct.add_argument("--downloads_dir", type=Path, required=True)
    p_construct.add_argument("--out_candidates", type=Path, required=True)
    p_construct.add_argument("--out_edges", type=Path, required=True)
    p_construct.add_argument("--k", type=int, default=6)
    p_construct.add_argument("--min_len", type=int, default=50)
    p_construct.add_argument("--max_postings", type=int, default=100)
    p_construct.add_argument("--min_shared", type=int, default=6)
    p_construct.add_argument("--top_m", type=int, default=10)
    p_construct.add_argument("--q", type=float, default=0.9)
    p_construct.add_argument("--top_x", type=int, default=20)
    p_construct.add_argument("--allow_within_species", action="store_true")

    args = ap.parse_args()

    if args.cmd == "build-candidates":
        build_candidates(
            manifest=args.manifest,
            downloads_dir=args.downloads_dir,
            out_candidates=args.out_candidates,
            k=args.k,
            min_len=args.min_len,
            max_postings=args.max_postings,
            min_shared=args.min_shared,
            top_m=args.top_m,
            cross_species_only=not args.allow_within_species,
        )
        return

    if args.cmd == "prune-candidates":
        prune_candidates(
            in_candidates=args.in_candidates,
            out_edges=args.out_edges,
            q=args.q,
            top_x=args.top_x,
        )
        return

    if args.cmd == "construct-edges":
        build_candidates(
            manifest=args.manifest,
            downloads_dir=args.downloads_dir,
            out_candidates=args.out_candidates,
            k=args.k,
            min_len=args.min_len,
            max_postings=args.max_postings,
            min_shared=args.min_shared,
            top_m=args.top_m,
            cross_species_only=not args.allow_within_species,
        )
        prune_candidates(
            in_candidates=args.out_candidates,
            out_edges=args.out_edges,
            q=args.q,
            top_x=args.top_x,
        )
        return


if __name__ == "__main__":
    main()
