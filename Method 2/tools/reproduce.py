#!/usr/bin/env python3
"""Root-level reproducibility runner for end-to-end pipelines."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List


REPO_ROOT = Path(__file__).resolve().parent
SRC_DIR = REPO_ROOT / "src"
TOOLS_REPORTING_DIR = REPO_ROOT / "tools" / "reporting"


def run_step(args: List[str], cwd: Path) -> None:
    cmd = [sys.executable, *args]
    print(f"[RUN] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(cwd), check=False)
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def run_pipeline_and_reports(
    edges_path: Path,
    out_dir: Path,
    with_betweenness: bool,
    with_reports: bool,
    with_explanations: bool,
    component_ids: Iterable[int],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    pipeline_cmd = [
        "-m",
        "hgt_pipeline.pipeline",
        "--in_edges",
        str(edges_path.resolve()),
        "--out_dir",
        str(out_dir.resolve()),
    ]
    if not with_betweenness:
        pipeline_cmd.append("--no_betweenness")
    run_step(pipeline_cmd, cwd=SRC_DIR)

    if not with_reports:
        return

    results_dir = out_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    run_step(
        [
            str((TOOLS_REPORTING_DIR / "top_anomaly_edges.py").resolve()),
            "--edges",
            str((out_dir / "edge_features.tsv").resolve()),
            "--top_n",
            "25",
            "--out_dir",
            str(results_dir.resolve()),
        ],
        cwd=REPO_ROOT,
    )
    run_step(
        [
            str((TOOLS_REPORTING_DIR / "summarize_global_stats.py").resolve()),
            "--component_features",
            str((out_dir / "component_features.tsv").resolve()),
            "--protein_features",
            str((out_dir / "protein_features.tsv").resolve()),
            "--edge_features",
            str((out_dir / "edge_features.tsv").resolve()),
            "--hgt_candidates",
            str((out_dir / "hgt_candidates.tsv").resolve()),
            "--out_prefix",
            str((results_dir / "global_stats").resolve()),
        ],
        cwd=REPO_ROOT,
    )

    if not with_explanations:
        return

    expl_dir = out_dir / "explanations"
    expl_dir.mkdir(parents=True, exist_ok=True)

    for cid in component_ids:
        out_file = expl_dir / f"comp_{cid}_explained.txt"
        with out_file.open("w", encoding="utf-8", newline="\n") as f:
            cmd = [
                sys.executable,
                str((TOOLS_REPORTING_DIR / "explain_component.py").resolve()),
                "--component_id",
                str(cid),
                "--edges",
                str((out_dir / "edge_features.tsv").resolve()),
                "--protein_features",
                str((out_dir / "protein_features.tsv").resolve()),
                "--component_features",
                str((out_dir / "component_features.tsv").resolve()),
                "--hgt_candidates",
                str((out_dir / "hgt_candidates.tsv").resolve()),
                "--top_nodes",
                "20",
                "--top_edges",
                "25",
            ]
            print(f"[RUN] {' '.join(cmd)} > {out_file}")
            result = subprocess.run(cmd, cwd=str(REPO_ROOT), stdout=f, check=False)
            if result.returncode != 0:
                raise SystemExit(result.returncode)

    top_file = expl_dir / "top_candidates_explained.txt"
    with top_file.open("w", encoding="utf-8", newline="\n") as f:
        cmd = [
            sys.executable,
            str((TOOLS_REPORTING_DIR / "explain_top_candidates.py").resolve()),
            "--edges",
            str((out_dir / "edge_features.tsv").resolve()),
            "--protein_features",
            str((out_dir / "protein_features.tsv").resolve()),
            "--component_features",
            str((out_dir / "component_features.tsv").resolve()),
            "--hgt_candidates",
            str((out_dir / "hgt_candidates.tsv").resolve()),
            "--top_n",
            "20",
            "--top_k_neighbors",
            "12",
        ]
        print(f"[RUN] {' '.join(cmd)} > {top_file}")
        result = subprocess.run(cmd, cwd=str(REPO_ROOT), stdout=f, check=False)
        if result.returncode != 0:
            raise SystemExit(result.returncode)


def main() -> None:
    ap = argparse.ArgumentParser(description="End-to-end reproducibility recipes.")
    sub = ap.add_subparsers(dest="recipe", required=True)

    p_edges = sub.add_parser("from-edges", help="Run HGT pipeline + reports starting from an edges TSV.")
    p_edges.add_argument("--in_edges", type=Path, required=True)
    p_edges.add_argument("--out_dir", type=Path, required=True)
    p_edges.add_argument("--with_betweenness", action="store_true")
    p_edges.add_argument("--with_reports", action="store_true")
    p_edges.add_argument("--with_explanations", action="store_true")
    p_edges.add_argument("--component_ids", type=int, nargs="*", default=[5, 8, 32])

    p_manifest = sub.add_parser(
        "from-manifest",
        help="Run graph construction (manifest+downloads -> edges) then HGT pipeline + reports.",
    )
    p_manifest.add_argument("--manifest", type=Path, required=True)
    p_manifest.add_argument("--downloads_dir", type=Path, required=True)
    p_manifest.add_argument("--work_dir", type=Path, required=True)
    p_manifest.add_argument("--k", type=int, default=6)
    p_manifest.add_argument("--min_len", type=int, default=50)
    p_manifest.add_argument("--max_postings", type=int, default=100)
    p_manifest.add_argument("--min_shared", type=int, default=6)
    p_manifest.add_argument("--top_m", type=int, default=10)
    p_manifest.add_argument("--q", type=float, default=0.9)
    p_manifest.add_argument("--top_x", type=int, default=20)
    p_manifest.add_argument("--with_betweenness", action="store_true")
    p_manifest.add_argument("--with_reports", action="store_true")
    p_manifest.add_argument("--with_explanations", action="store_true")
    p_manifest.add_argument("--component_ids", type=int, nargs="*", default=[5, 8, 32])

    p_summary = sub.add_parser(
        "from-assembly-summary",
        help="Full E2E from assembly_summary_refseq.txt + species list (downloads included).",
    )
    p_summary.add_argument("--assembly_summary", type=Path, required=True)
    p_summary.add_argument("--species_list", type=Path, required=True)
    p_summary.add_argument("--work_dir", type=Path, required=True)
    p_summary.add_argument("--max_assemblies_per_species", type=int, default=2)
    p_summary.add_argument("--with_betweenness", action="store_true")
    p_summary.add_argument("--with_reports", action="store_true")
    p_summary.add_argument("--with_explanations", action="store_true")
    p_summary.add_argument("--component_ids", type=int, nargs="*", default=[5, 8, 32])

    args = ap.parse_args()

    if args.recipe == "from-edges":
        run_pipeline_and_reports(
            edges_path=args.in_edges,
            out_dir=args.out_dir,
            with_betweenness=args.with_betweenness,
            with_reports=args.with_reports,
            with_explanations=args.with_explanations,
            component_ids=args.component_ids,
        )
        return

    if args.recipe == "from-manifest":
        work_dir: Path = args.work_dir
        work_dir.mkdir(parents=True, exist_ok=True)
        candidates = work_dir / "candidates.tsv"
        edges = work_dir / "edges_pruned.tsv"
        run_step(
            [
                "-m",
                "graph_construction.orchestrator",
                "construct-edges",
                "--manifest",
                str(args.manifest.resolve()),
                "--downloads_dir",
                str(args.downloads_dir.resolve()),
                "--out_candidates",
                str(candidates.resolve()),
                "--out_edges",
                str(edges.resolve()),
                "--k",
                str(args.k),
                "--min_len",
                str(args.min_len),
                "--max_postings",
                str(args.max_postings),
                "--min_shared",
                str(args.min_shared),
                "--top_m",
                str(args.top_m),
                "--q",
                str(args.q),
                "--top_x",
                str(args.top_x),
            ],
            cwd=SRC_DIR,
        )
        run_pipeline_and_reports(
            edges_path=edges,
            out_dir=work_dir / "pipeline",
            with_betweenness=args.with_betweenness,
            with_reports=args.with_reports,
            with_explanations=args.with_explanations,
            component_ids=args.component_ids,
        )
        return

    if args.recipe == "from-assembly-summary":
        work_dir = args.work_dir
        gc_out = work_dir / "graph_construction"
        gc_out.mkdir(parents=True, exist_ok=True)
        run_step(
            [
                "-m",
                "graph_construction.refseq_fetch_proteins",
                "--assembly_summary",
                str(args.assembly_summary.resolve()),
                "--species_list",
                str(args.species_list.resolve()),
                "--out_dir",
                str(gc_out.resolve()),
                "--max_assemblies_per_species",
                str(args.max_assemblies_per_species),
                "--require_latest",
                "--download",
            ],
            cwd=SRC_DIR,
        )
        manifest = gc_out / "manifest.tsv"
        downloads = gc_out / "downloads"
        cmd = [
            sys.executable,
            str((REPO_ROOT / "reproduce.py").resolve()),
            "from-manifest",
            "--manifest",
            str(manifest.resolve()),
            "--downloads_dir",
            str(downloads.resolve()),
            "--work_dir",
            str((work_dir / "run").resolve()),
        ]
        if args.with_betweenness:
            cmd.append("--with_betweenness")
        if args.with_reports:
            cmd.append("--with_reports")
        if args.with_explanations:
            cmd.append("--with_explanations")
        if args.component_ids:
            cmd.extend(["--component_ids", *[str(x) for x in args.component_ids]])
        print(f"[RUN] {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=str(REPO_ROOT), check=False)
        if result.returncode != 0:
            raise SystemExit(result.returncode)
        return


if __name__ == "__main__":
    main()

