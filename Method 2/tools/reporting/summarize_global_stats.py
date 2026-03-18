#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Dict, Any, List, Optional

import numpy as np
import pandas as pd


def qstats(x: pd.Series, qs=(0.05, 0.25, 0.5, 0.75, 0.95)) -> Dict[str, float]:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) == 0:
        return {}
    out = {
        "n": int(len(x)),
        "min": float(x.min()),
        "max": float(x.max()),
        "mean": float(x.mean()),
        "std": float(x.std(ddof=1)) if len(x) > 1 else 0.0,
    }
    quant = x.quantile(list(qs))
    for q, v in quant.items():
        out[f"q{int(q*100):02d}"] = float(v)
    return out


def pct(cond: pd.Series) -> float:
    cond = cond.fillna(False)
    if len(cond) == 0:
        return float("nan")
    return 100.0 * float(cond.mean())


def pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def safe_read_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # normalize column names slightly (keep originals too)
    df.columns = [c.strip() for c in df.columns]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--component_features", required=True)
    ap.add_argument("--protein_features", required=True)
    ap.add_argument("--edge_features", required=True)
    ap.add_argument("--hgt_candidates", required=True)
    ap.add_argument("--out_prefix", default="global_stats")
    args = ap.parse_args()

    comp = safe_read_tsv(Path(args.component_features))
    prot = safe_read_tsv(Path(args.protein_features))
    edge = safe_read_tsv(Path(args.edge_features))
    cand = safe_read_tsv(Path(args.hgt_candidates))

    # ---- Column inference (handles minor naming differences) ----
    comp_size = pick_col(comp, ["size", "n_nodes", "component_size"])
    comp_edges = pick_col(comp, ["edges", "n_edges", "component_edges"])
    comp_density = pick_col(comp, ["density"])
    comp_species = pick_col(comp, ["n_species", "num_species", "species_count"])
    comp_H = pick_col(comp, ["H_norm", "entropy_norm", "normalized_entropy"])
    comp_highz = pick_col(comp, ["high_z_frac", "frac_high_z", "highz_frac"])
    comp_maxz = pick_col(comp, ["max_z", "z_max", "max_z_robust"])

    # Edge columns
    edge_z = pick_col(edge, ["z_robust", "z", "robust_z"])
    edge_jac = pick_col(edge, ["jaccard", "jac"])
    edge_shared = pick_col(edge, ["shared", "shared_kmers", "n_shared"])

    # Candidate score column: in your explain logs it prints "score="
    cand_score = pick_col(cand, ["score", "final_score", "S", "hgt_score"])

    # ---- Global stats container ----
    report: Dict[str, Any] = {"files": {k: getattr(args, k) for k in
                                        ["component_features", "protein_features", "edge_features", "hgt_candidates"]}}

    # ---- Component-level overview ----
    report["components"] = {
        "n_components": int(len(comp)),
        "size": qstats(comp[comp_size]) if comp_size else {},
        "edges": qstats(comp[comp_edges]) if comp_edges else {},
        "density": qstats(comp[comp_density]) if comp_density else {},
        "n_species": qstats(comp[comp_species]) if comp_species else {},
        "H_norm": qstats(comp[comp_H]) if comp_H else {},
        "high_z_frac": qstats(comp[comp_highz]) if comp_highz else {},
        "max_z": qstats(comp[comp_maxz]) if comp_maxz else {},
    }

    # Useful threshold summaries (edit thresholds here if you want)
    if comp_species:
        report["components"]["pct_mixed_ge3_species"] = pct(pd.to_numeric(comp[comp_species], errors="coerce") >= 3)
        report["components"]["pct_mixed_ge5_species"] = pct(pd.to_numeric(comp[comp_species], errors="coerce") >= 5)

    if comp_H:
        H = pd.to_numeric(comp[comp_H], errors="coerce")
        report["components"]["pct_entropy_gt_0p7"] = pct(H > 0.7)
        report["components"]["pct_entropy_gt_0p9"] = pct(H > 0.9)

    if comp_highz:
        hz = pd.to_numeric(comp[comp_highz], errors="coerce")
        report["components"]["pct_highz_gt_0p1"] = pct(hz > 0.1)
        report["components"]["pct_highz_gt_0p3"] = pct(hz > 0.3)

    # ---- Protein / gating overview ----
    report["proteins"] = {
        "n_proteins": int(len(prot)),
    }

    # candidates file typically includes only proteins in candidate set; treat score==0 as gated out
    if cand_score:
        s = pd.to_numeric(cand[cand_score], errors="coerce")
        report["candidates"] = {
            "n_candidates_rows": int(len(cand)),
            "score": qstats(s),
            "pct_scored_gt_0": pct(s > 0),
            "pct_scored_ge_1": pct(s >= 1),
            "pct_scored_ge_5": pct(s >= 5),
            "pct_scored_ge_10": pct(s >= 10),
        }
    else:
        report["candidates"] = {"n_candidates_rows": int(len(cand)), "note": "score column not found"}

    # ---- Edge-level overview ----
    report["edges"] = {
        "n_edges_rows": int(len(edge)),
        "z_robust": qstats(edge[edge_z]) if edge_z else {},
        "jaccard": qstats(edge[edge_jac]) if edge_jac else {},
        "shared": qstats(edge[edge_shared]) if edge_shared else {},
    }

    # ---- Pretty console output (paper-ready numbers) ----
    def fmt(v, nd=3):
        if v is None or (isinstance(v, float) and (np.isnan(v) or np.isinf(v))):
            return "NA"
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        return f"{float(v):.{nd}f}"

    print("\n=== GLOBAL SUMMARY ===")
    print(f"Components: {report['components']['n_components']}")
    print(f"Edges (rows): {report['edges']['n_edges_rows']}")
    print(f"Proteins (rows): {report['proteins']['n_proteins']}")
    if "n_candidates_rows" in report["candidates"]:
        print(f"Candidates (rows): {report['candidates']['n_candidates_rows']}")

    if comp_H and comp_highz and comp_species:
        print("\n--- Component mixing / anomaly concentration ---")
        Hs = report["components"]["H_norm"]
        Hz = report["components"]["high_z_frac"]
        Ns = report["components"]["n_species"]
        print(f"H_norm median [IQR]: {fmt(Hs.get('q50'))} [{fmt(Hs.get('q25'))}, {fmt(Hs.get('q75'))}]")
        print(f"high_z_frac median [IQR]: {fmt(Hz.get('q50'))} [{fmt(Hz.get('q25'))}, {fmt(Hz.get('q75'))}]")
        print(f"n_species median [IQR]: {fmt(Ns.get('q50'),0)} [{fmt(Ns.get('q25'),0)}, {fmt(Ns.get('q75'),0)}]")
        print(f"% comps with >=3 species: {fmt(report['components'].get('pct_mixed_ge3_species'),1)}%")
        print(f"% comps with H_norm>0.9: {fmt(report['components'].get('pct_entropy_gt_0p9'),1)}%")
        print(f"% comps with high_z_frac>0.3: {fmt(report['components'].get('pct_highz_gt_0p3'),1)}%")

    if cand_score:
        print("\n--- Candidate scoring / gating ---")
        cs = report["candidates"]["score"]
        print(f"Score median [IQR]: {fmt(cs.get('q50'))} [{fmt(cs.get('q25'))}, {fmt(cs.get('q75'))}]")
        print(f"Score 95th pct: {fmt(cs.get('q95'))}   max: {fmt(cs.get('max'))}")
        print(f"% candidates with score>0 (pass gates): {fmt(report['candidates'].get('pct_scored_gt_0'),1)}%")
        print(f"% candidates with score>=5: {fmt(report['candidates'].get('pct_scored_ge_5'),1)}%")
        print(f"% candidates with score>=10: {fmt(report['candidates'].get('pct_scored_ge_10'),1)}%")

    # ---- Save outputs ----
    out_json = Path(f"{args.out_prefix}.json")
    out_tsv = Path(f"{args.out_prefix}.tsv")

    out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")

    # Flatten for TSV (top-level + some key nested fields)
    rows = []
    def add_row(section, key, value):
        rows.append({"section": section, "metric": key, "value": value})

    add_row("components", "n_components", report["components"]["n_components"])
    for k in ["pct_mixed_ge3_species", "pct_mixed_ge5_species", "pct_entropy_gt_0p7", "pct_entropy_gt_0p9",
              "pct_highz_gt_0p1", "pct_highz_gt_0p3"]:
        if k in report["components"]:
            add_row("components", k, report["components"][k])

    if cand_score:
        add_row("candidates", "pct_scored_gt_0", report["candidates"]["pct_scored_gt_0"])
        add_row("candidates", "pct_scored_ge_5", report["candidates"]["pct_scored_ge_5"])
        add_row("candidates", "pct_scored_ge_10", report["candidates"]["pct_scored_ge_10"])

    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    print(f"\nWrote: {out_json} and {out_tsv}")


if __name__ == "__main__":
    main()
