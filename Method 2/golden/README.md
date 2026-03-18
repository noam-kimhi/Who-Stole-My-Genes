# Golden Artifacts

This directory is the canonical regression source for tidy/refactor work.

## Canonical Input

- `reference_inputs/edges_PRUNED_JACCARD_92790.tsv`

## Canonical No-Betweenness Pipeline Baseline

Generated on `2026-03-18` with:

```bash
python hgt_pipeline/pipeline.py \
  --in_edges ./data/out_refseq/edges_PRUNED_JACCARD_92790.tsv \
  --out_dir ./golden/no_bw_pipeline/rerun_pruned \
  --no_betweenness
```

Outputs:

- `no_bw_pipeline/rerun_pruned/edge_features.tsv`
- `no_bw_pipeline/rerun_pruned/protein_features.tsv`
- `no_bw_pipeline/rerun_pruned/component_features.tsv`
- `no_bw_pipeline/rerun_pruned/hgt_candidates.tsv`
- `no_bw_pipeline/rerun_pruned/all_scores.tsv`

Derived report artifacts:

```bash
python reporting/top_anomaly_edges.py \
  --edges ./golden/no_bw_pipeline/rerun_pruned/edge_features.tsv \
  --top_n 25 \
  --out_dir ./golden/no_bw_pipeline/results

python reporting/summarize_global_stats.py \
  --component_features ./golden/no_bw_pipeline/rerun_pruned/component_features.tsv \
  --protein_features ./golden/no_bw_pipeline/rerun_pruned/protein_features.tsv \
  --edge_features ./golden/no_bw_pipeline/rerun_pruned/edge_features.tsv \
  --hgt_candidates ./golden/no_bw_pipeline/rerun_pruned/hgt_candidates.tsv \
  --out_prefix ./golden/no_bw_pipeline/results/global_stats
```

Derived outputs:

- `no_bw_pipeline/results/top_anomaly_edges.tsv`
- `no_bw_pipeline/results/top_anomaly_edges.png`
- `no_bw_pipeline/results/global_stats.tsv`
- `no_bw_pipeline/results/global_stats.json`

## Canonical Betweenness-On Pipeline Baseline

Generated on `2026-03-18` with:

```bash
python hgt_pipeline/pipeline.py \
  --in_edges ./golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv \
  --out_dir ./golden/bw_pipeline/rerun_pruned
```

Derived report artifacts:

```bash
python reporting/top_anomaly_edges.py \
  --edges ./golden/bw_pipeline/rerun_pruned/edge_features.tsv \
  --top_n 25 \
  --out_dir ./golden/bw_pipeline/results

python reporting/summarize_global_stats.py \
  --component_features ./golden/bw_pipeline/rerun_pruned/component_features.tsv \
  --protein_features ./golden/bw_pipeline/rerun_pruned/protein_features.tsv \
  --edge_features ./golden/bw_pipeline/rerun_pruned/edge_features.tsv \
  --hgt_candidates ./golden/bw_pipeline/rerun_pruned/hgt_candidates.tsv \
  --out_prefix ./golden/bw_pipeline/results/global_stats

python reporting/explain_component.py \
  --component_id 5 \
  --edges ./golden/bw_pipeline/rerun_pruned/edge_features.tsv \
  --protein_features ./golden/bw_pipeline/rerun_pruned/protein_features.tsv \
  --component_features ./golden/bw_pipeline/rerun_pruned/component_features.tsv \
  --hgt_candidates ./golden/bw_pipeline/rerun_pruned/hgt_candidates.tsv \
  --top_nodes 20 --top_edges 25 \
  > ./golden/bw_pipeline/explanations/comp_5_explained.txt

python reporting/explain_component.py \
  --component_id 8 \
  --edges ./golden/bw_pipeline/rerun_pruned/edge_features.tsv \
  --protein_features ./golden/bw_pipeline/rerun_pruned/protein_features.tsv \
  --component_features ./golden/bw_pipeline/rerun_pruned/component_features.tsv \
  --hgt_candidates ./golden/bw_pipeline/rerun_pruned/hgt_candidates.tsv \
  --top_nodes 20 --top_edges 25 \
  > ./golden/bw_pipeline/explanations/comp_8_explained.txt

python reporting/explain_component.py \
  --component_id 32 \
  --edges ./golden/bw_pipeline/rerun_pruned/edge_features.tsv \
  --protein_features ./golden/bw_pipeline/rerun_pruned/protein_features.tsv \
  --component_features ./golden/bw_pipeline/rerun_pruned/component_features.tsv \
  --hgt_candidates ./golden/bw_pipeline/rerun_pruned/hgt_candidates.tsv \
  --top_nodes 20 --top_edges 25 \
  > ./golden/bw_pipeline/explanations/comp_32_explained.txt

python reporting/explain_top_candidates.py \
  --edges ./golden/bw_pipeline/rerun_pruned/edge_features.tsv \
  --protein_features ./golden/bw_pipeline/rerun_pruned/protein_features.tsv \
  --component_features ./golden/bw_pipeline/rerun_pruned/component_features.tsv \
  --hgt_candidates ./golden/bw_pipeline/rerun_pruned/hgt_candidates.tsv \
  --top_n 20 --top_k_neighbors 12 \
  > ./golden/bw_pipeline/explanations/top_candidates_explained.txt
```

## Preserved Historical References

- `hackathon_report_refs/` (copied from `artifacts/historical_results/results_hackathon_run/`)
- `no_bw_pipeline/legacy_rerun_pruned_reference/` (copied from historical rerun outputs)
