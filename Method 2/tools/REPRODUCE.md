# Reproduce End-to-End (Module CLIs)

This document uses only the two module entrypoints used in the codebase:
- `graph_construction.*` (including the graph-construction orchestrator)
- `hgt_pipeline.pipeline`

No `reproduce.py` wrapper is required.

## Environment

From repository root (PowerShell):

```powershell
python -m pip install -e .
```

## A) Start From `config/species.txt` + `data/assembly_summary_refseq.txt`

Canonical manifest is tracked in the repo at:
- `data/out_refseq/manifest.tsv`

1. Build manifest from assembly summary and species config:

```powershell
python -m graph_construction.refseq_fetch_proteins `
  --assembly_summary data/assembly_summary_refseq.txt `
  --species_list config/species.txt `
  --out_dir tmp/e2e_review/graph_construction `
  --max_assemblies_per_species 2 `
  --require_latest
```

2. Construct graph inputs with the graph-construction orchestrator:

```powershell
python -m graph_construction.orchestrator construct-edges `
  --manifest tmp/e2e_review/graph_construction/manifest.tsv `
  --downloads_dir data/out_refseq/downloads `
  --out_candidates tmp/e2e_review/run/candidates.tsv `
  --out_edges tmp/e2e_review/run/edges_pruned.tsv `
  --k 6 --min_len 50 --max_postings 100 --min_shared 6 --top_m 10 --q 0.9 --top_x 20
```

3. Run the analysis pipeline (no betweenness, then betweenness):

```powershell
python -m hgt_pipeline.pipeline `
  --in_edges tmp/e2e_review/run/edges_pruned.tsv `
  --out_dir tmp/e2e_review/run/pipeline_no_bw `
  --no_betweenness

python -m hgt_pipeline.pipeline `
  --in_edges tmp/e2e_review/run/edges_pruned.tsv `
  --out_dir tmp/e2e_review/run/pipeline_bw
```

## B) Reporting / Explanations

No-betweenness reports:

```powershell
python tools/reporting/top_anomaly_edges.py `
  --edges tmp/e2e_review/run/pipeline_no_bw/edge_features.tsv `
  --top_n 25 `
  --out_dir tmp/e2e_review/run/pipeline_no_bw/results

python tools/reporting/summarize_global_stats.py `
  --component_features tmp/e2e_review/run/pipeline_no_bw/component_features.tsv `
  --protein_features tmp/e2e_review/run/pipeline_no_bw/protein_features.tsv `
  --edge_features tmp/e2e_review/run/pipeline_no_bw/edge_features.tsv `
  --hgt_candidates tmp/e2e_review/run/pipeline_no_bw/hgt_candidates.tsv `
  --out_prefix tmp/e2e_review/run/pipeline_no_bw/results/global_stats
```

Betweenness-on reports + explanations:

```powershell
python tools/reporting/top_anomaly_edges.py `
  --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv `
  --top_n 25 `
  --out_dir tmp/e2e_review/run/pipeline_bw/results

python tools/reporting/summarize_global_stats.py `
  --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv `
  --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv `
  --edge_features tmp/e2e_review/run/pipeline_bw/edge_features.tsv `
  --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv `
  --out_prefix tmp/e2e_review/run/pipeline_bw/results/global_stats
```

```powershell
python tools/reporting/explain_component.py --component_id 5 --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_component.py --component_id 8 --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_component.py --component_id 32 --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_top_candidates.py --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_n 20 --top_k_neighbors 12
```

## C) Exact-Match Targets

Expected exact targets:
- `golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv`
- `golden/no_bw_pipeline/rerun_pruned/*.tsv`
- `golden/no_bw_pipeline/results/*.tsv`
- `golden/bw_pipeline/rerun_pruned/*.tsv`
- `golden/bw_pipeline/results/*.tsv`
- `golden/bw_pipeline/explanations/*.txt`

Note:
- `data/assembly_summary_refseq.txt` is large and not tracked by git.
- Download from NCBI when needed:

```powershell
curl.exe -L -o data/assembly_summary_refseq.txt https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```
