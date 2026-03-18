# Reproduce End-to-End

This repository provides a single root runner for exact, scriptable E2E recipes:

```powershell
python reproduce.py --help
```

## 1) Deterministic Reproduction From Golden Edges

Recreate pipeline outputs + reports from the canonical pruned graph input:

```powershell
python reproduce.py from-edges `
  --in_edges golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv `
  --out_dir rerun_from_golden `
  --with_reports `
  --with_explanations
```

Optional full run with betweenness:

```powershell
python reproduce.py from-edges `
  --in_edges golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv `
  --out_dir rerun_from_golden_bw `
  --with_betweenness `
  --with_reports `
  --with_explanations
```

## 2) Reproduce From Manifest + Downloaded FASTAs

Construct candidate graph and pruned edges, then run the analysis pipeline:

```powershell
python reproduce.py from-manifest `
  --manifest data/out_refseq/manifest.tsv `
  --downloads_dir data/out_refseq/downloads `
  --work_dir rerun_from_manifest `
  --with_reports `
  --with_explanations
```

## 3) Full Reproduction From `assembly_summary_refseq.txt`

Starting from assembly summary + species list, including downloads:

```powershell
python reproduce.py from-assembly-summary `
  --assembly_summary data/assembly_summary_refseq.txt `
  --species_list config/species.txt `
  --work_dir rerun_full `
  --with_reports `
  --with_explanations
```

This runs:
1. `graph_construction.refseq_fetch_proteins` (manifest + FASTA downloads)
2. `graph_construction.orchestrator construct-edges` (candidates + pruned edges)
3. `hgt_pipeline.pipeline` (feature extraction + scoring)
4. `tools/reporting/*` (report tables and explanations)

## Assembly Summary Source

`assembly_summary_refseq.txt` is currently not tracked in git because it is very large.

Download from NCBI:

```powershell
curl.exe -L -o data/assembly_summary_refseq.txt https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```

