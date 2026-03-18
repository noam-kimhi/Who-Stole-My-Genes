# Regression Tests

Run from method root:

```powershell
python tests/test_regression_baselines.py
```

Modes:
- `graph` (default): no-betweenness pipeline regression + fast checks
- `graph_bw`: includes betweenness-on regression
- `full`: includes heavy full k-mer regression (only if full local data exists)

Examples:

```powershell
python tests/test_regression_baselines.py --mode graph
python tests/test_regression_baselines.py --mode graph_bw
python tests/test_regression_baselines.py --mode full
```

