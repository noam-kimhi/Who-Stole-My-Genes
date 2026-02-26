import pandas as pd
import sys

# -------------------------
# 1. Load CSV
# -------------------------
df = pd.read_csv('mergedWithSeqs.csv')

# Lowercase protein names once (efficient)
df['protein_name_lc'] = df['protein_name'].str.lower()

# -------------------------
# 2. Query via CLI
# -------------------------
# Check if an argument was provided; otherwise use a default
if len(sys.argv) > 1:
    # Joins all arguments into one string in case query isn't quoted
    query = " ".join(sys.argv[1:])
else:
    print("Usage: python script_name.py <protein_name>")
    sys.exit(1)

query_lc = query.lower()

# Find rows where query is a substring of protein_name
hits = df[df['protein_name_lc'].str.contains(query_lc, na=False)]

# Get unique cluster IDs
cluster_ids = sorted(hits['cluster_id'].unique())

print(f"Results for: '{query}'")
print("Matching cluster IDs:")
print(cluster_ids)
