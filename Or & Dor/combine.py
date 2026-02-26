import pandas as pd
import re

# -------------------------
# 1. Load cluster mapping
# -------------------------
clusters = pd.read_csv(
    "combined_out_cluster.tsv",
    sep="\t",
    header=None,
    names=["cluster_rep", "protein_id"]
)

# Assign numeric cluster IDs
clusters["cluster_id"] = (
    clusters["cluster_rep"]
    .astype("category")
    .cat.codes + 1
)

clusters = clusters[["protein_id", "cluster_id"]]

# -------------------------
# 2. Parse FASTA headers
# -------------------------
records = []

with open("combined.fasta") as f:
    for line in f:
        if line.startswith(">"):
            line = line.strip()

            # Protein ID
            protein_id = line.split()[0][1:]

            # Organism (inside [])
            org_match = re.search(r"\[(.*)\]$", line)
            organism = org_match.group(1) if org_match else None

            # Protein name (between ID and [])
            protein_name = line[len(protein_id) + 2:]
            protein_name = re.sub(r"\s*\[.*\]$", "", protein_name)

            records.append({
                "protein_id": protein_id,
                "protein_name": protein_name,
                "organism": organism
            })

fasta_df = pd.DataFrame(records)

# -------------------------
# 3. Merge + save
# -------------------------
final_df = fasta_df.merge(clusters, on="protein_id", how="left")

final_df.to_csv("merged.csv", index=False)

print(final_df.head())
