import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from hgt_pipeline import pipeline as hgt


# Keep the original sweep's pruning behavior for backward-comparable results.
def keep_q_percentile_edges(df, q=0.9):
    copy_df = df.copy()
    copy_df["species_pair"] = copy_df.apply(
        lambda r: tuple(sorted([r["species_u"], r["species_v"]])), axis=1
    )
    thresholds = copy_df.groupby("species_pair")["jaccard"].transform(lambda x: x.quantile(q))
    filtered_df = copy_df[copy_df["jaccard"] >= thresholds]
    return filtered_df.drop(columns=["species_pair"])


def keep_top_X_edges_per_node(df, X=20):
    filtered = df.groupby("v").filter(lambda x: len(x) <= X)
    top_edges = filtered.groupby("u").apply(lambda x: x.nlargest(15, "jaccard")).reset_index(drop=True)
    return top_edges


def generate_unpruned_world_at_sim(num_species, proteins_per_species, hgt_sim):
    species_list = [f"Sp_{chr(65+i)}" for i in range(num_species)]
    nodes = {sp: [f"{sp}_P{j}" for j in range(proteins_per_species)] for sp in species_list}
    edges = []

    # Background noise (the decoys)
    for i in range(num_species):
        for j in range(i + 1, num_species):
            sp_u, sp_v = species_list[i], species_list[j]
            for _ in range(250):
                u, v = random.choice(nodes[sp_u]), random.choice(nodes[sp_v])
                jaccard = random.betavariate(1, 15) * 0.3
                edges.append(
                    {
                        "u": u,
                        "v": v,
                        "shared_kmers": int(jaccard * 1000),
                        "jaccard": jaccard,
                        "species_u": sp_u,
                        "species_v": sp_v,
                    }
                )

    # Inject HGTs at the specific target similarity.
    tp_nodes = set()
    for _ in range(30):
        sp_u, sp_v = random.sample(species_list, 2)
        u, v = random.choice(nodes[sp_u]), random.choice(nodes[sp_v])
        tp_nodes.update([u, v])
        edges.append(
            {
                "u": u,
                "v": v,
                "shared_kmers": int(hgt_sim * 1000),
                "jaccard": hgt_sim,
                "species_u": sp_u,
                "species_v": sp_v,
            }
        )

    return pd.DataFrame(edges), tp_nodes


def evaluate_pipeline_performance(hgt_sim):
    """Run one full cycle and return AUC."""
    df_raw, tp_nodes = generate_unpruned_world_at_sim(10, 50, hgt_sim)

    # Apply filtering.
    df_q = keep_q_percentile_edges(df_raw, q=0.9)
    df_pruned = keep_top_X_edges_per_node(df_q, X=20)

    if df_pruned.empty:
        return 0.0

    pruned_edges = [hgt.Edge(**row) for row in df_pruned.to_dict("records")]
    pruned_edges = hgt.dedupe_edges(pruned_edges)

    # Pipeline logic.
    pair_stats = hgt.compute_pair_robust_stats(pruned_edges, weight="jaccard", min_pair_edges=5)
    edge_feats = hgt.compute_edge_features(pruned_edges, pair_stats, min_pair_edges_for_z=5)
    G = hgt.build_graph(pruned_edges, use_networkx=True)
    hgt.attach_edge_z_to_graph(G, edge_feats)
    cid_of, comps = hgt.compute_components(G)
    comp_feats = hgt.compute_component_features(G, comps, z0=3.0)
    node_feats = hgt.compute_node_features(G, cid_of, compute_clustering=True, compute_betweenness=True)

    ranked = hgt.score_hgt_likeness(node_feats, comp_feats, density_eta=8.0, cluster_eta=4.0)

    score_dict = {u: s for u, s in ranked}
    y_true, y_scores = [], []
    for node in list(G.nodes):
        y_true.append(1 if node in tp_nodes else 0)
        y_scores.append(score_dict.get(node, 0.0))

    if len(set(y_true)) < 2:
        return 0.0  # No TPs survived pruning.

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    return auc(fpr, tpr)


if __name__ == "__main__":
    sim_range = np.linspace(0.1, 0.95, 15)
    auc_results = []

    print("Running sensitivity sweep across Jaccard similarities...")
    for sim in sim_range:
        current_auc = evaluate_pipeline_performance(sim)
        auc_results.append(current_auc)
        print(f"Jaccard: {sim:.2f} | AUC: {current_auc:.4f}")

    plt.figure(figsize=(10, 6))
    plt.plot(sim_range, auc_results, marker="o", linestyle="-", color="darkred", linewidth=2)
    plt.xlabel("HGT Jaccard Similarity (Signal Strength)", fontsize=12)
    plt.ylabel("Pipeline ROC-AUC", fontsize=12)
    plt.title("ROC-AUC vs. HGT Signal Strength", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()
