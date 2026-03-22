"""Stress test for ancient/ameliorated HGT under strong conserved-hub confounders."""

import random
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from hgt_pipeline import pipeline as hgt
from graph_construction.graph_pruning import keep_q_percentile_edges, keep_top_X_edges_per_node


def generate_difficult_world(num_species=12, proteins_per_species=100):
    print("Generating DIFFICULT world: Low-signal HGTs vs. high-similarity hubs...")
    species_list = [f"Sp_{i}" for i in range(num_species)]
    nodes = {sp: [f"{sp}_P{j}" for j in range(proteins_per_species)] for sp in species_list}
    edges = []

    # 1) Trap: universally conserved hubs (TN).
    num_hubs = 10
    hub_nodes = set()
    for h in range(num_hubs):
        orthologs = [nodes[sp][h] for sp in species_list]
        hub_nodes.update(orthologs)
        for i in range(len(orthologs)):
            for j in range(i + 1, len(orthologs)):
                jacc = random.uniform(0.90, 0.99)
                edges.append(
                    {
                        "u": orthologs[i],
                        "v": orthologs[j],
                        "jaccard": jacc,
                        "species_u": orthologs[i].split("_")[0],
                        "species_v": orthologs[j].split("_")[0],
                    }
                )

    # 2) Challenge: ancient/weak HGTs (TP).
    tp_nodes = set()
    for _ in range(50):
        sp_u, sp_v = random.sample(species_list, 2)
        u = random.choice([n for n in nodes[sp_u] if n not in hub_nodes and n not in tp_nodes])
        v = random.choice([n for n in nodes[sp_v] if n not in hub_nodes and n not in tp_nodes])
        tp_nodes.update([u, v])
        jacc = random.uniform(0.30, 0.60)
        edges.append({"u": u, "v": v, "jaccard": jacc, "species_u": sp_u, "species_v": sp_v})

    # 3) Background noise.
    for i in range(num_species):
        for j in range(i + 1, num_species):
            for _ in range(200):
                u, v = random.choice(nodes[species_list[i]]), random.choice(nodes[species_list[j]])
                if u in hub_nodes or v in hub_nodes or u in tp_nodes or v in tp_nodes:
                    continue
                jacc = random.betavariate(1, 20) * 0.5
                edges.append({"u": u, "v": v, "jaccard": jacc, "species_u": species_list[i], "species_v": species_list[j]})

    return pd.DataFrame(edges), tp_nodes, hub_nodes


def run_stress_test():
    df_raw, tp_nodes, hub_nodes = generate_difficult_world()

    print("Applying upstream pruning...")
    df_filtered = keep_q_percentile_edges(df_raw, q=0.1)
    df_pruned = keep_top_X_edges_per_node(df_filtered, X=15)

    pruned_edges = [
        hgt.Edge(
            u=r.u,
            v=r.v,
            shared_kmers=int(r.jaccard * 1000),
            jaccard=r.jaccard,
            species_u=r.species_u,
            species_v=r.species_v,
        )
        for r in df_pruned.itertuples()
    ]

    G = hgt.build_graph(pruned_edges, use_networkx=True)
    pair_stats = hgt.compute_pair_robust_stats(pruned_edges, weight="jaccard")
    edge_feats = hgt.compute_edge_features(pruned_edges, pair_stats)
    hgt.attach_edge_z_to_graph(G, edge_feats)
    cid_of, comps = hgt.compute_components(G)
    comp_feats = hgt.compute_component_features(G, comps)
    node_feats = hgt.compute_node_features(G, cid_of, compute_clustering=True, compute_betweenness=True)

    ranked = hgt.score_hgt_likeness(node_feats, comp_feats, density_eta=10.0, cluster_eta=5.0)

    score_map = {u: s for u, s in ranked}
    y_true, y_scores = [], []
    test_nodes = list(tp_nodes) + list(hub_nodes)
    for node in test_nodes:
        y_true.append(1 if node in tp_nodes else 0)
        y_scores.append(score_map.get(node, 0.0))

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    auc_val = auc(fpr, tpr)
    print(f"STRESS TEST AUC: {auc_val:.4f}")

    plt.plot(fpr, tpr, label=f"Stress Test (AUC={auc_val:.2f})")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve: Ancient HGTs vs. Conserved Hubs")
    plt.show()


if __name__ == "__main__":
    run_stress_test()
