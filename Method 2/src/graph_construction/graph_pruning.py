import pandas as pd
import matplotlib.pyplot as plt
import argparse
from typing import Optional


def keep_q_percentile_edges(df, q=0.9):
    """Keep only the top q percent of edges based on jaccard, for each pair (u,v)"""
    copy = df.copy()
    copy["species_pair"] = copy.apply(
        lambda r: tuple(sorted([r["species_u"], r["species_v"]])),
        axis=1
    )

    q = 0.9  # keep top 10% per species pair
    # Compute per-pair quantile threshold
    thresholds = copy.groupby("species_pair")["jaccard"].transform(
        lambda x: x.quantile(q)
    )

    filtered_df = copy[copy["jaccard"] >= thresholds]
    return filtered_df.drop(columns=["species_pair"])


def keep_top_X_edges_per_node(
    df,
    X: int = 20,
    max_in_edges_per_v: Optional[int] = None,
    top_u_edges: Optional[int] = None,
):
    """
    Prune candidate edges with configurable degree controls.

    Parameters:
      X:
        Backward-compatible primary knob.
      max_in_edges_per_v:
        Optional cap for edges per target protein `v` (legacy behavior).
      top_u_edges:
        Optional cap for top edges kept per source protein `u`.
        If omitted, falls back to `X`.
    """
    top_k_u = X if top_u_edges is None else top_u_edges
    top_edges = df

    if max_in_edges_per_v is not None:
        top_edges = top_edges.groupby("v").filter(lambda x: len(x) <= max_in_edges_per_v)

    top_edges = top_edges.groupby("u").apply(lambda x: x.nlargest(top_k_u, "jaccard")).reset_index(drop=True)
    return top_edges


def plots(df_before, df_after):
    # check before and after pruning
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    df_before.groupby('u').size().hist(bins=50, ax=axes[0])
    axes[0].set_xlabel('Number of edges per node u')
    axes[0].set_ylabel('Number of nodes u')
    axes[0].set_title('Distribution of edges per node u before pruning')

    df_after.groupby('v').size().hist(bins=50, ax=axes[1])
    axes[1].set_xlabel('Number of edges per node v')
    axes[1].set_title('Distribution of edges per node v after pruning')
    plt.show()
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Graph pruning")

    parser.add_argument(
        "--path",
        type=str,
        required=True,
        help="Path to taxonomic distances file (txt with columns: species_u,species_v,distance)"
    )

    parser.add_argument(
        "--tax_dist_thr",
        type=int,
        required=False
    )

    parser.add_argument(
        "--jaccard_thr",
        type=int,
        required=False
    )

    parser.add_argument(
        "--plot",
        type=str,
        required=False,
        help="Whether to plot distributions before/after pruning"
    )

    args = parser.parse_args()
    path_to_taxonomic_distances = args.path

    df = pd.read_csv('data/out_refseq/candidates.tsv', sep='\t') # read roee's graph

    filtered_df = keep_q_percentile_edges(df, q=0.9)

    top_edges = keep_top_X_edges_per_node(filtered_df, X=20)

    # add filtering based on taxonomic distance
    top_edges['tax_distance'] = None  # Placeholder for taxonomic distances
    with open(path_to_taxonomic_distances, 'r') as f:
        for line in f.readlines():
            u, v, dist = line.split(',')
            top_edges.loc[(top_edges['species_u'] == u) & (top_edges['species_v'] == v), 'tax_distance'] = int(dist)
            top_edges.loc[(top_edges['species_u'] == v) & (top_edges['species_v'] == u), 'tax_distance'] = int(dist)

    # filter out distant edges with low jaccard
    tax_distance_thr = 4 if args.tax_dist_thr is None else args.tax_dist_thr
    jaccard_thr = 0.1 if args.jaccard_thr is None else args.jaccard_thr
    top_edges =  top_edges[~((top_edges.tax_distance > tax_distance_thr) & (top_edges.jaccard < jaccard_thr))]
    # filter out close edges
    top_edges = top_edges[(top_edges.tax_distance >= 2)]

    # Remove redundant proteins. Should make sure we don't lose any HGTs.
    top_edges = top_edges[top_edges.u.str.contains('W') & top_edges.v.str.contains('W')]
    top_edges.to_csv('jaccard_pruned_edges.tsv', sep='\t', index=False)

    if args.plot:
        plots(df, top_edges)


if __name__ == "__main__":
    main()
