#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import os
import argparse
import sys

def run_clustering(df: pd.DataFrame):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df[['fluidity_calc', 'openness']])

    results = {}
    for k in [3, 4, 5]:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X_scaled)
        score = silhouette_score(X_scaled, labels)
        centroids_original = scaler.inverse_transform(kmeans.cluster_centers_)
        results[k] = {'labels': labels, 'silhouette': score, 'centroids': centroids_original}

    best_k = max(results, key=lambda k: results[k]['silhouette'])
    best_labels = results[best_k]['labels']
    df['cluster'] = best_labels

    centroids_df = pd.DataFrame(results[best_k]['centroids'],
                                columns=['fluidity_calc', 'openness'])

    cluster_to_category = {}
    open_cluster = centroids_df.assign(idx=centroids_df.index).sort_values(
        by=['fluidity_calc','openness'], ascending=[False, True]).iloc[0]['idx']
    cluster_to_category[open_cluster] = 'Open'

    closed_cluster = centroids_df.assign(idx=centroids_df.index).sort_values(
        by=['fluidity_calc','openness'], ascending=[True, False]).iloc[0]['idx']
    cluster_to_category[closed_cluster] = 'Closed'

    for c in centroids_df.index:
        if c not in cluster_to_category:
            cluster_to_category[c] = 'Moderate' if best_k >= 3 else 'Other'

    df['category'] = df['cluster'].map(cluster_to_category)
    return df, centroids_df, best_k

import matplotlib.pyplot as plt

def plot_clusters(df, centroids_df, outdir):
    """
    Scatter plot of fluidity vs openness with cluster assignments.
    Centroids are marked with larger black Xs.
    """

    plt.figure(figsize=(8,6))

    # Color by category
    categories = df['category'].unique()
    colors = {'Open':'tab:blue', 'Closed':'tab:red',
              'Moderate':'tab:green', 'Other':'tab:gray'}

    for cat in categories:
        subset = df[df['category'] == cat]
        plt.scatter(subset['fluidity_calc'], subset['openness'],
                    c=colors.get(cat,'tab:gray'), label=cat, alpha=0.7, s=40)

    # Plot centroids
    plt.scatter(centroids_df['fluidity_calc'], centroids_df['openness'],
                c='black', marker='x', s=200, linewidths=3, label='Centroids')

    plt.xlabel("Genomic fluidity")
    plt.ylabel("Pangenome openness (a)")
    plt.title("Cluster assignments of pangenomes")
    plt.legend()
    plt.tight_layout()

    fig_path = os.path.join(outdir, "cluster_assignments.pdf")
    plt.savefig(fig_path)
    plt.close()

    print(f"Cluster figure saved to {fig_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge fluidity/openness and cluster pangenomes")
    parser.add_argument("--matched", required=True, help="Path to matched_all.csv")
    parser.add_argument("--alpha", required=True, help="Path to alpha_fluidity_all.csv")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    matched_df = pd.read_csv(args.matched)

    if {'fluidity_calc', 'openness'}.issubset(matched_df.columns):
        print(f"{args.matched} already contains fluidity_calc and openness columns. "
              "categorise.py should not be run again. Exiting.")
        sys.exit(0)

    alpha_df = pd.read_csv(args.alpha)
    alpha_df = alpha_df.rename(columns={"fluidity": "fluidity_calc"})

    matched_df['species_taxid'] = matched_df['species_taxid'].astype(str)
    alpha_df['species_taxid'] = alpha_df['species_taxid'].astype(str)

    merged = matched_df.merge(alpha_df[['species_taxid','fluidity_calc','openness']],
                              on='species_taxid', how='left')

    # Reorder fluidity_calc/openness after pangenome_fluidity
    cols = list(merged.columns)
    for c in ['fluidity_calc','openness']:
        if c in cols:
            cols.remove(c)
    insert_at = cols.index("pangenome_fluidity") + 1
    cols = cols[:insert_at] + ['fluidity_calc','openness'] + cols[insert_at:]
    merged = merged[cols]

    # Run clustering
    clustered, centroids_df, best_k = run_clustering(merged.copy())
    plot_clusters(clustered, centroids_df, args.outdir)

    # Insert category as first column
    merged['category'] = clustered['category']
    merged = merged[['category'] + [c for c in merged.columns if c != 'category']]

    os.makedirs(args.outdir, exist_ok=True)
    out_csv = os.path.join(args.outdir, "matched_all.csv")
    merged.to_csv(out_csv, index=False)

    # Uncomment if for extra files
    #clustered[['species_taxid','fluidity_calc','openness','category']].to_csv(
    #    os.path.join(args.outdir, "species_category_assignments.csv"), index=False)
    #centroids_df.to_csv(os.path.join(args.outdir, "cluster_centroids.csv"), index_label='cluster')

    print(f"Best k based on silhouette score: {best_k}")
    print(f"Updated matched_all.csv written with {len(merged)} rows.")
