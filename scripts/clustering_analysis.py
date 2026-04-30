"""
Phase 2: Cross-Factor Effect Matrix & Clustering
=================================================
Build the locus x factor BETA matrix from CDG2025 data,
then run hierarchical, k-means, GMM, and DBSCAN clustering.

HOW TO USE:
  - Uncomment one step at a time, top to bottom
  - Run the full script after each uncomment
  - Each step prints its results and saves figures
  - Later steps depend on earlier ones, so keep them uncommented

Usage: python scripts/clustering_analysis.py
"""
import sys
sys.path.insert(0, ".")

from src.utils import setup_output, teardown_output
output_path = setup_output("clustering_output.txt")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

from src.config import CROSS_DISORDER
from src.data_utils import load_sumstats

RESULTS_DIR = Path("data/results")
FIGURES_DIR = Path("figures")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# ============================================
# STEP 1: Load hits file and extract unique SNPs
# ============================================
print("="*60)
print("  Step 1: Loading hits file")
print("="*60)

df_hits = load_sumstats(CROSS_DISORDER["cdg2025_hits"]["path"])
print(f"  Total loci: {len(df_hits)}")
print(f"  Columns: {df_hits.columns.tolist()}")

# The same SNP can appear in multiple factor hits (e.g. significant
# for both F2 and Hier). We want unique SNPs for the lookup.
hit_snps = df_hits["SNP"].unique()
print(f"  Unique lead SNPs: {len(hit_snps)}")

# Quick look at how many duplicates there are:
print(f"  Total rows (719) vs unique SNPs ({len(hit_snps)}) = {719 - len(hit_snps)} duplicated across factors")


# ============================================
# STEP 2: Build cross-factor effect matrix
# ============================================
# For each factor (F1-PFactor), load the full file (~2.8M SNPs),
# filter to our hit SNPs, and extract the BETA value.
# End result: a DataFrame with rows=SNPs, columns=factors

print("\n" + "="*60)
print("  Step 2: Building cross-factor effect matrix")
print("="*60)

FACTOR_KEYS = {
    "F1_Compulsive":     "cdg2025_F1_compulsive",
    "F2_SCZ_BIP":        "cdg2025_F2_scz_bip",
    "F3_Neurodev":       "cdg2025_F3_neurodev",
    "F4_Internalizing":  "cdg2025_F4_internalizing",
    "F5_Substance":      "cdg2025_F5_substance",
    "PFactor":           "cdg2025_pfactor",
}

factor_betas = {}

for factor_name, config_key in FACTOR_KEYS.items():
    print(f"\n  Loading {factor_name}...")
    
    df_factor = load_sumstats(CROSS_DISORDER[config_key]["path"])

    # filter to only hit snps
    filtered = df_factor[df_factor["SNP"].isin(hit_snps)]

    # Set SNP as index and store the BETA column
    factor_betas[factor_name] = factor_betas[factor_name].set_index("SNP")["BETA"]

    print(f"    Found {len(filtered)} / {len(hit_snps)} hit SNPs")
    del df_factor

# Combine into a single DataFrame
effect_matrix = pd.DataFrame(factor_betas)

print(f"\n  Effect matrix shape: {effect_matrix.shape}")
print(f"  Expected: ({len(hit_snps)}, 5)")
print(f"  Missing values per column:\n{effect_matrix.isnull().sum()}")

# VERIFICATION: pick a random SNP and manually check one value
sample_snp = effect_matrix.index[0]
print(f"\n  Verification — {sample_snp}:")
print(f"    Effect matrix F1 BETA: {effect_matrix.loc[sample_snp, 'F1_Compulsive']}")
print(f"    (Check this against the F1 factor file manually)")


# # ============================================
# # STEP 3: Sanity checks on the effect matrix
# # ============================================
# print("\n" + "="*60)
# print("  Step 3: Effect matrix sanity checks")
# print("="*60)
#
# print("\n  Descriptive statistics:")
# print(effect_matrix.describe().round(6).to_string())
#
# print("\n  Inter-factor correlations:")
# print(effect_matrix.corr().round(4).to_string())
# # If any pair is >0.8, they carry redundant info for clustering.
# # Genomic SEM factors should be low-to-moderate correlation.
#
# # Save a quick correlation heatmap
# fig, ax = plt.subplots(figsize=(8, 6))
# sns.heatmap(effect_matrix.corr(), annot=True, fmt=".3f",
#             cmap="RdBu_r", center=0, ax=ax, vmin=-1, vmax=1)
# ax.set_title("Inter-Factor Correlation (Hit Loci BETAs)")
# fig.savefig(FIGURES_DIR / "factor_correlation.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"\n  Saved: {FIGURES_DIR / 'factor_correlation.png'}")


# # ============================================
# # STEP 4: Standardize the effect matrix
# # ============================================
# # Standardize so factors with larger BETAs don't dominate
# # distance calculations in clustering.
#
# print("\n" + "="*60)
# print("  Step 4: Standardization")
# print("="*60)
#
# from sklearn.preprocessing import StandardScaler
#
# # YOUR CODE: Fit and transform the effect matrix
# # Hint: scaler = StandardScaler()
# #       X_scaled = scaler.fit_transform(effect_matrix)
# scaler = ???
# X_scaled = ???
#
# print(f"  Scaled matrix shape: {X_scaled.shape}")
# print(f"  Column means (should be ~0): {X_scaled.mean(axis=0).round(6)}")
# print(f"  Column stds (should be ~1):  {X_scaled.std(axis=0).round(6)}")


# # ============================================
# # STEP 5: Hierarchical clustering
# # ============================================
# # Exploratory — the dendrogram shows natural groupings
# # without forcing a specific k. Look for large vertical gaps.
#
# print("\n" + "="*60)
# print("  Step 5: Hierarchical clustering")
# print("="*60)
#
# from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
# from scipy.spatial.distance import pdist
#
# # YOUR CODE: Compute pairwise distances and linkage
# # Hint: distances = pdist(X_scaled, metric="euclidean")
# #       Z_linkage = linkage(distances, method="ward")
# distances = ???
# Z_linkage = ???
#
# # Plot the dendrogram
# fig, ax = plt.subplots(figsize=(14, 6))
# dendrogram(Z_linkage, ax=ax, truncate_mode="lastp", p=30,
#            leaf_rotation=90, leaf_font_size=8)
# ax.set_title("Hierarchical Clustering Dendrogram (Ward, Euclidean)")
# ax.set_xlabel("Cluster (size)")
# ax.set_ylabel("Distance")
# fig.savefig(FIGURES_DIR / "dendrogram_ward.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"  Saved: {FIGURES_DIR / 'dendrogram_ward.png'}")
#
# # Also try cutting the dendrogram at a few heights to see cluster counts
# for n_clust in [2, 3, 4, 5, 6]:
#     labels_hier = fcluster(Z_linkage, n_clust, criterion="maxclust")
#     sizes = pd.Series(labels_hier).value_counts().sort_index()
#     print(f"    k={n_clust}: cluster sizes = {sizes.tolist()}")


# # ============================================
# # STEP 6: K-Means with silhouette analysis
# # ============================================
# # Try k=2 through 10. Silhouette score measures how well-separated
# # clusters are. Range: -1 to 1, higher = better.
#
# print("\n" + "="*60)
# print("  Step 6: K-Means clustering")
# print("="*60)
#
# from sklearn.cluster import KMeans
# from sklearn.metrics import silhouette_score
#
# k_range = range(2, 11)
# kmeans_results = {}
#
# # YOUR CODE: Loop through k values, fit KMeans, compute silhouette
# # Hint: for k in k_range:
# #           km = KMeans(n_clusters=k, n_init=50, random_state=42)
# #           labels = km.fit_predict(X_scaled)
# #           sil = silhouette_score(X_scaled, labels)
# #           kmeans_results[k] = {"model": km, "labels": labels, "silhouette": sil}
# #           print(f"    k={k}: silhouette={sil:.4f}")
# for k in k_range:
#     pass  # replace with your code
#
# # Plot silhouette scores
# fig, ax = plt.subplots(figsize=(8, 5))
# sil_scores = [kmeans_results[k]["silhouette"] for k in k_range]
# ax.plot(list(k_range), sil_scores, "bo-", linewidth=2, markersize=8)
# ax.set_xlabel("Number of clusters (k)")
# ax.set_ylabel("Silhouette Score")
# ax.set_title("K-Means: Silhouette Score vs k")
# ax.set_xticks(list(k_range))
# # Highlight the best k
# best_sil_k = list(k_range)[np.argmax(sil_scores)]
# ax.axvline(x=best_sil_k, color="red", linestyle="--", alpha=0.5,
#            label=f"Best k={best_sil_k}")
# ax.legend()
# fig.savefig(FIGURES_DIR / "kmeans_silhouette.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"\n  Best silhouette k={best_sil_k} (score={max(sil_scores):.4f})")
# print(f"  Saved: {FIGURES_DIR / 'kmeans_silhouette.png'}")


# # ============================================
# # STEP 7: Gaussian Mixture Models
# # ============================================
# # Soft clustering — each point gets a probability per cluster.
# # BIC (lower = better) helps pick k. Also compute silhouette
# # on hard assignments for comparison with k-means.
#
# print("\n" + "="*60)
# print("  Step 7: Gaussian Mixture Models")
# print("="*60)
#
# from sklearn.mixture import GaussianMixture
#
# gmm_results = {}
#
# # YOUR CODE: Loop through k values, fit GMM, compute BIC + silhouette
# # Hint: for k in k_range:
# #           gmm = GaussianMixture(n_components=k, n_init=10,
# #                                  covariance_type="full", random_state=42)
# #           gmm.fit(X_scaled)
# #           labels = gmm.predict(X_scaled)
# #           bic = gmm.bic(X_scaled)
# #           sil = silhouette_score(X_scaled, labels)
# #           gmm_results[k] = {"model": gmm, "labels": labels,
# #                              "bic": bic, "silhouette": sil}
# #           print(f"    k={k}: BIC={bic:.1f}, silhouette={sil:.4f}")
# for k in k_range:
#     pass  # replace with your code
#
# # Plot BIC and silhouette side by side
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
#
# bics = [gmm_results[k]["bic"] for k in k_range]
# ax1.plot(list(k_range), bics, "go-", linewidth=2, markersize=8)
# ax1.set_xlabel("Number of components (k)")
# ax1.set_ylabel("BIC (lower = better)")
# ax1.set_title("GMM: BIC vs k")
# ax1.set_xticks(list(k_range))
#
# gmm_sils = [gmm_results[k]["silhouette"] for k in k_range]
# ax2.plot(list(k_range), gmm_sils, "ro-", linewidth=2, markersize=8)
# ax2.set_xlabel("Number of components (k)")
# ax2.set_ylabel("Silhouette Score")
# ax2.set_title("GMM: Silhouette vs k")
# ax2.set_xticks(list(k_range))
#
# fig.tight_layout()
# fig.savefig(FIGURES_DIR / "gmm_bic_silhouette.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"\n  Best BIC k={list(k_range)[np.argmin(bics)]}")
# print(f"  Best silhouette k={list(k_range)[np.argmax(gmm_sils)]}")
# print(f"  Saved: {FIGURES_DIR / 'gmm_bic_silhouette.png'}")


# # ============================================
# # STEP 8: DBSCAN (density-based sanity check)
# # ============================================
# # DBSCAN finds clusters by density — no k needed.
# # Points in sparse regions = "noise" (label -1).
# # If DBSCAN finds similar structure to k-means/GMM, clusters are real.
#
# print("\n" + "="*60)
# print("  Step 8: DBSCAN")
# print("="*60)
#
# from sklearn.cluster import DBSCAN
# from sklearn.neighbors import NearestNeighbors
#
# # First: k-distance graph to estimate eps
# # YOUR CODE: Fit nearest neighbors and get distances to 5th neighbor
# # Hint: nn = NearestNeighbors(n_neighbors=5)
# #       nn.fit(X_scaled)
# #       distances_nn, _ = nn.kneighbors(X_scaled)
# #       k_distances = np.sort(distances_nn[:, -1])
# nn = ???
# distances_nn, _ = ???
# k_distances = ???
#
# fig, ax = plt.subplots(figsize=(8, 5))
# ax.plot(k_distances)
# ax.set_xlabel("Points (sorted by distance)")
# ax.set_ylabel("5th Nearest Neighbor Distance")
# ax.set_title("K-Distance Graph (for eps selection)")
# fig.savefig(FIGURES_DIR / "dbscan_kdistance.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"  Saved: {FIGURES_DIR / 'dbscan_kdistance.png'}")
# print(f"  Look at the plot — the elbow suggests your eps value\n")
#
# # Try several eps values around where the elbow appears
# # YOUR CODE: Loop through eps values and run DBSCAN
# # Hint: for eps_val in [0.5, 1.0, 1.5, 2.0, 2.5]:
# #           db = DBSCAN(eps=eps_val, min_samples=5)
# #           labels = db.fit_predict(X_scaled)
# #           n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
# #           n_noise = (labels == -1).sum()
# #           print(f"    eps={eps_val}: {n_clusters} clusters, {n_noise} noise points")
# for eps_val in [0.5, 1.0, 1.5, 2.0, 2.5]:
#     pass  # replace with your code


# # ============================================
# # STEP 9: Pick best k and visualize
# # ============================================
# # Look at Steps 5-8 results and pick the k where:
# #   - Silhouette is high (k-means and GMM)
# #   - Dendrogram shows natural breaks
# #   - DBSCAN finds similar number of clusters
# #   - Clusters are interpretable
#
# print("\n" + "="*60)
# print("  Step 9: Visualization with best k")
# print("="*60)
#
# # YOUR DECISION: set this based on your results
# BEST_K = ???
# best_labels = kmeans_results[BEST_K]["labels"]
# print(f"  Chosen k={BEST_K}")
# print(f"  Cluster sizes: {pd.Series(best_labels).value_counts().sort_index().tolist()}")


# # --- 9A: Heatmap sorted by cluster ---
# # The money figure. Distinct horizontal bands = real clusters.
#
# df_plot = effect_matrix.copy()
# df_plot["cluster"] = best_labels
# df_plot = df_plot.sort_values("cluster")
#
# sorted_index = df_plot.index
# X_sorted = X_scaled[effect_matrix.index.get_indexer(sorted_index)]
#
# fig, ax = plt.subplots(figsize=(10, 14))
# sns.heatmap(
#     pd.DataFrame(X_sorted, columns=FACTOR_KEYS.keys()),
#     cmap="RdBu_r", center=0, ax=ax,
#     yticklabels=False,
#     cbar_kws={"label": "Standardized BETA"}
# )
# cluster_sizes = df_plot["cluster"].value_counts().sort_index()
# cumulative = 0
# for size in cluster_sizes.values[:-1]:
#     cumulative += size
#     ax.axhline(y=cumulative, color="black", linewidth=2)
# ax.set_title(f"Cross-Factor Effect Profiles (k={BEST_K}, sorted by cluster)")
# ax.set_ylabel("Loci")
# fig.savefig(FIGURES_DIR / "heatmap_clustered.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"  Saved: {FIGURES_DIR / 'heatmap_clustered.png'}")


# # --- 9B: UMAP projection colored by cluster ---
#
# import umap
#
# reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
# embedding = reducer.fit_transform(X_scaled)
#
# fig, ax = plt.subplots(figsize=(10, 8))
# scatter = ax.scatter(embedding[:, 0], embedding[:, 1],
#                      c=best_labels, cmap="tab10", s=15, alpha=0.7)
# ax.set_xlabel("UMAP 1")
# ax.set_ylabel("UMAP 2")
# ax.set_title(f"UMAP Projection of Pleiotropic Loci (k={BEST_K})")
# plt.colorbar(scatter, label="Cluster")
# fig.savefig(FIGURES_DIR / "umap_clusters.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"  Saved: {FIGURES_DIR / 'umap_clusters.png'}")


# # --- 9C: Cluster profile bar charts ---
# # Mean BETA per factor for each cluster — shows what each
# # cluster "looks like" biologically.
#
# df_profiles = effect_matrix.copy()
# df_profiles["cluster"] = best_labels
# cluster_means = df_profiles.groupby("cluster").mean()
# print("\n  Cluster mean BETAs:")
# print(cluster_means.round(4).to_string())
#
# fig, axes = plt.subplots(1, BEST_K, figsize=(4*BEST_K, 5), sharey=True)
# if BEST_K == 1:
#     axes = [axes]
# for i, ax in enumerate(axes):
#     means = cluster_means.loc[i]
#     colors = ["red" if v < 0 else "steelblue" for v in means]
#     ax.bar(means.index, means.values, color=colors)
#     ax.set_title(f"Cluster {i} (n={sum(best_labels == i)})")
#     ax.axhline(y=0, color="black", linewidth=0.5)
#     ax.tick_params(axis="x", rotation=45)
# fig.suptitle("Mean Cross-Factor Effect Profiles by Cluster")
# fig.tight_layout()
# fig.savefig(FIGURES_DIR / "cluster_profiles.png", dpi=300, bbox_inches="tight")
# plt.close()
# print(f"  Saved: {FIGURES_DIR / 'cluster_profiles.png'}")


# # ============================================
# # STEP 10: Save results
# # ============================================
# print("\n" + "="*60)
# print("  Step 10: Saving results")
# print("="*60)
#
# df_output = effect_matrix.copy()
# df_output["cluster_kmeans"] = kmeans_results[BEST_K]["labels"]
# df_output["cluster_gmm"] = gmm_results[BEST_K]["labels"]
# df_output.to_csv(RESULTS_DIR / "effect_matrix_clustered.csv")
# print(f"  Saved: {RESULTS_DIR / 'effect_matrix_clustered.csv'}")
#
# # Save summary text
# with open(RESULTS_DIR / "clustering_summary.txt", "w") as f:
#     f.write("Phase 2: Clustering Analysis Summary\n")
#     f.write("="*60 + "\n\n")
#     f.write(f"Hit loci: {len(hit_snps)} unique SNPs\n")
#     f.write(f"Effect matrix: {effect_matrix.shape}\n")
#     f.write(f"Chosen k: {BEST_K}\n\n")
#     f.write("K-Means Silhouette Scores:\n")
#     for k in k_range:
#         f.write(f"  k={k}: {kmeans_results[k]['silhouette']:.4f}\n")
#     f.write(f"\nGMM BIC Values:\n")
#     for k in k_range:
#         f.write(f"  k={k}: {gmm_results[k]['bic']:.1f}\n")
#     f.write(f"\nCluster sizes (k-means, k={BEST_K}):\n")
#     for i in range(BEST_K):
#         n = sum(kmeans_results[BEST_K]["labels"] == i)
#         f.write(f"  Cluster {i}: {n} loci\n")
#     f.write(f"\nCluster mean BETAs:\n")
#     f.write(cluster_means.round(4).to_string())
# print(f"  Saved: {RESULTS_DIR / 'clustering_summary.txt'}")
#
# print("\n  Phase 2: Clustering Analysis Complete")

if __name__ == "__main__":
    teardown_output(output_path)
    print("Phase 2: Clustering Analysis Complete")
