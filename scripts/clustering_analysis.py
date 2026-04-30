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
import gc

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

# TODO: Load the hits file using load_sumstats
# Hint: CROSS_DISORDER["cdg2025_hits"]["path"]
df_hits = load_sumstats(CROSS_DISORDER["cdg2025_hits"]["path"])

print(f"  Total loci: {len(df_hits)}")
print(f"  Columns: {df_hits.columns.tolist()}")

# TODO: Extract the unique lead SNP rsIDs from the hits file
# Hint: df_hits["SNP"].unique() — but note that the same SNP
# can appear multiple times if it's significant for multiple factors.
# We want unique SNPs because we'll look each one up once in the factor files.
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
    filtered = df_factor[df_factor["SNP"].isin(hit_snps)]
    filtered.set_index("SNP")
    factor_betas[factor_name] = filtered["BETA"]
    del df_factor
    gc.collect()

# TODO: Combine into effect matrix
# Hint: effect_matrix = pd.DataFrame(factor_betas)
effect_matrix = pd.DataFrame(factor_betas)


print(f"\n  Effect matrix shape: {effect_matrix.shape}")
print(f"  Expected: ({len(hit_snps)}, 6)")

# TODO: Check for missing values — if any SNPs weren't found in a factor file
# Hint: effect_matrix.isnull().sum()
# If there are NaNs, investigate why. Likely won't be any since all factor
# files share the same SNP set.
effect_matrix.isnull().sum()

# ============================================
# STEP 3: Quick sanity checks on the effect matrix
# ============================================
print("\n" + "="*60)
print("  Step 3: Effect matrix sanity checks")
print("="*60)

# TODO: Print basic statistics for each factor column
# Hint: effect_matrix.describe()
effect_matrix.describe()

# TODO: Check the correlation between factors
# Hint: effect_matrix.corr()
# Question to think about: should the factors be correlated?
# Genomic SEM extracts orthogonal-ish factors, but they're not
# perfectly uncorrelated. High correlation between two factor
# columns would mean they carry redundant information for clustering.
effect_matrix.corr()


# ============================================
# STEP 4: Standardize the effect matrix
# ============================================
# Before clustering, we should standardize so that factors with
# larger BETAs don't dominate the distance calculations.

print("\n" + "="*60)
print("  Step 4: Standardization")
print("="*60)

    # filter to only hit snps
    filtered = df_factor[df_factor["SNP"].isin(hit_snps)]

    # Set SNP as index and store the BETA column
    factor_betas[factor_name] = factor_betas[factor_name].set_index("SNP")["BETA"]

scaler = StandardScaler()
X_scaled = scaler.fit_transform(effect_matrix)
em_indices = effect_matrix.index

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


# TODO: Compute the linkage matrix and plot a dendrogram
# Hint:
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

# Compute pairwise distances
distances = pdist(X_scaled, metric="euclidean")

# Compute linkage (try "ward" — minimizes within-cluster variance)
Z_linkage = linkage(distances, method="ward")

# Plot dendrogram
fig, ax = plt.subplots(figsize=(14, 6))
dendrogram(Z_linkage, ax=ax, truncate_mode="lastp", p=30,
            leaf_rotation=90, leaf_font_size=8)
ax.set_title("Hierarchical Clustering Dendrogram (Ward, Euclidean)")
ax.set_xlabel("Cluster (size)")
ax.set_ylabel("Distance")
fig.savefig(FIGURES_DIR / "dendrogram_ward.png", dpi=300, bbox_inches="tight")
#
# The dendrogram will suggest natural cluster counts — look for
# large vertical gaps. These indicate where merging clusters
# requires a big jump in distance (natural breakpoints).


# ============================================
# STEP 6: K-Means Clustering with Silhouette Analysis
# ============================================
# Try k = 2 through 10. For each k, compute silhouette score.
# The silhouette score measures how similar each point is to its
# own cluster vs. the nearest neighboring cluster. Range: -1 to 1.
# Higher = better defined clusters.

print("\n" + "="*60)
print("  Step 6: K-Means clustering")
print("="*60)

# TODO: Run k-means for k=2 to k=10, store silhouette scores
# Hint:
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

k_range = range(2, 11)
kmeans_results = {}  # k -> {"model": fitted_model, "silhouette": score}

for k in k_range:
    km = KMeans(n_clusters=k, n_init=50, random_state=42)
    labels = km.fit_predict(X_scaled)
    sil = silhouette_score(X_scaled, labels)
    kmeans_results[k] = {"model": km, "labels": labels, "silhouette": sil}
    print(f"    k={k}: silhouette={sil:.4f}")

# n_init=50 runs 50 random initializations and picks the best.
# This matters because k-means is sensitive to initialization.

# your code here

# TODO: Plot silhouette scores vs k (elbow plot)
# Hint:
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(list(k_range), [kmeans_results[k]["silhouette"] for k in k_range],
        "bo-", linewidth=2)
ax.set_xlabel("Number of clusters (k)")
ax.set_ylabel("Silhouette Score")
ax.set_title("K-Means: Silhouette Score vs k")
fig.savefig(FIGURES_DIR / "kmeans_silhouette.png", dpi=300, bbox_inches="tight")


# ============================================
# STEP 7: Gaussian Mixture Models
# ============================================
# GMM is like a soft version of k-means — each point gets a
# probability of belonging to each cluster, not a hard assignment.
# We evaluate with BIC (Bayesian Information Criterion) — lower is better.
# Also compute silhouette on the hard assignments for comparison.


# TODO: Run GMM for k=2 to k=10, store BIC and silhouette
# Hint:
from sklearn.mixture import GaussianMixture

gmm_results = {}
for k in k_range:
    gmm = GaussianMixture(n_components=k, n_init=10,
                            covariance_type="full", random_state=42)
    gmm.fit(X_scaled)
    labels = gmm.predict(X_scaled)
    bic = gmm.bic(X_scaled)
    sil = silhouette_score(X_scaled, labels)
    gmm_results[k] = {"model": gmm, "labels": labels,
                        "bic": bic, "silhouette": sil}
    print(f"    k={k}: BIC={bic:.1f}, silhouette={sil:.4f}")


# TODO: Plot BIC vs k (look for the "elbow" — where BIC stops dropping fast)
# Hint: same pattern as silhouette plot but with BIC on y-axis


# ============================================
# STEP 8: DBSCAN (density-based sanity check)
# ============================================
# DBSCAN doesn't require specifying k — it finds clusters based
# on density. Points in sparse regions become "noise" (label -1).
# This is our reality check: if DBSCAN finds similar structure
# to k-means/GMM, our clusters are likely real.
#
# The tricky part is choosing eps (neighborhood radius).

print("\n" + "="*60)
print("  Step 8: DBSCAN")
print("="*60)

# TODO: Run DBSCAN with a few eps values
# Hint:
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

# First, use k-nearest neighbors to estimate a good eps
# Plot the k-distance graph (sorted distances to 5th nearest neighbor)
nn = NearestNeighbors(n_neighbors=5)
nn.fit(X_scaled)
distances, _ = nn.kneighbors(X_scaled)
k_distances = np.sort(distances[:, -1])  # distance to 5th neighbor

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(k_distances)
ax.set_xlabel("Points (sorted)")
ax.set_ylabel("5th Nearest Neighbor Distance")
ax.set_title("K-Distance Graph (for eps selection)")
fig.savefig(FIGURES_DIR / "dbscan_kdistance.png", dpi=300, bbox_inches="tight")

# Look for the "elbow" in this plot — that's your eps
# Then try a few values around that elbow:
for eps_val in [0.5, 1.0, 1.5, 2.0]:
    db = DBSCAN(eps=eps_val, min_samples=5)
    labels = db.fit_predict(X_scaled)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = (labels == -1).sum()
    print(f"    eps={eps_val}: {n_clusters} clusters, {n_noise} noise points")

# ============================================
# STEP 9: Pick best k and visualize
# ============================================
# Based on Steps 5-8, pick the k that has:
#   - High silhouette score
#   - Consistent across methods (k-means and GMM agree)
#   - Supported by the dendrogram's natural breaks
#   - DBSCAN finds similar number of clusters

print("\n" + "="*60)
print("  Step 9: Visualization with best k")
print("="*60)

# TODO: Set your chosen k based on the results above
# BEST_K = # your choice here — look at the silhouette plots

# TODO: Get the cluster labels from k-means at BEST_K
# Hint: best_labels = kmeans_results[BEST_K]["labels"]
# best_labels = # your code here

# --- 9A: Heatmap of effect matrix sorted by cluster ---
# This is the money figure. It shows how loci within each cluster
# share similar cross-factor effect profiles.
#
# TODO:
#   # Add cluster labels to the effect matrix
#   df_plot = effect_matrix.copy()
#   df_plot["cluster"] = best_labels
#   df_plot = df_plot.sort_values("cluster")
#
#   # Create a heatmap of the scaled values, sorted by cluster
#   # Use the UNSCALED effect_matrix values but sorted by cluster order
#   sorted_index = df_plot.index
#   X_sorted = X_scaled[effect_matrix.index.get_indexer(sorted_index)]
#
#   fig, ax = plt.subplots(figsize=(10, 14))
#   sns.heatmap(
#       pd.DataFrame(X_sorted, columns=FACTOR_KEYS.keys()),
#       cmap="RdBu_r", center=0, ax=ax,
#       yticklabels=False,  # too many loci to label individually
#       cbar_kws={"label": "Standardized BETA"}
#   )
#   # Add horizontal lines between clusters
#   cluster_sizes = df_plot["cluster"].value_counts().sort_index()
#   cumulative = 0
#   for size in cluster_sizes.values[:-1]:
#       cumulative += size
#       ax.axhline(y=cumulative, color="black", linewidth=2)
#   ax.set_title(f"Cross-Factor Effect Profiles (k={BEST_K}, sorted by cluster)")
#   ax.set_ylabel("Loci")
#   fig.savefig(FIGURES_DIR / "heatmap_clustered.png", dpi=300, bbox_inches="tight")

# your code here


# --- 9B: UMAP projection colored by cluster ---
# UMAP reduces 5 dimensions to 2 for visualization.
# Each point = one locus, colored by cluster assignment.
#
# TODO:
#   import umap  # this is umap-learn from our conda env
#
#   reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
#   embedding = reducer.fit_transform(X_scaled)
#
#   fig, ax = plt.subplots(figsize=(10, 8))
#   scatter = ax.scatter(embedding[:, 0], embedding[:, 1],
#                        c=best_labels, cmap="tab10", s=15, alpha=0.7)
#   ax.set_xlabel("UMAP 1")
#   ax.set_ylabel("UMAP 2")
#   ax.set_title(f"UMAP Projection of Pleiotropic Loci (k={BEST_K})")
#   plt.colorbar(scatter, label="Cluster")
#   fig.savefig(FIGURES_DIR / "umap_clusters.png", dpi=300, bbox_inches="tight")

# your code here


# --- 9C: Cluster profile summary ---
# For each cluster, show the mean BETA across factors.
# This tells us what each cluster "looks like" biologically.
#
# TODO:
#   df_profiles = effect_matrix.copy()
#   df_profiles["cluster"] = best_labels
#   cluster_means = df_profiles.groupby("cluster").mean()
#   print("\n  Cluster mean BETAs:")
#   print(cluster_means.round(4).to_string())
#
#   # Bar plot of cluster profiles
#   fig, axes = plt.subplots(1, BEST_K, figsize=(4*BEST_K, 5), sharey=True)
#   if BEST_K == 1:
#       axes = [axes]
#   for i, ax in enumerate(axes):
#       means = cluster_means.loc[i]
#       colors = ["red" if v < 0 else "steelblue" for v in means]
#       ax.bar(means.index, means.values, color=colors)
#       ax.set_title(f"Cluster {i} (n={sum(best_labels == i)})")
#       ax.axhline(y=0, color="black", linewidth=0.5)
#       ax.tick_params(axis="x", rotation=45)
#   fig.suptitle("Mean Cross-Factor Effect Profiles by Cluster")
#   fig.tight_layout()
#   fig.savefig(FIGURES_DIR / "cluster_profiles.png", dpi=300, bbox_inches="tight")

# your code here


# ============================================
# STEP 10: Save results
# ============================================
print("\n" + "="*60)
print("  Step 10: Saving results")
print("="*60)

# TODO: Save the effect matrix with cluster labels
# Hint:
#   df_output = effect_matrix.copy()
#   df_output["cluster_kmeans"] = best_labels
#   df_output["cluster_gmm"] = gmm_results[BEST_K]["labels"]
#   df_output.to_csv(RESULTS_DIR / "effect_matrix_clustered.csv")
#   print(f"  Saved: {RESULTS_DIR / 'effect_matrix_clustered.csv'}")


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
