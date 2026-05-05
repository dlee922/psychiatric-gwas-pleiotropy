"""
Phase 2: Cross-Factor Effect Matrix & Clustering
=================================================
Build the locus x factor z-score matrix from CDG2025 data,
then run k-means (primary), hierarchical, GMM, and DBSCAN clustering.

Approach: Sign-align all SNPs to positive F4 (Internalizing),
preserving directional patterns across factors.

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

# Primary cluster count — change this single value to update all downstream outputs
BEST_K = 3


# ============================================
# STEP 1: Load hits file and extract unique SNPs
# ============================================
print("="*60)
print("  Step 1: Loading hits file")
print("="*60)

df_hits = load_sumstats(CROSS_DISORDER["cdg2025_hits"]["path"])
print(f"  Total loci: {len(df_hits)}")
print(f"  Columns: {df_hits.columns.tolist()}")

hit_snps = df_hits["SNP"].unique()
print(f"  Unique lead SNPs: {len(hit_snps)}")
print(f"  Total rows ({len(df_hits)}) vs unique SNPs ({len(hit_snps)}) = {len(df_hits) - len(hit_snps)} duplicated across factors")


# ============================================
# STEP 2: Build cross-factor effect matrix
# ============================================
print("\n" + "="*60)
print("  Step 2: Building cross-factor z-score matrix")
print("="*60)

FACTOR_KEYS = {
    "F1_Compulsive":     "cdg2025_F1_compulsive",
    "F2_SCZ_BIP":        "cdg2025_F2_scz_bip",
    "F3_Neurodev":       "cdg2025_F3_neurodev",
    "F4_Internalizing":  "cdg2025_F4_internalizing",
    "F5_Substance":      "cdg2025_F5_substance",
}

factor_zscores = {}

for factor_name, config_key in FACTOR_KEYS.items():
    print(f"\n  Loading {factor_name}...")
    df_factor = load_sumstats(CROSS_DISORDER[config_key]["path"])
    filtered = df_factor[df_factor["SNP"].isin(hit_snps)]
    factor_zscores[factor_name] = filtered.set_index("SNP").eval("BETA / SE")
    print(f"    Found {len(filtered)} / {len(hit_snps)} hit SNPs")
    del df_factor
    gc.collect()

effect_matrix = pd.DataFrame(factor_zscores)

print(f"\n  Raw effect matrix shape: {effect_matrix.shape}")
print(f"  Expected: ({len(hit_snps)}, 5)")
print(f"  Missing values per column:\n{effect_matrix.isnull().sum()}")

# Sign-alignment: flip each SNP so F4_Internalizing is always positive.
# This removes the arbitrary allele direction while preserving the
# relative pattern of effects across factors.
sign = np.sign(effect_matrix["F4_Internalizing"])
# Handle any exact zeros (unlikely but safe)
sign = sign.replace(0, 1)
effect_matrix_aligned = effect_matrix.multiply(sign, axis=0)

n_flipped = (sign == -1).sum()
print(f"\n  Sign-alignment to F4_Internalizing:")
print(f"    SNPs flipped: {n_flipped} / {len(effect_matrix)} ({n_flipped/len(effect_matrix)*100:.1f}%)")
print(f"    F4 now all positive: {(effect_matrix_aligned['F4_Internalizing'] >= 0).all()}")

# Verification
sample_snp = effect_matrix_aligned.index[0]
print(f"\n  Verification — {sample_snp}:")
print(f"    Raw z-scores:     {effect_matrix.loc[sample_snp].round(4).to_dict()}")
print(f"    Aligned z-scores: {effect_matrix_aligned.loc[sample_snp].round(4).to_dict()}")


# ============================================
# STEP 3: Quick sanity checks
# ============================================
print("\n" + "="*60)
print("  Step 3: Effect matrix sanity checks (sign-aligned)")
print("="*60)

print("\n  Descriptive statistics:")
print(effect_matrix_aligned.describe().round(4).to_string())

print("\n  Inter-factor correlations:")
print(effect_matrix_aligned.corr().round(4).to_string())

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sns.heatmap(effect_matrix.corr(), annot=True, fmt=".3f",
            cmap="RdBu_r", center=0, ax=axes[0], vmin=-1, vmax=1)
axes[0].set_title("Raw Z-scores (before alignment)")

sns.heatmap(effect_matrix_aligned.corr(), annot=True, fmt=".3f",
            cmap="RdBu_r", center=0, ax=axes[1], vmin=-1, vmax=1)
axes[1].set_title("Sign-aligned Z-scores")

fig.tight_layout()
fig.savefig(FIGURES_DIR / "factor_correlation_comparison.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\n  Saved: {FIGURES_DIR / 'factor_correlation_comparison.png'}")


# ============================================
# STEP 4: Standardize
# ============================================
print("\n" + "="*60)
print("  Step 4: Standardization")
print("="*60)

from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
X_scaled = scaler.fit_transform(effect_matrix_aligned)

print(f"  Scaled matrix shape: {X_scaled.shape}")
print(f"  Column means (should be ~0): {X_scaled.mean(axis=0).round(6)}")
print(f"  Column stds (should be ~1):  {X_scaled.std(axis=0).round(6)}")


# ============================================
# STEP 5: K-Means clustering (PRIMARY METHOD)
# ============================================
print("\n" + "="*60)
print("  Step 5: K-Means clustering (primary method)")
print("="*60)

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

k_range = range(2, 11)
kmeans_results = {}

for k in k_range:
    km = KMeans(n_clusters=k, n_init=50, random_state=42)
    labels = km.fit_predict(X_scaled)
    sil = silhouette_score(X_scaled, labels)
    kmeans_results[k] = {"model": km, "labels": labels, "silhouette": sil}
    print(f"    k={k}: silhouette={sil:.4f}")

fig, ax = plt.subplots(figsize=(8, 5))
sil_scores = [kmeans_results[k]["silhouette"] for k in k_range]
ax.plot(list(k_range), sil_scores, "bo-", linewidth=2, markersize=8)
ax.set_xlabel("Number of clusters (k)")
ax.set_ylabel("Silhouette Score")
ax.set_title("K-Means: Silhouette Score vs k (Sign-aligned)")
ax.set_xticks(list(k_range))
best_sil_k = list(k_range)[np.argmax(sil_scores)]
ax.axvline(x=best_sil_k, color="red", linestyle="--", alpha=0.5, label=f"Best k={best_sil_k}")
ax.legend()
fig.savefig(FIGURES_DIR / "kmeans_silhouette.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\n  Best silhouette k={best_sil_k} (score={max(sil_scores):.4f})")
print(f"  Saved: {FIGURES_DIR / 'kmeans_silhouette.png'}")


# ============================================
# STEP 6: Hierarchical clustering (comparison method)
# ============================================
print("\n" + "="*60)
print("  Step 6: Hierarchical clustering (comparison method)")
print("="*60)

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

distances = pdist(X_scaled, metric="euclidean")
Z_linkage = linkage(distances, method="ward")

fig, ax = plt.subplots(figsize=(14, 6))
dendrogram(Z_linkage, ax=ax, truncate_mode="lastp", p=30,
           leaf_rotation=90, leaf_font_size=8)
ax.set_title("Hierarchical Clustering Dendrogram (Ward, Euclidean) — Sign-aligned")
ax.set_xlabel("Cluster (size)")
ax.set_ylabel("Distance")
fig.savefig(FIGURES_DIR / "dendrogram_ward.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGURES_DIR / 'dendrogram_ward.png'}")

for n_clust in [2, 3, 4, 5, 6]:
    labels_hier = fcluster(Z_linkage, n_clust, criterion="maxclust")
    sizes = pd.Series(labels_hier).value_counts().sort_index()
    print(f"    k={n_clust}: cluster sizes = {sizes.tolist()}")


# ============================================
# STEP 7: Gaussian Mixture Models (comparison method)
# ============================================
print("\n" + "="*60)
print("  Step 7: Gaussian Mixture Models (comparison method)")
print("="*60)

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

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
bics = [gmm_results[k]["bic"] for k in k_range]
ax1.plot(list(k_range), bics, "go-", linewidth=2, markersize=8)
ax1.set_xlabel("Number of components (k)")
ax1.set_ylabel("BIC (lower = better)")
ax1.set_title("GMM: BIC vs k (Sign-aligned)")
ax1.set_xticks(list(k_range))

gmm_sils = [gmm_results[k]["silhouette"] for k in k_range]
ax2.plot(list(k_range), gmm_sils, "ro-", linewidth=2, markersize=8)
ax2.set_xlabel("Number of components (k)")
ax2.set_ylabel("Silhouette Score")
ax2.set_title("GMM: Silhouette vs k (Sign-aligned)")
ax2.set_xticks(list(k_range))

fig.tight_layout()
fig.savefig(FIGURES_DIR / "gmm_bic_silhouette.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\n  Best BIC k={list(k_range)[np.argmin(bics)]}")
print(f"  Best silhouette k={list(k_range)[np.argmax(gmm_sils)]}")
print(f"  Saved: {FIGURES_DIR / 'gmm_bic_silhouette.png'}")


# ============================================
# STEP 8: DBSCAN (comparison method)
# ============================================
print("\n" + "="*60)
print("  Step 8: DBSCAN (comparison method)")
print("="*60)

from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

nn = NearestNeighbors(n_neighbors=5)
nn.fit(X_scaled)
distances_nn, _ = nn.kneighbors(X_scaled)
k_distances = np.sort(distances_nn[:, -1])

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(k_distances)
ax.set_xlabel("Points (sorted)")
ax.set_ylabel("5th Nearest Neighbor Distance")
ax.set_title("K-Distance Graph (Sign-aligned)")
fig.savefig(FIGURES_DIR / "dbscan_kdistance.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGURES_DIR / 'dbscan_kdistance.png'}")

for eps_val in [0.5, 1.0, 1.5, 2.0, 2.5]:
    db = DBSCAN(eps=eps_val, min_samples=5)
    labels = db.fit_predict(X_scaled)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = (labels == -1).sum()
    print(f"    eps={eps_val}: {n_clusters} clusters, {n_noise} noise points")


# ============================================
# STEP 9: Visualization with primary k
# ============================================
print("\n" + "="*60)
print(f"  Step 9: Visualization (k-means, k={BEST_K})")
print("="*60)

best_labels = kmeans_results[BEST_K]["labels"]
print(f"  Primary method: K-Means, k={BEST_K}")
print(f"  Silhouette score: {kmeans_results[BEST_K]['silhouette']:.4f}")
print(f"  Cluster sizes: {pd.Series(best_labels).value_counts().sort_index().tolist()}")

# --- 9A: Heatmap sorted by cluster ---
df_plot = effect_matrix_aligned.copy()
df_plot["cluster"] = best_labels
df_plot = df_plot.sort_values("cluster")

sorted_index = df_plot.index
X_sorted = X_scaled[effect_matrix_aligned.index.get_indexer(sorted_index)]

fig, ax = plt.subplots(figsize=(10, 14))
sns.heatmap(
    pd.DataFrame(X_sorted, columns=FACTOR_KEYS.keys()),
    cmap="RdBu_r", center=0, ax=ax,
    yticklabels=False,
    cbar_kws={"label": "Standardized Z-score (sign-aligned)"}
)
cluster_sizes = df_plot["cluster"].value_counts().sort_index()
cumulative = 0
for size in cluster_sizes.values[:-1]:
    cumulative += size
    ax.axhline(y=cumulative, color="black", linewidth=2)
ax.set_title(f"Cross-Factor Effect Profiles — Sign-aligned (K-Means, k={BEST_K})")
ax.set_ylabel("Loci")
fig.savefig(FIGURES_DIR / "heatmap_clustered.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGURES_DIR / 'heatmap_clustered.png'}")

# --- 9B: UMAP projection ---
import umap

reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding = reducer.fit_transform(X_scaled)

fig, ax = plt.subplots(figsize=(10, 8))
scatter = ax.scatter(embedding[:, 0], embedding[:, 1],
                     c=best_labels, cmap="tab10", s=15, alpha=0.7)
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title(f"UMAP Projection — Sign-aligned (K-Means, k={BEST_K})")
plt.colorbar(scatter, label="Cluster")
fig.savefig(FIGURES_DIR / "umap_clusters.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGURES_DIR / 'umap_clusters.png'}")

# --- 9C: Cluster profiles ---
df_profiles = effect_matrix_aligned.copy()
df_profiles["cluster"] = best_labels
cluster_means = df_profiles.groupby("cluster").mean()
print(f"\n  Cluster mean z-scores (sign-aligned, k-means k={BEST_K}):")
print(cluster_means.round(4).to_string())

cluster_stds = df_profiles.groupby("cluster").std()
print(f"\n  Cluster std z-scores:")
print(cluster_stds.round(4).to_string())

fig, axes = plt.subplots(1, BEST_K, figsize=(4*BEST_K, 5), sharey=True)
if BEST_K == 1:
    axes = [axes]
for i, ax in enumerate(axes):
    means = cluster_means.loc[i]
    colors = ["red" if v < 0 else "steelblue" for v in means]
    ax.bar(means.index, means.values, color=colors)
    ax.set_title(f"Cluster {i} (n={sum(best_labels == i)})")
    ax.axhline(y=0, color="black", linewidth=0.5)
    ax.tick_params(axis="x", rotation=45)
fig.suptitle(f"Mean Cross-Factor Z-score Profiles by Cluster (K-Means, k={BEST_K})")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "cluster_profiles.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGURES_DIR / 'cluster_profiles.png'}")


# ============================================
# STEP 10: Cross-method comparison at k=BEST_K
# ============================================
print("\n" + "="*60)
print(f"  Step 10: Cross-method comparison")
print("="*60)

# K-Means across k values
print("\n  K-Means cluster profiles across k=2, 3, 4:")
for k_val in [2, 3, 4]:
    labels_k = kmeans_results[k_val]["labels"]
    df_temp = effect_matrix_aligned.copy()
    df_temp["cluster"] = labels_k
    means = df_temp.groupby("cluster").mean()
    sizes = df_temp.groupby("cluster").size()
    print(f"\n    k={k_val} (silhouette={kmeans_results[k_val]['silhouette']:.4f}):")
    for i in range(k_val):
        profile = means.loc[i].round(2).to_dict()
        print(f"      Cluster {i} (n={sizes[i]}): {profile}")

# Hierarchical at k=BEST_K for comparison
print(f"\n  Hierarchical (Ward) cluster profiles at k={BEST_K}:")
labels_hier = fcluster(Z_linkage, BEST_K, criterion="maxclust")
# Convert to 0-indexed for consistency
labels_hier = labels_hier - 1
df_temp = effect_matrix_aligned.copy()
df_temp["cluster"] = labels_hier
means_hier = df_temp.groupby("cluster").mean()
sizes_hier = df_temp.groupby("cluster").size()
for i in range(BEST_K):
    profile = means_hier.loc[i].round(2).to_dict()
    print(f"    Cluster {i} (n={sizes_hier[i]}): {profile}")

# GMM at k=BEST_K for comparison
print(f"\n  GMM cluster profiles at k={BEST_K}:")
labels_gmm = gmm_results[BEST_K]["labels"]
df_temp = effect_matrix_aligned.copy()
df_temp["cluster"] = labels_gmm
means_gmm = df_temp.groupby("cluster").mean()
sizes_gmm = df_temp.groupby("cluster").size()
for i in range(BEST_K):
    profile = means_gmm.loc[i].round(2).to_dict()
    print(f"    Cluster {i} (n={sizes_gmm[i]}): {profile}")


# ============================================
# STEP 11: Save results
# ============================================
print("\n" + "="*60)
print("  Step 11: Saving results")
print("="*60)

df_output = effect_matrix_aligned.copy()
df_output[f"cluster_kmeans_k{BEST_K}"] = kmeans_results[BEST_K]["labels"]
df_output[f"cluster_hier_k{BEST_K}"] = labels_hier
df_output[f"cluster_gmm_k{BEST_K}"] = gmm_results[BEST_K]["labels"]
df_output.to_csv(RESULTS_DIR / "effect_matrix_clustered.csv")
print(f"  Saved: {RESULTS_DIR / 'effect_matrix_clustered.csv'}")

with open(RESULTS_DIR / "clustering_summary.txt", "w") as f:
    f.write(f"Phase 2: Clustering Analysis Summary (Sign-aligned)\n")
    f.write("="*60 + "\n\n")
    f.write(f"Sign-alignment reference: F4_Internalizing (flipped to positive)\n")
    f.write(f"Hit loci: {len(hit_snps)} unique SNPs\n")
    f.write(f"Effect matrix: {effect_matrix_aligned.shape}\n")
    f.write(f"Primary method: K-Means\n")
    f.write(f"Primary k: {BEST_K}\n\n")
    f.write("K-Means Silhouette Scores:\n")
    for k in k_range:
        marker = " <-- primary" if k == BEST_K else ""
        f.write(f"  k={k}: {kmeans_results[k]['silhouette']:.4f}{marker}\n")
    f.write(f"\nGMM BIC Values:\n")
    for k in k_range:
        f.write(f"  k={k}: {gmm_results[k]['bic']:.1f}\n")
    f.write(f"\nPrimary cluster sizes (k-means, k={BEST_K}):\n")
    for i in range(BEST_K):
        n = sum(kmeans_results[BEST_K]["labels"] == i)
        f.write(f"  Cluster {i}: {n} loci\n")
    f.write(f"\nCluster mean z-scores (sign-aligned, k-means k={BEST_K}):\n")
    f.write(cluster_means.round(4).to_string())
    f.write(f"\n\nCluster std z-scores:\n")
    f.write(cluster_stds.round(4).to_string())
print(f"  Saved: {RESULTS_DIR / 'clustering_summary.txt'}")

print(f"\n  Phase 2 Complete (primary: k-means, k={BEST_K})")


if __name__ == "__main__":
    teardown_output(output_path)
    print("Phase 2: Clustering Analysis Complete")