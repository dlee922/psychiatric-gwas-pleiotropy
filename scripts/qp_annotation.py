"""
Phase 2B: Q_P Annotation of Clusters
=====================================
Annotate the k=3 sign-aligned clusters with QSNP heterogeneity
values from the factor files, and test whether _Q hits from the
hits file separate differently across clusters.

Usage: python scripts/qp_annotation.py
"""
import sys
sys.path.insert(0, ".")

from src.utils import setup_output, teardown_output
output_path = setup_output("qp_annotation_output.txt")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import gc

from pathlib import Path
from scipy import stats

from src.config import CROSS_DISORDER
from src.data_utils import load_sumstats

RESULTS_DIR = Path("data/results")
FIGURES_DIR = Path("figures")


# ============================================
# STEP 1: Load the clustered effect matrix
# ============================================
print("="*60)
print("  Step 1: Loading clustered effect matrix")
print("="*60)

df_clustered = pd.read_csv(RESULTS_DIR / "effect_matrix_clustered.csv", index_col=0)
print(f"  Shape: {df_clustered.shape}")
print(f"  Columns: {df_clustered.columns.tolist()}")

# We need k=3 labels — re-run k-means with k=3 on the sign-aligned data
# First reconstruct the sign-aligned matrix from the stored data
# (the saved file has sign-aligned z-scores + k=2 labels)
factor_cols = ["F1_Compulsive", "F2_SCZ_BIP", "F3_Neurodev", "F4_Internalizing", "F5_Substance"]
effect_matrix_aligned = df_clustered[factor_cols]

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

scaler = StandardScaler()
X_scaled = scaler.fit_transform(effect_matrix_aligned)

km3 = KMeans(n_clusters=3, n_init=50, random_state=42)
labels_k3 = km3.fit_predict(X_scaled)
df_clustered["cluster_k3"] = labels_k3

# Print cluster sizes
print(f"\n  k=3 cluster sizes:")
for i in range(3):
    n = (labels_k3 == i).sum()
    mean_profile = effect_matrix_aligned.loc[labels_k3 == i].mean().round(2).to_dict()
    print(f"    Cluster {i} (n={n}): {mean_profile}")


# ============================================
# STEP 2: Extract Q_P values for hit SNPs
# ============================================
print("\n" + "="*60)
print("  Step 2: Extracting Q_P values from factor files")
print("="*60)

FACTOR_KEYS = {
    "F1_Compulsive":     "cdg2025_F1_compulsive",
    "F2_SCZ_BIP":        "cdg2025_F2_scz_bip",
    "F3_Neurodev":       "cdg2025_F3_neurodev",
    "F4_Internalizing":  "cdg2025_F4_internalizing",
    "F5_Substance":      "cdg2025_F5_substance",
}

hit_snps = df_clustered.index.values
factor_qp = {}

for factor_name, config_key in FACTOR_KEYS.items():
    print(f"\n  Loading Q_P from {factor_name}...")
    df_factor = load_sumstats(CROSS_DISORDER[config_key]["path"])
    filtered = df_factor[df_factor["SNP"].isin(hit_snps)]
    factor_qp[f"QP_{factor_name}"] = filtered.set_index("SNP")["Q_P"]
    print(f"    Found {len(filtered)} / {len(hit_snps)} hit SNPs")
    del df_factor
    gc.collect()

qp_matrix = pd.DataFrame(factor_qp)
print(f"\n  Q_P matrix shape: {qp_matrix.shape}")
print(f"  Missing values:\n{qp_matrix.isnull().sum()}")

# Merge with clustered data
df_annotated = df_clustered.join(qp_matrix)
print(f"\n  Annotated matrix shape: {df_annotated.shape}")


# ============================================
# STEP 3: Q_P distribution by cluster
# ============================================
print("\n" + "="*60)
print("  Step 3: Q_P values by cluster")
print("="*60)

qp_cols = [c for c in df_annotated.columns if c.startswith("QP_")]

print("\n  Mean Q_P per cluster (low Q_P = heterogeneous, high = homogeneous):")
qp_by_cluster = df_annotated.groupby("cluster_k3")[qp_cols].mean()
print(qp_by_cluster.round(4).to_string())

print("\n  Median Q_P per cluster:")
qp_median = df_annotated.groupby("cluster_k3")[qp_cols].median()
print(qp_median.round(4).to_string())

# Count SNPs with significant heterogeneity (Q_P < 0.05) per cluster per factor
print("\n  SNPs with significant heterogeneity (Q_P < 0.05) per cluster:")
for qp_col in qp_cols:
    print(f"\n  {qp_col}:")
    for cluster_id in sorted(df_annotated["cluster_k3"].unique()):
        subset = df_annotated[df_annotated["cluster_k3"] == cluster_id]
        n_het = (subset[qp_col] < 0.05).sum()
        n_total = len(subset)
        print(f"    Cluster {cluster_id}: {n_het}/{n_total} ({n_het/n_total*100:.1f}%)")

# Heatmap of mean Q_P by cluster
fig, ax = plt.subplots(figsize=(10, 5))
sns.heatmap(qp_by_cluster, annot=True, fmt=".3f", cmap="RdYlGn",
            ax=ax, vmin=0, vmax=1)
ax.set_title("Mean Q_P (Heterogeneity) by Cluster and Factor")
ax.set_ylabel("Cluster (k=3)")
fig.savefig(FIGURES_DIR / "qp_by_cluster_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\n  Saved: {FIGURES_DIR / 'qp_by_cluster_heatmap.png'}")

# Box plots of Q_P distributions per cluster for each factor
fig, axes = plt.subplots(1, 5, figsize=(20, 5), sharey=True)
for i, qp_col in enumerate(qp_cols):
    factor_label = qp_col.replace("QP_", "")
    data_to_plot = [df_annotated[df_annotated["cluster_k3"] == c][qp_col].dropna().values
                    for c in sorted(df_annotated["cluster_k3"].unique())]
    bp = axes[i].boxplot(data_to_plot, labels=[f"C{c}" for c in sorted(df_annotated["cluster_k3"].unique())])
    axes[i].set_title(factor_label, fontsize=10)
    axes[i].set_xlabel("Cluster")
    if i == 0:
        axes[i].set_ylabel("Q_P value")
    axes[i].axhline(y=0.05, color="red", linestyle="--", alpha=0.5, linewidth=1)
fig.suptitle("Q_P Distribution by Cluster (red line = 0.05 significance threshold)")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "qp_boxplots_by_cluster.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGURES_DIR / 'qp_boxplots_by_cluster.png'}")


# ============================================
# STEP 4: Statistical tests
# ============================================
print("\n" + "="*60)
print("  Step 4: Statistical tests — Q_P differences across clusters")
print("="*60)

# Kruskal-Wallis test (non-parametric) for each Q_P column across clusters
print("\n  Kruskal-Wallis test (do clusters differ in Q_P?):")
for qp_col in qp_cols:
    groups = [df_annotated[df_annotated["cluster_k3"] == c][qp_col].dropna().values
              for c in sorted(df_annotated["cluster_k3"].unique())]
    stat, pval = stats.kruskal(*groups)
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    print(f"    {qp_col}: H={stat:.2f}, p={pval:.4e} {sig}")


# ============================================
# STEP 5: Cross-reference with hits file factor labels
# ============================================
print("\n" + "="*60)
print("  Step 5: Hits file factor labels vs clusters")
print("="*60)

df_hits = load_sumstats(CROSS_DISORDER["cdg2025_hits"]["path"])

# For each SNP, get the list of factors it was significant for
snp_factors = df_hits.groupby("SNP")["gwas_name"].apply(list).to_dict()

# Separate _Q hits from regular hits
df_annotated["is_Q_hit"] = df_annotated.index.map(
    lambda snp: any("_Q" in f for f in snp_factors.get(snp, []))
)
df_annotated["is_regular_hit"] = df_annotated.index.map(
    lambda snp: any("_Q" not in f for f in snp_factors.get(snp, []))
)
df_annotated["is_both"] = df_annotated["is_Q_hit"] & df_annotated["is_regular_hit"]

print(f"\n  SNPs that are _Q hits: {df_annotated['is_Q_hit'].sum()}")
print(f"  SNPs that are regular hits: {df_annotated['is_regular_hit'].sum()}")
print(f"  SNPs that are both: {df_annotated['is_both'].sum()}")

# How do _Q hits distribute across clusters?
print("\n  _Q hit distribution across clusters:")
q_by_cluster = pd.crosstab(df_annotated["cluster_k3"], df_annotated["is_Q_hit"],
                            margins=True)
q_by_cluster.columns = ["Regular only", "Has _Q hit", "Total"]
print(q_by_cluster.to_string())

# Proportions
print("\n  Proportion of _Q hits per cluster:")
for cluster_id in sorted(df_annotated["cluster_k3"].unique()):
    subset = df_annotated[df_annotated["cluster_k3"] == cluster_id]
    n_q = subset["is_Q_hit"].sum()
    n_total = len(subset)
    print(f"    Cluster {cluster_id}: {n_q}/{n_total} ({n_q/n_total*100:.1f}%)")

# Chi-square test: are _Q hits unevenly distributed across clusters?
contingency = pd.crosstab(df_annotated["cluster_k3"], df_annotated["is_Q_hit"])
chi2, pval, dof, expected = stats.chi2_contingency(contingency)
print(f"\n  Chi-square test (_Q hits vs cluster): chi2={chi2:.2f}, p={pval:.4e}, dof={dof}")

# Which specific factor _Q hits land in which clusters?
print("\n  Detailed: which factor's _Q hits land in which cluster:")
for factor_q in ["F1_Q", "F2_Q", "F3_Q", "F4_Q", "F5_Q", "Hier_Q"]:
    snps_in_factor = df_hits[df_hits["gwas_name"] == factor_q]["SNP"].unique()
    snps_in_matrix = [s for s in snps_in_factor if s in df_annotated.index]
    if len(snps_in_matrix) > 0:
        cluster_dist = df_annotated.loc[snps_in_matrix, "cluster_k3"].value_counts().sort_index()
        print(f"    {factor_q} ({len(snps_in_matrix)} SNPs): {cluster_dist.to_dict()}")


# ============================================
# STEP 6: Factor-specific hit labels by cluster
# ============================================
print("\n" + "="*60)
print("  Step 6: Which factors' hits populate each cluster?")
print("="*60)

# For each regular (non-Q) factor, where do its hits land?
for factor_label in ["F1", "F2", "F3", "F4", "F5", "Hier"]:
    snps_in_factor = df_hits[df_hits["gwas_name"] == factor_label]["SNP"].unique()
    snps_in_matrix = [s for s in snps_in_factor if s in df_annotated.index]
    if len(snps_in_matrix) > 0:
        cluster_dist = df_annotated.loc[snps_in_matrix, "cluster_k3"].value_counts().sort_index()
        print(f"  {factor_label} hits ({len(snps_in_matrix)} SNPs): {cluster_dist.to_dict()}")


# ============================================
# STEP 7: Save annotated results
# ============================================
print("\n" + "="*60)
print("  Step 7: Saving annotated results")
print("="*60)

df_annotated.to_csv(RESULTS_DIR / "effect_matrix_qp_annotated.csv")
print(f"  Saved: {RESULTS_DIR / 'effect_matrix_qp_annotated.csv'}")

# Summary
print("\n" + "="*60)
print("  Summary")
print("="*60)
print("""
  Key questions answered:
  1. Do clusters differ in Q_P heterogeneity? (Kruskal-Wallis tests)
  2. Do _Q hits concentrate in specific clusters? (Chi-square test)
  3. Which factors' hits populate each cluster? (Factor label distribution)

  If the antagonistic pleiotropy cluster (Cluster 2, 74 loci) is enriched
  for _Q hits, that validates our clustering — it means the loci our ML
  identified as "breaking the rules" are the same ones the QSNP test
  flagged as heterogeneous.
""")


if __name__ == "__main__":
    teardown_output(output_path)
    print("Q_P Annotation Complete")