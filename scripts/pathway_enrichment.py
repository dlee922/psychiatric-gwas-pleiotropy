"""
Pathway Enrichment Analysis
======================================
For each cluster, map lead SNPs to nearest genes,
then run GO and pathway enrichment using g:Profiler.

Usage: python scripts/pathway_enrichment.py
"""
import sys
sys.path.insert(0, ".")

from src.utils import setup_output, teardown_output
output_path = setup_output("pathway_enrichment_output.txt")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from gprofiler import GProfiler

from src.config import CROSS_DISORDER
from src.data_utils import load_sumstats

RESULTS_DIR = Path("data/results")
FIGURES_DIR = Path("figures")


# ============================================
# 1. Load annotated data and factor files
# ============================================
print("="*60)
print("  Step 1: Loading data")
print("="*60)

df_annotated = pd.read_csv(RESULTS_DIR / "effect_matrix_qp_annotated.csv", index_col=0)
print(f"  Annotated matrix shape: {df_annotated.shape}")

# Re-derive k=3 labels (same as qp_annotation.py)
factor_cols = ["F1_Compulsive", "F2_SCZ_BIP", "F3_Neurodev", "F4_Internalizing", "F5_Substance"]
effect_matrix = df_annotated[factor_cols]

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

scaler = StandardScaler()
X_scaled = scaler.fit_transform(effect_matrix)
km3 = KMeans(n_clusters=3, n_init=50, random_state=42)
labels_k3 = km3.fit_predict(X_scaled)
df_annotated["cluster_k3"] = labels_k3

hit_snps = df_annotated.index.values

# Print cluster profiles as reminder
print("\n  k=3 cluster profiles:")
for i in range(3):
    n = (labels_k3 == i).sum()
    profile = effect_matrix.loc[labels_k3 == i].mean().round(2).to_dict()
    print(f"    Cluster {i} (n={n}): {profile}")


# ============================================
# 2. Map SNPs to genes using factor file positions
# ============================================
print("\n" + "="*60)
print("  Step 2: Getting SNP positions")
print("="*60)

# Load one factor file to get CHR and BP for each SNP
df_factor = load_sumstats(CROSS_DISORDER["cdg2025_F1_compulsive"]["path"])
snp_positions = df_factor[df_factor["SNP"].isin(hit_snps)].set_index("SNP")[["CHR", "BP"]]
del df_factor

print(f"  SNP positions retrieved: {len(snp_positions)}")

# Use the SNP rsIDs directly with g:Profiler 
# Accepts rsIDs and maps them to nearby genes internally.

# Try querying g:Profiler with rsIDs directly
# g:Profiler's convert function can map SNPs to genes
print("\n  Testing g:Profiler SNP to gene mapping...")

gp = GProfiler(return_dataframe=True)

# Test with a small batch
test_snps = hit_snps[:10].tolist()
try:
    test_convert = gp.convert(
        organism="hsapiens",
        query=test_snps,
        target_namespace="ENSG"
    )
    print(f"  Test conversion result: {len(test_convert)} mappings for {len(test_snps)} SNPs")
    print(f"  Columns: {test_convert.columns.tolist()}")
    if len(test_convert) > 0:
        print(f"  Sample:\n{test_convert.head(3).to_string()}")
    snp_to_gene_method = "convert"
except Exception as e:
    print(f"  g:Profiler convert failed: {e}")
    print(f"  Will query with rsIDs directly in enrichment")
    snp_to_gene_method = "direct"
# Confimed working

# ============================================
# 3. Map all SNPs to genes
# ============================================
print("\n" + "="*60)
print("  Step 3: Mapping all SNPs to genes")
print("="*60)

if snp_to_gene_method == "convert":
    # Convert all hit SNPs to genes
    all_conversions = gp.convert(
        organism="hsapiens",
        query=hit_snps.tolist(),
        target_namespace="ENSG"
    )
    print(f"  Total conversions: {len(all_conversions)}")

    # Filter to only successful conversions
    # Handle both None and "None" cases
    valid = all_conversions[
        all_conversions["converted"].notna() &
        (all_conversions["converted"] != "None") &
        (all_conversions["converted"].str.startswith("ENSG", na=False))
    ]
    print(f"  Valid gene mappings: {len(valid)}")

    # Build SNP -> gene mapping
    # The 'incoming' column may have the original rsID or a positional format
    # Need to map back to original rsIDs
    # g:Profiler returns results in the same order as the query
    snp_to_genes = {}
    for _, row in valid.iterrows():
        # The query column groups results by query batch
        # incoming is the identifier g:Profiler used
        incoming = row["incoming"]
        gene = row["converted"]

        # Check if incoming is one of the rsIDs directly
        if incoming in hit_snps:
            if incoming not in snp_to_genes:
                snp_to_genes[incoming] = []
            snp_to_genes[incoming].append(gene)
        else:
            # g:Profiler may convert rsID to position format
            # Use the n_incoming index to map back
            idx = int(row["n_incoming"]) - 1  # 1-indexed to 0-indexed
            if idx < len(hit_snps):
                snp = hit_snps[idx]
                if snp not in snp_to_genes:
                    snp_to_genes[snp] = []
                snp_to_genes[snp].append(gene)

    # Map each SNP to its gene(s)
    df_annotated["genes"] = df_annotated.index.map(
        lambda snp: snp_to_genes.get(snp, [])
    )

    mapped = df_annotated["genes"].apply(len).gt(0).sum()
    total_genes = len(set(g for genes in df_annotated["genes"] for g in genes))
    print(f"  SNPs mapped to genes: {mapped} / {len(df_annotated)}")
    print(f"  Total unique genes: {total_genes}")

    if mapped == 0:
        # Fallback: skip conversion, use rsIDs directly in enrichment
        print("\n  WARNING: Gene mapping failed, falling back to direct rsID queries")
        snp_to_gene_method = "direct"

if snp_to_gene_method == "direct" or mapped == 0:
    # g:Profiler's profile function can accept rsIDs directly
    # It will map them to nearby genes internally
    print("  Using rsIDs directly for enrichment (g:Profiler handles mapping)")

# Build gene lists per cluster regardless of method
cluster_genes = {}
for cluster_id in range(3):
    if snp_to_gene_method == "convert" and mapped > 0:
        genes = []
        for gene_list in df_annotated[df_annotated["cluster_k3"] == cluster_id]["genes"]:
            genes.extend(gene_list)
        cluster_genes[cluster_id] = list(set(genes))
    else:
        # Use rsIDs directly
        cluster_genes[cluster_id] = df_annotated[df_annotated["cluster_k3"] == cluster_id].index.tolist()

    print(f"  Cluster {cluster_id}: {len(cluster_genes[cluster_id])} {'genes' if snp_to_gene_method == 'convert' and mapped > 0 else 'SNPs'}")


# ============================================
# 4. Run enrichment analysis per cluster
# ============================================
print("\n" + "="*60)
print("  Step 4: Running pathway enrichment per cluster")
print("="*60)

# Background: all genes from all clusters combined
all_genes = []
for genes in cluster_genes.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))
print(f"  Background gene set: {len(all_genes)} genes")

cluster_enrichment = {}

for cluster_id in range(3):
    query_genes = cluster_genes[cluster_id]
    n_genes = len(query_genes)
    print(f"\n  Cluster {cluster_id} ({n_genes} genes/SNPs):")

    if n_genes < 5:
        print(f"    Too few genes for enrichment analysis, skipping")
        continue

    try:
        # If querying with rsIDs, don't use a custom background
        # (g:Profiler needs gene-level background, not SNP-level)
        if snp_to_gene_method == "direct" or mapped == 0:
            results = gp.profile(
                organism="hsapiens",
                query=query_genes,
                sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
                significance_threshold_method="fdr",
                user_threshold=0.05,
                no_evidences=False,
            )
        else:
            results = gp.profile(
                organism="hsapiens",
                query=query_genes,
                background=all_genes,
                sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
                significance_threshold_method="fdr",
                user_threshold=0.05,
                no_evidences=False,
            )

        if len(results) > 0:
            cluster_enrichment[cluster_id] = results
            print(f"    Significant terms: {len(results)}")
            print(f"    Sources: {results['source'].value_counts().to_dict()}")

            # Top 10 terms
            top = results.nsmallest(10, "p_value")[["source", "native", "name", "p_value", "intersection_size"]]
            print(f"\n    Top 10 enriched terms:")
            for _, row in top.iterrows():
                print(f"      [{row['source']}] {row['name']} (p={row['p_value']:.2e}, genes={row['intersection_size']})")
        else:
            print(f"    No significant enrichment found")
            cluster_enrichment[cluster_id] = pd.DataFrame()

    except Exception as e:
        print(f"    Error in enrichment: {e}")
        cluster_enrichment[cluster_id] = pd.DataFrame()


# ============================================
# 5. Compare enrichment across clusters
# ============================================
print("\n" + "="*60)
print("  Step 5: Cross-cluster enrichment comparison")
print("="*60)

# Collect all significant terms and check overlap
all_terms = {}
for cluster_id, results in cluster_enrichment.items():
    if len(results) > 0:
        terms = set(results["native"].tolist())
        all_terms[cluster_id] = terms
        print(f"  Cluster {cluster_id}: {len(terms)} significant terms")

if len(all_terms) >= 2:
    cluster_ids = sorted(all_terms.keys())
    for i in range(len(cluster_ids)):
        for j in range(i+1, len(cluster_ids)):
            ci, cj = cluster_ids[i], cluster_ids[j]
            overlap = all_terms[ci] & all_terms[cj]
            unique_i = all_terms[ci] - all_terms[cj]
            unique_j = all_terms[cj] - all_terms[ci]
            print(f"\n  Cluster {ci} vs Cluster {cj}:")
            print(f"    Shared terms: {len(overlap)}")
            print(f"    Unique to Cluster {ci}: {len(unique_i)}")
            print(f"    Unique to Cluster {cj}: {len(unique_j)}")

    # Terms unique to each cluster
    for cluster_id in cluster_ids:
        other_terms = set()
        for other_id in cluster_ids:
            if other_id != cluster_id:
                other_terms |= all_terms.get(other_id, set())
        unique = all_terms[cluster_id] - other_terms

        if len(unique) > 0:
            unique_results = cluster_enrichment[cluster_id][
                cluster_enrichment[cluster_id]["native"].isin(unique)
            ].nsmallest(10, "p_value")
            print(f"\n  Top terms UNIQUE to Cluster {cluster_id}:")
            for _, row in unique_results.iterrows():
                print(f"    [{row['source']}] {row['name']} (p={row['p_value']:.2e})")


# ============================================
# 6. Visualization
# ============================================
print("\n" + "="*60)
print("  Step 6: Enrichment visualization")
print("="*60)

# Bar plot of top terms per cluster
fig, axes = plt.subplots(3, 1, figsize=(12, 15))

for cluster_id in range(3):
    ax = axes[cluster_id]
    if cluster_id in cluster_enrichment and len(cluster_enrichment[cluster_id]) > 0:
        top = cluster_enrichment[cluster_id].nsmallest(15, "p_value")
        y_pos = range(len(top))
        bars = ax.barh(y_pos, -np.log10(top["p_value"]), color=f"C{cluster_id}")
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top["name"], fontsize=8)
        ax.set_xlabel("-log10(p-value)")
        ax.invert_yaxis()

        n_loci = (labels_k3 == cluster_id).sum()
        profile = effect_matrix.loc[labels_k3 == cluster_id].mean().round(1).to_dict()
        ax.set_title(f"Cluster {cluster_id} (n={n_loci}): Top Enriched Pathways")
    else:
        ax.text(0.5, 0.5, f"Cluster {cluster_id}: No significant enrichment",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_title(f"Cluster {cluster_id}")

fig.tight_layout()
fig.savefig(FIGURES_DIR / "pathway_enrichment_by_cluster.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGURES_DIR / 'pathway_enrichment_by_cluster.png'}")


# ============================================
# 7. Save results
# ============================================
print("\n" + "="*60)
print("  Step 7: Saving enrichment results")
print("="*60)

for cluster_id, results in cluster_enrichment.items():
    if len(results) > 0:
        filepath = RESULTS_DIR / f"enrichment_cluster{cluster_id}.csv"
        results.to_csv(filepath, index=False)
        print(f"  Saved: {filepath}")

# Summary table
summary_rows = []
for cluster_id in range(3):
    n_loci = (labels_k3 == cluster_id).sum()
    n_terms = len(cluster_enrichment.get(cluster_id, pd.DataFrame()))
    profile = effect_matrix.loc[labels_k3 == cluster_id].mean().round(2).to_dict()

    if cluster_id in cluster_enrichment and len(cluster_enrichment[cluster_id]) > 0:
        top_term = cluster_enrichment[cluster_id].nsmallest(1, "p_value").iloc[0]
        top_pathway = f"{top_term['name']} (p={top_term['p_value']:.2e})"
    else:
        top_pathway = "None"

    summary_rows.append({
        "Cluster": cluster_id,
        "N_Loci": n_loci,
        "N_Enriched_Terms": n_terms,
        "Top_Pathway": top_pathway,
        **{f"Mean_Z_{k}": v for k, v in profile.items()}
    })

df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(RESULTS_DIR / "cluster_summary_with_enrichment.csv", index=False)
print(f"  Saved: {RESULTS_DIR / 'cluster_summary_with_enrichment.csv'}")

print("\n" + "="*60)
print("  Enrichment Analysis Complete")
print("="*60)


if __name__ == "__main__":
    teardown_output(output_path)
    print("Pathway Enrichment Complete")