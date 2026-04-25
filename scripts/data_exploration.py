"""
Phase 1: Data Exploration
=========================
Run a standardized exploration across all datasets.

Usage: python scripts/data_exploration.py
"""
import sys
sys.path.insert(0, ".")

from pathlib import Path
output_file = Path("data/results/exploration_output.txt")
if output_file.exists():
    output_file.unlink()
sys.stdout = open(output_file, "w")

import matplotlib
matplotlib.use("Agg")

from src.config import CROSS_DISORDER, REPLICATION, LATEST
from src.data_utils import (
    explore_dataset, load_sumstats, identify_columns,
    load_cdg2019, gwas_summary, quick_look
)

# ============================================
# 1. CDG2019 Cross-Disorder (8-disorder)
# ============================================
print("="*60)
print("  1. CDG2019: 8-Disorder Cross-Disorder Meta-Analysis")
print("="*60)

df_cdg = load_cdg2019(CROSS_DISORDER["cdg2019"]["path"])
col_map_cdg = identify_columns(df_cdg)
gwas_summary(df_cdg, col_map_cdg, name="cdg2019")

# Note: m-value columns are all NA in this file.
# Significant loci and m-values are in the supplementary data.
# This file is still useful for LDSC (uses all SNPs).
m_cols = [c for c in df_cdg.columns if c.startswith("m.")]
m_populated = sum(df_cdg[col].notna().sum() for col in m_cols)
print(f"  M-value columns present but all NA (total non-null: {m_populated})")
print(f"  Note: significant loci with m-values are in CDG3_Hits or supplementary files")
print(f"  Note: P-values in this file are truncated at ~1e-6")

# Free memory
del df_cdg

# ============================================
# 2. CDG2025 Hits File (significant loci)
# ============================================
print("\n" + "="*60)
print("  2. CDG2025: Significant Loci (Hits File)")
print("="*60)

df_hits = load_sumstats(CROSS_DISORDER["cdg2025_hits"]["path"], verbose=True)
quick_look(df_hits, name="cdg2025_hits")

# ============================================
# 3. CDG2025 Factor Files (14-disorder)
# ============================================
factor_keys = [
    "cdg2025_F1_compulsive",
    "cdg2025_F2_scz_bip",
    "cdg2025_F3_neurodev",
    "cdg2025_F4_internalizing",
    "cdg2025_F5_substance",
    "cdg2025_pfactor",
]

for key in factor_keys:
    print("\n" + "="*60)
    print(f"  3. CDG2025 Factor: {CROSS_DISORDER[key]['description']}")
    print("="*60)

    df_factor = load_sumstats(CROSS_DISORDER[key]["path"], nrows=1000, verbose=True)
    quick_look(df_factor, name=key)
    col_map_factor = identify_columns(df_factor)
    del df_factor

# ============================================
# 4. Replication Datasets (8 disorders, matching 2019 paper)
# ============================================
print("\n" + "="*60)
print("  4. Replication Datasets: Column Comparison")
print("="*60)

for name, info in REPLICATION.items():
    print(f"\n{'='*40}")
    print(f"  {name}: {info['description']}")
    try:
        df_tmp = load_sumstats(info["path"], nrows=5, verbose=True)
        col_map_tmp = identify_columns(df_tmp)
        del df_tmp
    except Exception as e:
        print(f"  ERROR loading: {e}")

# ============================================
# 5. Full exploration of each replication dataset
# ============================================
print("\n" + "="*60)
print("  5. Full Replication Dataset Summaries")
print("="*60)

for name, info in REPLICATION.items():
    print(f"\n{'='*40}")
    print(f"  {name}: {info['description']}")
    try:
        df_tmp = load_sumstats(info["path"], verbose=True)
        col_map_tmp = identify_columns(df_tmp, verbose=False)
        gwas_summary(df_tmp, col_map_tmp, name=name)
        del df_tmp
    except Exception as e:
        print(f"  ERROR loading: {e}")

# ============================================
# 6. Latest Datasets: Column Comparison
# ============================================
print("\n" + "="*60)
print("  6. Latest Datasets: Column Comparison")
print("="*60)

for name, info in LATEST.items():
    print(f"\n{'='*40}")
    print(f"  {name}: {info['description']}")
    try:
        df_tmp = load_sumstats(info["path"], nrows=5, verbose=True)
        col_map_tmp = identify_columns(df_tmp)
        del df_tmp
    except Exception as e:
        print(f"  ERROR loading: {e}")

# ============================================
# 7. Summary of Findings
# ============================================
print("\n" + "="*60)
print("  7. Key Findings")
print("="*60)
print("""
  CDG2019:
    - 6.74M SNPs after cleaning (removed ~15K malformed separator rows)
    - P-values truncated at ~1e-6 in this public release
    - M-value columns present but entirely NA
    - Still useful for LDSC genetic correlation analysis

  CDG2025:
    - Hits file contains significant pleiotropic loci for clustering
    - Factor files (F1-F5 + PFactor) contain Genomic SEM results
    - These are the primary data for our novel ML analysis

  Next Steps:
    - Use individual disorder GWAS files for LDSC genetic correlations
    - Use CDG2025 hits file for locus-level clustering analysis
    - Harmonize column names across all files before pipeline begins
""")

print("\n  Hits per factor:")
print(df_hits["gwas_name"].value_counts().to_string())
print(f"\n  P-value range per factor:")
for factor in sorted(df_hits["gwas_name"].unique()):
    subset = df_hits[df_hits["gwas_name"] == factor]
    print(f"    {factor}: {len(subset)} hits, P range [{subset['P'].min():.2e}, {subset['P'].max():.2e}]")

if __name__ == "__main__":
    sys.stdout.close()
    sys.stdout = sys.__stdout__
    print("Done — output saved to data/results/exploration_output.txt")
    print("Phase 1: Data Exploration Complete")