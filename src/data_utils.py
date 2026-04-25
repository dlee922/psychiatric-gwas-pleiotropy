"""
Utility functions for loading and exploring GWAS summary statistics.
Handles different file formats, separators, and column naming conventions.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from pathlib import Path

def detect_separator(filepath, n_lines=5):
    """
    Detect whether a file is tab-separated, whitespace-separated, or comma-separated
    by reading the first few lines. Handles VCF-format files with ## headers.
    """
    filepath = Path(filepath)

    if filepath.suffix == ".gz":
        import gzip
        opener = gzip.open
    else:
        opener = open

    with opener(filepath, "rt") as f:
        # Skip ## header lines (VCF format)
        for line in f:
            if not line.startswith("##"):
                header = line
                break

    if "\t" in header:
        return "\t"
    elif "," in header:
        return ","
    else:
        return r"\s+"
        # Whitespace-separated
        return r"\s+"


def load_sumstats(filepath, nrows=None, usecols=None, verbose=True):
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    compression = "gzip" if filepath.suffix == ".gz" else "infer"
    sep = detect_separator(filepath)

    # Count ## header lines to skip (VCF format)
    skip_rows = 0
    if filepath.suffix == ".gz":
        import gzip
        opener = gzip.open
    else:
        opener = open

    with opener(filepath, "rt") as f:
        for line in f:
            if line.startswith("##"):
                skip_rows += 1
            else:
                break

    if verbose:
        print(f"Loading: {filepath.name}")
        print(f"  Separator: {'tab' if sep == chr(9) else 'whitespace' if sep == r's+' else sep}")
        print(f"  Compression: {compression}")
        if skip_rows > 0:
            print(f"  Skipping {skip_rows} VCF header lines")
        if nrows:
            print(f"  Rows: first {nrows}")

    start = time.time()

    df = pd.read_csv(
        filepath,
        sep=sep,
        compression=compression,
        nrows=nrows,
        usecols=usecols,
        skiprows=skip_rows,
        na_values=["NA", ".", ""],
    )

    elapsed = time.time() - start

    if verbose:
        print(f"  Loaded: {df.shape[0]:,} rows x {df.shape[1]} columns in {elapsed:.1f}s")

    return df


def quick_look(df, name="dataset"):
    """
    Print a quick summary of a GWAS summary statistics DataFrame.
    """
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    print(f"  Shape: {df.shape[0]:,} rows x {df.shape[1]} columns")
    print(f"  Columns: {df.columns.tolist()}")
    print(f"  Dtypes:\n{df.dtypes.to_string()}")
    print(f"\n  Missing values:")
    missing = df.isnull().sum()
    if missing.sum() > 0:
        print(missing[missing > 0].to_string())
    else:
        print("    None")
    print(f"\n  First 3 rows:")
    print(df.head(3).to_string())
    print(f"{'='*60}\n")


def find_column(df, candidates):
    """
    Find a column by trying multiple possible names.
    GWAS files are inconsistent — SNP might be called 'SNP', 'MarkerName', 'ID', 'rsid', etc.

    Parameters
    ----------
    df : pd.DataFrame
    candidates : list of str
        Possible column names to try (case-insensitive)

    Returns
    -------
    str or None
        The matching column name, or None if not found
    """
    col_lower = {c.lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in col_lower:
            return col_lower[candidate.lower()]
    return None


# Common column name mappings across GWAS formats
COLUMN_CANDIDATES = {
    "snp": ["SNP", "MarkerName", "ID", "rsid", "RSID", "rsID", "variant_id", "SNPID"],
    "chr": ["CHR", "chr", "chromosome", "CHROM", "#CHROM", "hg19chrc"],
    "bp": ["BP", "bp", "pos", "POS", "position", "base_pair_location"],
    "a1": ["A1", "a1", "allele1", "effect_allele", "ALT", "Allele1"],
    "a2": ["A2", "a2", "allele2", "other_allele", "REF", "Allele2"],
    "beta": ["BETA", "beta", "b", "Effect", "effect", "EFFECT"],
    "or": ["OR", "or", "odds_ratio"],
    "se": ["SE", "se", "stderr", "StdErr", "standard_error"],
    "p": ["P", "p", "pval", "p_value", "P-value", "PVAL", "p.value"],
    "z": ["Z", "z", "zscore", "Z_OR", "Zscore"],
    "freq": ["FRQ", "MAF", "freq", "Freq", "FRQ_A_*", "EAF", "frequency"],
    "info": ["INFO", "info", "imputation_quality"],
    "n": ["N", "n", "NMISS", "n_total", "Neff"],
    "nca": ["Nca", "nca", "N_cases", "Ncas", "n_cases"],
    "nco": ["Nco", "nco", "N_controls", "Ncon", "n_controls"],
}


def identify_columns(df, verbose=True):
    """
    Attempt to identify standard GWAS columns regardless of naming convention.

    Returns
    -------
    dict
        Mapping of standard name -> actual column name in this DataFrame
    """
    found = {}
    for standard_name, candidates in COLUMN_CANDIDATES.items():
        match = find_column(df, candidates)
        if match:
            found[standard_name] = match

    if verbose:
        print("Column mapping:")
        for standard, actual in found.items():
            print(f"  {standard:>6} -> {actual}")

        # Flag unmapped columns
        mapped_cols = set(found.values())
        unmapped = [c for c in df.columns if c not in mapped_cols]
        if unmapped:
            print(f"  Unmapped columns: {unmapped}")

    return found


def compute_zscore(df, col_map):
    """
    Compute z-scores from available columns.
    Tries BETA/SE first, falls back to log(OR)/SE.

    Parameters
    ----------
    df : pd.DataFrame
    col_map : dict
        Output of identify_columns()

    Returns
    -------
    pd.Series
        Z-scores
    """
    if "z" in col_map:
        print("  Z-scores already present in data")
        return df[col_map["z"]]
    elif "beta" in col_map and "se" in col_map:
        print("  Computing Z = BETA / SE")
        return df[col_map["beta"]] / df[col_map["se"]]
    elif "or" in col_map and "se" in col_map:
        print("  Computing Z = log(OR) / SE")
        return np.log(df[col_map["or"]]) / df[col_map["se"]]
    else:
        print("  WARNING: Cannot compute z-scores — missing required columns")
        return None


def gwas_summary(df, col_map, name="dataset"):
    """
    Print a GWAS-specific summary: SNP counts, significance, p-value range, etc.
    """
    print(f"\n{'='*60}")
    print(f"  GWAS Summary: {name}")
    print(f"{'='*60}")
    print(f"  Total SNPs: {len(df):,}")

    # if "chr" in col_map:
    #     chroms = df[col_map['chr']].dropna().unique()
    #     print(f"  Chromosomes: {sorted(chroms, key=str)}")

    if "p" in col_map:
        p = df[col_map["p"]]
        print(f"  P-value range: [{p.min():.2e}, {p.max():.2e}]")
        n_sig = (p < 5e-8).sum()
        print(f"  Genome-wide significant (P < 5e-8): {n_sig:,}")
        n_suggestive = ((p >= 5e-8) & (p < 1e-5)).sum()
        print(f"  Suggestive (1e-5 < P < 5e-8): {n_suggestive:,}")

    if "freq" in col_map:
        freq = df[col_map["freq"]]
        print(f"  Allele frequency range: [{freq.min():.4f}, {freq.max():.4f}]")

    if "info" in col_map:
        info = df[col_map["info"]]
        print(f"  INFO score range: [{info.min():.4f}, {info.max():.4f}]")
        print(f"  SNPs with INFO > 0.9: {(info > 0.9).sum():,}")

    print(f"{'='*60}\n")


def plot_exploration(df, col_map, name="dataset", save_dir=None):
    """
    Generate standard exploration plots for GWAS summary statistics.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Data Exploration: {name}", fontsize=14)

    # 1. P-value histogram
    if "p" in col_map:
        axes[0].hist(df[col_map["p"]], bins=50, edgecolor="black", alpha=0.7)
        axes[0].set_xlabel("P-value")
        axes[0].set_ylabel("Count")
        axes[0].set_title("P-value Distribution")

    # 2. Z-score histogram
    z = compute_zscore(df, col_map)
    if z is not None:
        axes[1].hist(z.dropna(), bins=100, edgecolor="black", alpha=0.7)
        axes[1].set_xlabel("Z-score")
        axes[1].set_ylabel("Count")
        axes[1].set_title("Z-score Distribution")

    # 3. Manhattan-style plot
    if "p" in col_map:
        neg_log_p = -np.log10(df[col_map["p"]])
        axes[2].scatter(range(len(neg_log_p)), neg_log_p, s=0.5, alpha=0.3)
        axes[2].axhline(y=-np.log10(5e-8), color="red", linestyle="--",
                        linewidth=1, label="Genome-wide sig.")
        axes[2].set_xlabel("SNP Index")
        axes[2].set_ylabel("-log10(P)")
        axes[2].set_title("Manhattan Plot (by index)")
        axes[2].legend()

    plt.tight_layout()

    if save_dir:
        save_path = Path(save_dir) / f"exploration_{name}.png"
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_path}")

    plt.show()


def explore_dataset(filepath, name="dataset", nrows=None, save_dir=None):
    """
    One-call full exploration of any GWAS summary statistics file.
    Load, identify columns, print summary, generate plots.

    Parameters
    ----------
    filepath : str or Path
    name : str
        Label for prints and plot titles
    nrows : int, optional
        Limit rows (useful for quick peeks)
    save_dir : str or Path, optional
        Directory to save plots
    """
    df = load_sumstats(filepath, nrows=nrows)
    quick_look(df, name=name)
    col_map = identify_columns(df)
    gwas_summary(df, col_map, name=name)
    plot_exploration(df, col_map, name=name, save_dir=save_dir)
    return df, col_map

def load_cdg2019(filepath, verbose=True):
    """
    Special loader for the CDG2019 cross-disorder file.
    This file has chromosome separator lines that corrupt standard parsing.
    We filter them out after loading.
    """
    df = load_sumstats(filepath, verbose=verbose)
    
    original_len = len(df)
    
    # Drop rows where core numeric columns are NaN (these are separator lines)
    df = df.dropna(subset=["PVAL", "BETA", "SE"])
    
    # Convert CHROM to clean integer values
    # Some CHROM values are full data lines that got misparsed
    df = df[df["CHROM"].apply(lambda x: len(str(x)) <= 2)]
    df["CHROM"] = df["CHROM"].astype(int)
    
    if verbose:
        print(f"  Cleaned: removed {original_len - len(df):,} malformed rows")
        print(f"  Final shape: {df.shape[0]:,} rows x {df.shape[1]} columns")
    
    return df
