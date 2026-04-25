"""
Central configuration for all dataset paths and metadata.
"""
from pathlib import Path

# Root data directory — adjust if your scripts live in different locations
DATA_DIR = Path(__file__).parent.parent / "data" / "raw"

# ============================================
# Cross-Disorder Datasets
# ============================================
CROSS_DISORDER = {
    "cdg2019": {
        "path": DATA_DIR / "cross_disorder" / "cdg2019" / "pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt",
        "description": "8-disorder cross-disorder meta-analysis (Lee et al., 2019)",
    },
    "cdg2025_F1_compulsive": {
        "path": DATA_DIR / "cross_disorder" / "cdg2025" / "F1_CompulsiveDisorders_2025.tsv.gz",
        "description": "CDG3 Factor 1: Compulsive Disorders",
    },
    "cdg2025_F2_scz_bip": {
        "path": DATA_DIR / "cross_disorder" / "cdg2025" / "F2_SchizophreniaBipolar_2025.tsv.gz",
        "description": "CDG3 Factor 2: Schizophrenia/Bipolar",
    },
    "cdg2025_F3_neurodev": {
        "path": DATA_DIR / "cross_disorder" / "cdg2025" / "F3_Neurodevelopmental_2025.tsv.gz",
        "description": "CDG3 Factor 3: Neurodevelopmental",
    },
    "cdg2025_F4_internalizing": {
        "path": DATA_DIR / "cross_disorder" / "cdg2025" / "F4_Internalizing_2025.tsv.gz",
        "description": "CDG3 Factor 4: Internalizing",
    },
    "cdg2025_F5_substance": {
        "path": DATA_DIR / "cross_disorder" / "cdg2025" / "F5_SubstanceUse_2025.tsv.gz",
        "description": "CDG3 Factor 5: Substance Use",
    },
    "cdg2025_pfactor": {
        "path": DATA_DIR / "cross_disorder" / "cdg2025" / "PFactor_2025.tsv.gz",
        "description": "CDG3 Hierarchical P-Factor",
    },
    "cdg2025_hits": {
        "path": DATA_DIR / "cross_disorder" / "cdg2025" / "CDG3_Hits_2025.txt",
        "description": "CDG3 significant loci hit list",
    },
}

# ============================================
# Replication Datasets (matching 2019 paper)
# ============================================
REPLICATION = {
    # "adhd2019": { # mac version
    #     "path": DATA_DIR / "replication" / "adhd2019" / "daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz",
    #     "description": "ADHD GWAS (Demontis et al., 2019)",
    # },
    "adhd2019": { 
        "path": DATA_DIR / "replication" / "adhd2019" / "daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta",
        "description": "ADHD GWAS (Demontis et al., 2019)",
    },
    "an2019": {
        "path": DATA_DIR / "replication" / "an2019" / "pgcAN2.2019-07.vcf.tsv.gz",
        "description": "Anorexia Nervosa GWAS (Watson et al., 2019)",
    },
    "asd2019": {
        "path": DATA_DIR / "replication" / "asd2019" / "iPSYCH-PGC_ASD_Nov2017.gz",
        "description": "Autism Spectrum Disorder GWAS (Grove et al., 2019)",
    },
    "bip2019": {
        "path": DATA_DIR / "replication" / "bip2019" / "daner_PGC_BIP32b_mds7a_0416a.gz",
        "description": "Bipolar Disorder GWAS (Stahl et al., 2019)",
    },
    "mdd2018": {
        "path": DATA_DIR / "replication" / "mdd2018" / "MDD2018_ex23andMe.gz",
        "description": "Major Depressive Disorder GWAS (Wray et al., 2018)",
    },
    "ocd2018": {
        "path": DATA_DIR / "replication" / "ocd2018" / "ocd_aug2017.gz",
        "description": "OCD GWAS (IOCDF-GC & OCGAS, 2018)",
    },
    # "scz2018": { # use this one for mac
    #     "path": DATA_DIR / "replication" / "scz2018" / "CLOZUK_PGC2noclo.METAL.assoc.dosage.fix.gz",
    #     "description": "Schizophrenia GWAS (Pardiñas et al., 2018)",
    # },
    "scz2018": { # use this one for windows. windows auto extract the gz file
        "path": DATA_DIR / "replication" / "scz2018" / "CLOZUK_PGC2noclo.METAL.assoc.dosage.fix",
        "description": "Schizophrenia GWAS (Pardiñas et al., 2018)",
    },
    "ts2019": {
        "path": DATA_DIR / "replication" / "ts2019" / "TS_Oct2018.gz",
        "description": "Tourette Syndrome GWAS (Yu et al., 2019)",
    },
}

# ============================================
# Latest Datasets (most recent GWAS)
# ============================================
LATEST = {
    "adhd2022": {
        "path": DATA_DIR / "latest" / "adhd2022" / "ADHD2022_iPSYCH_deCODE_PGC.meta.gz",
        "description": "ADHD GWAS (Demontis et al., 2022)",
    },
    "an2019": {
        "path": DATA_DIR / "latest" / "an2019" / "pgcAN2.2019-07.vcf.tsv.gz",
        "description": "Anorexia Nervosa GWAS (Watson et al., 2019)",
    },
    "asd2019": {
        "path": DATA_DIR / "latest" / "asd2019" / "iPSYCH-PGC_ASD_Nov2017.gz",
        "description": "Autism Spectrum Disorder GWAS (Grove et al., 2019)",
    },
    "bip2024": {
        "path": DATA_DIR / "latest" / "bip2024" / "bip2024_eur_no23andMe.gz",
        "description": "Bipolar Disorder GWAS (Mullins et al., 2024)",
    },
    "mdd2025": {
        "path": DATA_DIR / "latest" / "mdd2025" / "pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz",
        "description": "Major Depressive Disorder GWAS (PGC MDD, 2025)",
    },
    "ocd2025": {
        "path": DATA_DIR / "latest" / "ocd2025" / "daner_OCD_full_wo23andMe_190522.gz",
        "description": "OCD GWAS (PGC OCD, 2025)",
    },
    "scz2022": {
        "path": DATA_DIR / "latest" / "scz2022" / "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",  # TODO: verify filename
        "description": "Schizophrenia GWAS (Trubetskoy et al., 2022)",
    },
    "ts2019": {
        "path": DATA_DIR / "latest" / "ts2019" / "TS_Oct2018.gz",
        "description": "Tourette Syndrome GWAS (Yu et al., 2019)",
    },
}

# Convenience: all datasets in one dict
ALL_DATASETS = {
    **{f"cross_{k}": v for k, v in CROSS_DISORDER.items()},
    **{f"rep_{k}": v for k, v in REPLICATION.items()},
    **{f"latest_{k}": v for k, v in LATEST.items()},
}