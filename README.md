# Cross-Disorder Psychiatric GWAS Clustering

Locus-level clustering of pleiotropic genetic effects across five psychiatric
genomic factors derived from the Grotzinger et al. (2026) cross-disorder
Genomic SEM framework, using summary statistics from the Psychiatric Genomics
Consortium (PGC).

## Motivation

Grotzinger et al. (Nature, 2026) applied Genomic SEM to GWAS summary statistics
spanning 14 psychiatric disorders and identified five latent genomic factors:
compulsive (F1), schizophrenia-bipolar (F2), neurodevelopmental (F3),
internalizing (F4), and substance use (F5). Their analysis characterized
genome-wide significant loci on a factor-by-factor basis.

This project asks a different question: across all 611 genome-wide significant
loci, do the cross-factor effect profiles organize into distinct subtypes of
pleiotropy? Rather than examining which loci are significant for each factor
individually, this analysis clusters loci by their full five-dimensional
z-score vectors to identify structured patterns of shared and opposing genetic
influence across psychiatric domains.

## Data

All summary statistics were obtained from the
[PGC](https://pgc.unc.edu/for-researchers/download-results/) under their data
use terms. Raw data files are not included in this repository.

The primary analysis uses factor-level GWAS summary statistics from the 14-disorder
Genomic SEM model (Grotzinger et al., 2026), which provides effect estimates for
each locus on each of the five genomic factors. The input set consists of 611 loci
reaching genome-wide significance (p < 5e-8) on at least one factor or the
hierarchical p-factor.

## Methods

### Effect size standardization

Z-scores (beta / SE) served as the primary effect size metric. Z-scores are
field-standard for cross-GWAS comparisons because they normalize for differences
in sample size and allele frequency across disorder-specific GWAS. All loci were
filtered to INFO > 0.9, and SNPs were intersected across factor-level datasets
to ensure consistent allele harmonization.

### Clustering

K-means clustering (50 random initializations, Euclidean distance) was the
primary analysis. The number of clusters was selected by evaluating k = 2
through 10 using silhouette scores. Three comparison methods were run at
matched k values to assess robustness: hierarchical agglomerative clustering
(Ward linkage), Gaussian mixture models (evaluated by BIC and silhouette),
and DBSCAN as a density-based sanity check. All three comparison methods
recovered qualitatively similar cluster profiles at k = 3.

k = 3 was selected as the primary solution (silhouette = 0.256, nearly tied
with k = 2 at 0.254) based on biological interpretability, as k = 3 revealed
a distinct antagonistic pleiotropy cluster not visible at k = 2. k = 4 was
examined as a secondary solution to assess whether additional structure was
present.

### Validation

Clusters were validated using Q_P heterogeneity statistics from the Genomic SEM
output. Q_P functions as a p-value for a heterogeneity test, where lower Q_P
indicates greater deviation from the factor model's predicted effect pattern,
with Q_P < 0.05 indicating statistically significant heterogeneity. Q_P values
were used for post-hoc cluster annotation, not as clustering features.

Additionally, loci flagged as significant Q hits (QSNP p < 5e-8) in the original
CDG2025 results were tested for enrichment across clusters.

### Pathway enrichment

SNPs in each cluster were mapped to genes via g:Profiler's SNP-to-gene conversion.
Per-cluster enrichment was tested against Gene Ontology Biological Process (GO:BP)
and KEGG pathway databases, using the full set of 484 mapped genes as the
background.

## Results

### k = 3 cluster profiles

| Cluster | n loci | F1 Compulsive | F2 SCZ-BIP | F3 Neurodev | F4 Internalizing | F5 Substance |
|---------|--------|---------------|------------|-------------|------------------|--------------|
| C0      | 355    | 1.24          | 2.82       | 1.30        | 6.23             | 2.20         |
| C1      | 182    | 1.54          | 5.98       | -0.09       | 2.42             | 1.64         |
| C2      | 74     | -0.51         | -3.50      | 0.06        | 2.26             | -0.77        |

Values are mean z-scores per cluster.

**Cluster 0 (355 loci): Broadly pleiotropic, internalizing-dominant.** Positive
mean z-scores across all five factors, with the strongest signal on F4
(internalizing). This cluster captures loci with similar effects across
psychiatric domains.

**Cluster 1 (182 loci): Psychotic-spectrum specific.** Dominated by F2
(schizophrenia-bipolar, mean z = 5.98) with negligible neurodevelopmental signal
and moderate internalizing effects.

**Cluster 2 (74 loci): Antagonistic pleiotropy.** Positive internalizing (F4)
effects paired with negative schizophrenia-bipolar (F2), compulsive (F1), and
substance use (F5) effects. These loci increase risk for internalizing disorders
while decreasing risk for psychotic-spectrum and externalizing conditions.

### Q_P validation

Clusters showed a monotonic Q_P gradient across all five factors (Kruskal-Wallis
p < 0.001 for each): C0 had the highest mean Q_P (~0.45), C2 the lowest (~0.27).
Enrichment of QSNP-significant hits (_Q hits) across clusters was significant
(chi-squared p = 7e-54): C0 = 7.3%, C1 = 59.3%, C2 = 81.1%.

This gradient indicates that the clustering is capturing meaningful structure in
how loci relate to the factor model. C2 (antagonistic pleiotropy) contains the
highest concentration of loci whose effects deviate most from the Genomic SEM
factor predictions.

### Pathway enrichment

Cluster 2 was the only cluster with significant pathway enrichment:

| Term | Source | Adjusted p | Genes |
|------|--------|-----------|-------|
| Homophilic cell-cell adhesion | GO:BP | 1.26e-05 | 13 |
| Cadherin signaling pathway | KEGG | 1.35e-05 | 13 |
| Cell-cell adhesion | GO:BP | 1.37e-03 | 17 |
| Cell adhesion | GO:BP | 1.07e-02 | 18 |

Clusters 0 and 1 showed no significant enrichment against the within-study
gene background.

### k = 4 secondary analysis

Splitting to k = 4 subdivides the broad C0 cluster into an
internalizing-substance subgroup and an internalizing-neurodevelopmental
subgroup. The antagonistic pleiotropy cluster (C2) and psychotic-spectrum
cluster (C1) remain stable across k = 3 and k = 4.

## Project structure

    psychiatric-gwas-pleiotropy/
    ├── data/
    │   ├── raw/                      # PGC summary statistics (gitignored)
    │   │   ├── cross_disorder/
    │   │   │   ├── cdg2019/
    │   │   │   └── cdg2025/
    │   │   ├── replication/          # GWAS versions matching CDG2019
    │   │   └── latest/              # Most recent GWAS per disorder
    │   ├── processed/                # Harmonized and cleaned data
    │   └── results/                  # Analysis outputs
    ├── notebooks/                    # Jupyter notebooks
    ├── scripts/                      # Analysis scripts
    ├── src/                          # Reusable Python modules
    │   ├── config.py
    │   ├── data_processing.py
    │   ├── clustering.py
    │   ├── enrichment.py
    │   └── visualization.py
    ├── figures/
    ├── environment.yml
    ├── .gitignore
    └── README.md

## Reproducibility

### Prerequisites

- [Anaconda](https://www.anaconda.com/download) or
  [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Git
- PGC data access (see Data section above)

### Setup

```bash
git clone https://github.com/dlee922/psychiatric-gwas-pleiotropy.git
cd psychiatric-gwas-pleiotropy
conda env create -f environment.yml
conda activate psych-gwas
```

Raw summary statistics must be downloaded separately from PGC and placed in
`data/raw/` following the directory structure above.

## References

- Grotzinger, A.D., et al. (2026). Mapping the genetic landscape across
  14 psychiatric disorders. *Nature*, 649, 406-415.
  https://doi.org/10.1038/s41586-025-09820-3
- Lee, P.H., et al. (2019). Genomic relationships, novel loci, and pleiotropic
  mechanisms across eight psychiatric disorders. *Cell*, 179(7), 1469-1482.
- Grotzinger, A.D., et al. (2019). Genomic structural equation modelling provides
  insights into the multivariate genetic architecture of complex traits.
  *Nature Human Behaviour*, 3, 513-525.
- Marees, A.T., et al. (2018). A tutorial on conducting genome-wide association
  studies. *International Journal of Methods in Psychiatric Research*, 27(2), e1608.

## Status

This project is under active development. Planned extensions include LD score
regression analysis and individual disorder-level sensitivity analyses.

## License

For scientific research and educational purposes only, in accordance with PGC
data use terms.