# HNSCC scRNA-seq Analysis | HRPS and HOPS Phenotypes

**An R-based collection of exploratory workflows for investigating phenotypic heterogeneity, transcriptional programs, and cell–cell interactions in head and neck squamous cell carcinoma (HNSCC).**

This repository contains scripts for joint single-cell visualization, cluster characterization, sample-level differential expression, inferred ligand–receptor signaling, immune-cell fraction comparisons, and evaluation of an HOPS-associated gene-expression signature in an external dataset.

> **Scope:** These scripts represent separate analysis modules, **not** a fully automated end-to-end pipeline. They require user-supplied Seurat objects, sample metadata, or external count files as appropriate. No biological findings, performance claims, or clinically validated biomarkers are established by the repository alone.

## Analyses at a glance

| Script | Purpose | Main input | Main results |
|---|---|---|---|
| [`CombinedUMAP.R`](CombinedUMAP.R) | Merge HRPS and HOPS Seurat objects; run UMAP and unsupervised clustering; summarize cluster composition and marker expression | Two Seurat `.rds` objects | Combined Seurat object, UMAP, cluster proportions, cluster markers |
| [`Deseq2hypervshypo.R`](Deseq2hypervshypo.R) | Compare HRPS and HOPS gene expression using **biological-sample-level pseudobulk** counts | Two Seurat `.rds` objects with sample IDs | DESeq2 results, upregulated gene lists, volcano plot, heatmap |
| [`cellchat.R`](cellchat.R) | Infer putative ligand–receptor communication between cell groups | One normalized Seurat `.rds` object | Predicted interaction tables, count/weight matrices, network plots |
| [`cibersortx1.R`](cibersortx1.R) | Compare CIBERSORTx LM22 immune fractions between two HOPS/HRPS phenotype categories | CIBERSORTx results `.csv` and clinical metadata `.csv` | Group summaries, immune statistics with BH adjustment, plots |
| [`drugsensitivity3.R`](drugsensitivity3.R) | Evaluate a pre-defined HOPS expression signature in external GEO/HTSeq-like count data | Per-sample raw count files, metadata, gene list | HOPS signature scores and plot |

**Terminology:** The scripts use *HRPS*, *HOPS*, and phenotype labels such as `HOPS (Specialist)` and `HRPS (Architect)`. Their biological definitions and derivation should be documented alongside the associated study before interpreting these labels. In the deconvolution script, the precise input categories are `High HOPS / Low HRPS` and `High HRPS / Low HOPS`.

## Repository structure

```text
HNSCC_SCrna/
├── README.md
├── CombinedUMAP.R
├── Deseq2hypervshypo.R
├── cellchat.R
├── cibersortx1.R
└── drugsensitivity3.R
```

Input examples below use a local `data/` directory and save generated files under `results/`. These directories and input files are **not asserted to be included in this repository**.

## Requirements and installation

- R with `Rscript` available from the terminal.
- [Seurat](https://satijalab.org/seurat/) objects with an `RNA` assay for the single-cell workflows.
- Sufficient memory for sparse matrices, clustering, and (when applicable) dense pseudobulk count conversion.
- The correct biological sample identifiers and experimental metadata for statistical comparisons.

Install CRAN packages used across the scripts:

```r
install.packages(c(
  "Seurat", "ggplot2", "patchwork", "Matrix", "ggrepel",
  "pheatmap", "dplyr", "tidyr"
))
```

Install Bioconductor dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi"))
```

Install [CellChat](https://github.com/jinworks/CellChat) according to its **official installation instructions** and verify compatibility with your version of R and Seurat. The UMAP script can optionally use **MAST** for marker testing; otherwise it uses the Wilcoxon test:

```r
BiocManager::install("MAST")  # optional
```

To record software versions after an analysis:

```r
writeLines(capture.output(sessionInfo()), "session_info.txt")
```

> The repository does not currently supply an environment lockfile or tested version matrix. Package compatibility and successful execution must be checked in the user's R environment.

## Getting started

Clone the repository and run commands from its root directory:

```bash
git clone https://github.com/Dyumn256/HNSCC_SCrna.git
cd HNSCC_SCrna
mkdir -p results
```

The following commands are examples. Replace file paths and metadata column names with those of your dataset.

### 1. Joint UMAP and cluster analysis — `CombinedUMAP.R`

```bash
Rscript CombinedUMAP.R data/HRPS.rds data/HOPS.rds results/umap
```

**Inputs:** Two Seurat objects containing an `RNA` assay. The first is assigned the label `HRPS` and the second `HOPS`. The script merges them, joins RNA layers where supported, normalizes expression, selects variable features, scales, runs PCA, neighbor finding, clustering, and UMAP.

**Outputs:**

- `combined_seurat.rds` — merged and processed Seurat object.
- `umap.pdf` — condition and cluster UMAP panels.
- `cluster_composition.csv` and `cluster_composition.pdf` — cell proportions per cluster and per phenotype.
- `cluster_markers.csv` — positive cluster markers.
- `selected_gene_dotplot.pdf` — selected genes, when available in the input matrix.
- `analysis_notes.txt` — runtime notes and caveats.

**Interpretation:** The script does **not** apply batch integration. UMAP separation may reflect sequencing batch or other technical effects; cell-level marker screens do not provide donor-level inference.

### 2. Pseudobulk differential expression — `Deseq2hypervshypo.R`

```bash
Rscript Deseq2hypervshypo.R \
  data/HRPS.rds data/HOPS.rds results/deseq donor_id
```

**Inputs:** Two Seurat objects with raw `RNA` counts and a metadata column (here `donor_id`) identifying **independent biological samples** in **both** objects. The script needs at least **two independent samples per phenotype**; more are preferable. IDs must not overlap across conditions as currently implemented.

The script sums cell-level raw counts within each sample and applies DESeq2 to the resulting pseudobulk expression profiles. The contrast is explicitly **HRPS versus HOPS**:

- `log2FoldChange > 0`: higher expression in HRPS.
- `log2FoldChange < 0`: higher expression in HOPS.
- The script exports upregulated genes using **adjusted p-value < 0.05** and **absolute log2 fold change > 0.5** in the appropriate direction.

**Outputs:** `all_DESeq2_results.csv`, `HRPS_upregulated.csv`, `HOPS_upregulated.csv`, `HRPS_genes.txt`, `HOPS_genes.txt`, `sample_metadata.csv`, `volcano.pdf`, and a heatmap when sufficient genes are available.

**Interpretation:** Individual cells are not independent biological replicates. The current model uses `~ condition`; it does not explicitly adjust for donor pairing, batch, or other covariates. Confirm that this design is appropriate before interpreting significance.

### 3. Cell–cell communication inference — `cellchat.R`

```bash
Rscript cellchat.R data/HRPS.rds results/cellchat human cell_type
```

The fourth argument is optional; omit it to use the Seurat object's active identities:

```bash
Rscript cellchat.R data/HRPS.rds results/cellchat human
```

**Inputs:** One Seurat object with normalized data in the `RNA` assay and meaningful cell-group labels. Choose `human` or `mouse` to match the dataset's gene identifiers. A `cell_type` metadata column is only an example.

The script uses the species-specific CellChat database subset **Secreted Signaling**, computes predicted communications, applies a minimum cell-group size filter, and aggregates pathways and network weights.

**Outputs:** `cellchat_results.rds`, `cell_group_sizes.csv`, `interaction_count_matrix.csv`, `interaction_weight_matrix.csv`, `ligand_receptor_interactions.csv`, `detected_pathways.txt`, `overview_network.pdf`, and (if a pathway is detected) `example_pathway.pdf`, plus `analysis_notes.txt`.

**Interpretation:** CellChat interactions are computational predictions based on ligand/receptor expression, not direct experimental evidence of intercellular signaling. Results depend on cluster definitions, group sizes, expression preprocessing, and the selected database.

### 4. Immune-cell fractions — `cibersortx1.R`

```bash
Rscript cibersortx1.R \
  data/CIBERSORTx_results.csv data/clinical.csv results/cibersort
```

**Inputs:**

- CIBERSORTx CSV with `Mixture` sample identifiers and all **22 canonical LM22 immune-cell fraction columns**.
- Clinical metadata CSV with `SampleID` and `Clonal_Phenotype` columns.
- Exact sample-ID matches between `Mixture` and `SampleID`.

Example metadata structure:

```csv
SampleID,Clonal_Phenotype
sample_01,High HOPS / Low HRPS
sample_02,High HRPS / Low HOPS
```

These rows are illustrative only. The script requires at least two matched samples in each comparison group.

The intended analysis summarizes fractions per immune-cell category, performs Wilcoxon rank-sum tests, and applies Benjamini–Hochberg correction across the 22 comparisons.

**Expected outputs:** `sample_matching_audit.csv`, `group_sizes.csv`, `immune_statistics.csv`, `immune_boxplots.pdf`, and an immune-fraction bar plot.

**Known issue in the currently published code:** `cibersortx1.R` has a malformed statement at approximately line 27 (`... )ndat <- ...`) that causes a syntax error. **Correct this script separately before trying to run it. This README does not modify the R files.**

**Interpretation:** Immune fractions are compositional estimates and may be affected by deconvolution quality. The script does not automatically filter CIBERSORTx results by its fit statistics; plan any such filter before comparative testing.

### 5. External HOPS gene-signature evaluation — `drugsensitivity3.R`

```bash
Rscript drugsensitivity3.R \
  data/GSE200421_RAW data/fadu_metadata.csv \
  data/HOPS_genes.txt results/external_signature
```

**Inputs:**

1. A directory with one count file per sample (`.txt`, `.tsv`, and optionally `.gz`). Filenames must start with a unique GEO `GSM` identifier.
2. Each count file must contain **two tab-separated, headerless columns**: Ensembl gene identifier and raw integer count.
3. A metadata CSV with `ID` and `Group` columns, where IDs match count-file GSM identifiers.
4. A text file containing **at least two independently defined HGNC gene symbols**, one per line.

Example metadata:

```csv
ID,Group
GSM1234567,control 6h
GSM1234568,EGF 9nM 6h
```

The script maps human Ensembl IDs to gene symbols and computes the average **log2(CPM + 1)** expression across available signature genes for each sample.

**Outputs:** `signature_scores.csv`, `signature_scores.pdf`, `matched_signature_genes.txt`, `unmatched_signature_genes.txt`, and analysis notes.

**Important:** Despite its retained historical filename, `drugsensitivity3.R` **does not predict drug response, drug sensitivity, IC50, or treatment resistance**. It performs an exploratory expression signature analysis. The filename `GSE200421_RAW` is an example of the external count directory used in the script; verify data provenance, sample inclusion, and metadata before use.

## Suggested analysis order

```text
Prepared Seurat objects for HRPS and HOPS
      |
      +--> CombinedUMAP.R ----------> Joint visualization, clusters, markers
      |
      +--> Deseq2hypervshypo.R -----> Sample-level differential expression
      |
      +--> cellchat.R ---------------> Predicted cell–cell communication

CIBERSORTx results + matched clinical metadata
      |
      +--> cibersortx1.R -----------> LM22 immune fraction comparisons

External HTSeq counts + sample metadata + pre-defined HOPS genes
      |
      +--> drugsensitivity3.R ------> External expression signature scores
```

These branches are **not automatically chained**. Each requires its own appropriate input and biological interpretation.

## Reproducibility and quality checks

Before drawing biological conclusions:

1. Confirm that the input species, gene identifiers, RNA assay, and normalized/raw count layers match each tool's expectations.
2. Record cell and sample inclusion criteria, QC thresholds, batch structure, and biological replicate counts.
3. Confirm the `HRPS` and `HOPS` labels are mapped consistently between objects, metadata, plots, and outputs.
4. Inspect cluster composition and sample-level representation before interpreting pooled-cell results.
5. Report **effect sizes**, sample counts, statistical methods, and multiple-testing correction for differential analyses.
6. Review output warnings and `analysis_notes.txt` files; preserve `sessionInfo()` and input provenance for reproducibility.
7. Validate gene signature derivation independently of any external validation dataset to avoid data leakage.

## Data availability and privacy

The analysis code is provided here; the presence of scripts does **not** imply that patient-level data, processed Seurat objects, CIBERSORTx results, or GEO files are distributed with this repository. Check dataset licenses and data-use restrictions before sharing any biological or clinical metadata. Do not commit sensitive patient identifiers or confidential files.

## Citation and acknowledgements

If you use these scripts in a manuscript, cite the original methods and software as relevant: [Seurat](https://satijalab.org/seurat/), [DESeq2](https://bioconductor.org/packages/DESeq2/), [CellChat](https://github.com/jinworks/CellChat), [CIBERSORTx](https://cibersortx.stanford.edu/), and [Bioconductor](https://bioconductor.org/). Also cite the primary source of any HNSCC or external GEO data used in the analysis.

**Repository:** https://github.com/Dyumn256/HNSCC_SCrna

---

*Documentation describes the scripts currently available in the repository. It does not certify that they have passed end-to-end runtime or scientific validation.*
