# R workflows: single-cell communication, immunophenotyping and expression analysis

These five scripts are cleaned, documented adaptations of exploratory R notebooks. Their original filenames are retained for easier replacement in a GitHub repository. **No input datasets or biological results are included.** All examples use placeholder file paths. The scripts have not been runtime-validated in this environment (R is not installed); test with your own input files before interpreting outputs.

## Repository layout

```text
.
├── README.md
├── cellchat.R
├── cibersortx1.R
├── CombinedUMAP.R
├── Deseq2hypervshypo.R
├── drugsensitivity3.R
└── .gitignore
```

## Installation

Install a compatible current R release and required R packages in a fresh environment. Most packages can be installed with:

```r
install.packages(c("Seurat", "ggplot2", "patchwork", "Matrix", "ggrepel", "pheatmap", "dplyr", "tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi"))
```

Install **CellChat** following its official repository's instructions, as it may require GitHub installation and other dependencies. `MAST` is an optional package for cluster marker testing (the UMAP script falls back to a Wilcoxon test). `ggdist` is *not* required for the revised external-signature script. Check version compatibility for Seurat v4 versus v5, especially layered assays.

Record package versions for reproducibility:

```r
writeLines(capture.output(sessionInfo()), "session_info.txt")
```

## Quick-start commands

Run commands from the repository root; replace every example path with real paths.

| Script | Example command | Required input | Main output |
| --- | --- | --- | --- |
| `CombinedUMAP.R` | `Rscript CombinedUMAP.R data/HRPS.rds data/HOPS.rds results/umap` | Two Seurat RNA objects | Merged UMAP, cluster distribution, markers, `.rds` |
| `Deseq2hypervshypo.R` | `Rscript Deseq2hypervshypo.R data/HRPS.rds data/HOPS.rds results/deseq donor_id` | Seurat RNA objects with a biological-sample ID column in **both** objects | Pseudobulk DESeq2 tables, volcano plot, heatmap |
| `cellchat.R` | `Rscript cellchat.R data/HRPS.rds results/cellchat human cell_type` | Normalized Seurat RNA object; optional grouping metadata column | Inferred communication matrices, LR table, plots |
| `cibersortx1.R` | `Rscript cibersortx1.R data/CIBERSORTx_results.csv data/clinical.csv results/cibersort` | CIBERSORTx LM22 output and sample metadata | FDR-adjusted 22-cell comparisons, PDF plots |
| `drugsensitivity3.R` | `Rscript drugsensitivity3.R data/GSE200421_RAW data/fadu_metadata.csv data/HOPS_genes.txt results/external_signature` | Per-sample two-column HTSeq `.txt[.gz]` files, metadata, an independently defined gene list | HOPS expression signature scores and PDF plot |

**Clinical metadata example** (`data/clinical.csv`):

```csv
SampleID,Clonal_Phenotype
sample_01,High HOPS / Low HRPS
sample_02,High HRPS / Low HOPS
```

In practice, supply more than one independent biological sample for each category. IDs under `SampleID` must match `Mixture` in the CIBERSORTx results exactly. All **22 canonical LM22 columns** must be present.

**External metadata example** (`data/fadu_metadata.csv`):

```csv
ID,Group
GSM1234567,control 6h
GSM1234568,EGF 9nM 6h
```

The first column of each external HTSeq file must be the human Ensembl gene ID and the second must be raw integer count data (no header). File names must start with `GSM` followed by digits. `HOPS_genes.txt` must contain one HGNC symbol per line. Do not derive these same signature genes from the external dataset being used for validation.

## Methods, interpretation and limitations

- **UMAP**: Merges and reprocesses both Seurat objects without batch integration. If experimental condition is confounded with sequencing batch, joint UMAP separation is not evidence of a condition effect. Cluster-level marker tests are cell-level exploratory screens; cells are not independent biological replicates.
- **DESeq2**: Aggregates raw RNA UMI counts to one pseudobulk profile per biological sample, then compares HRPS vs HOPS (`log2FC > 0` indicates HRPS higher). Rejects fewer than two independent samples per group. Use real sample IDs; do **not** fabricate replicates by splitting cells. With confounded batch/donor/condition, unadjusted comparisons are not reliable. Report adjusted p-values and effect sizes.
- **CellChat**: Uses a single explicitly selected human or mouse `Secreted Signaling` database, with normalized RNA data. The output describes *predicted* ligand–receptor communication; it does not establish cell-to-cell signaling experimentally. Use the relevant species and appropriate expression matrix. The script plots only a detected pathway.
- **CIBERSORTx**: Recognizes the canonical LM22 panel rather than blindly using columns 2–23. Computes per-cell-type Wilcoxon tests and Benjamini–Hochberg adjusted p-values. Results depend on deconvolution quality and sample independence; fractions are compositional. The script does not remove samples using CIBERSORTx `P-value`/correlation/RMSE; introduce a prespecified quality threshold if your study warrants it.
- **External signature**: Calculates an exploratory mean log2(CPM+1) score from externally specified HOPS genes. Despite its legacy filename, this is **not a drug sensitivity model**; it does not estimate IC50, growth inhibition, response, or resistance. The pipeline does not automatically test 6 h/72 h contrasts because those must follow the actual biological replication and experimental design. Gene annotation is human-only.

## Reproducibility checklist

1. Verify species, input dimensions, raw-versus-normalized counts and metadata joins.
2. Confirm biological replicates and technical batches before statistical comparisons.
3. Record R/package versions, file provenance and GEO accession details.
4. Examine every `analysis_notes.txt`, matching audit and warning message.
5. Validate plots and statistical conclusions against the study design before citing them.

## Privacy and GitHub hygiene

Do not commit patient-level metadata, raw sequencing files, private accession keys, local `.rds` objects, downloaded results or restricted CIBERSORTx outputs. Review data licenses and privacy/consent requirements. Only commit scripts, documentation and **non-sensitive** sample templates.

## Acknowledgement

Workflow components use Seurat, DESeq2, CellChat, CIBERSORTx outputs and Bioconductor annotation resources. Cite the original software publications, specific tool versions, and primary public datasets in any manuscript.
