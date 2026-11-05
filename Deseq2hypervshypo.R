#!/usr/bin/env Rscript
# Deseq2hypervshypo.R — sample-level pseudobulk differential expression.
# Usage: Rscript Deseq2hypervshypo.R HRPS.rds HOPS.rds output_dir sample_id_column
# Both Seurat objects MUST contain a metadata column identifying independent
# biological samples (e.g. donor_id). Each condition needs >=2 samples, ideally >=3.
# Cell counts within one sample are NOT independent DESeq2 replicates.
# Positive log2FC = higher in HRPS; negative = higher in HOPS.
# Required: Seurat, Matrix, DESeq2, ggplot2, ggrepel, pheatmap.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript Deseq2hypervshypo.R HRPS.rds HOPS.rds output_dir sample_id_column")
req <- c("Seurat", "Matrix", "DESeq2", "ggplot2", "ggrepel", "pheatmap")
for (p in req) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(Seurat))
out <- args[3]; id_col <- args[4]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
a <- readRDS(args[1]); b <- readRDS(args[2])
if (!inherits(a, "Seurat") || !inherits(b, "Seurat")) stop("Input RDS must be Seurat objects")
if (!id_col %in% names(a[[]]) || !id_col %in% names(b[[]])) stop("sample_id_column must exist in BOTH Seurat objects")
if (!"RNA" %in% Assays(a) || !"RNA" %in% Assays(b)) stop("Both objects must have an RNA assay")
# Read raw integer counts. For Seurat v5, JoinLayers combines split count layers.
get_counts <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  if ("JoinLayers" %in% getNamespaceExports("Seurat")) obj <- JoinLayers(obj, assay = "RNA")
  Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
}
ca <- get_counts(a); cb <- get_counts(b)
genes <- intersect(rownames(ca), rownames(cb))
if (!length(genes)) stop("No common gene IDs")
ca <- ca[genes, , drop = FALSE]; cb <- cb[genes, , drop = FALSE]
ida <- as.character(a[[]][colnames(ca), id_col]); idb <- as.character(b[[]][colnames(cb), id_col])
if (anyNA(ida) || anyNA(idb) || any(!nzchar(c(ida, idb)))) stop("Missing/empty sample IDs")
# Avoid accidentally combining same-named samples from different conditions.
if (length(intersect(unique(ida), unique(idb)))) stop("Sample IDs overlap between phenotypes; make sample IDs globally unique")
# Construct a sparse cell-to-sample indicator and sum, never average, raw counts.
aggregate_counts <- function(counts, ids) {
  groups <- unique(ids)
  indicator <- Matrix::sparseMatrix(i = seq_along(ids), j = match(ids, groups), x = 1,
                                    dims = c(length(ids), length(groups)), dimnames = list(NULL, groups))
  counts %*% indicator
}
pa <- aggregate_counts(ca, ida); pb <- aggregate_counts(cb, idb)
if (ncol(pa) < 2 || ncol(pb) < 2) stop("DESeq2 requires biological replicates: >=2 independent samples per phenotype")
counts <- cbind(pa, pb)
metadata <- data.frame(condition = factor(c(rep("HRPS", ncol(pa)), rep("HOPS", ncol(pb))), levels = c("HOPS", "HRPS")), row.names = colnames(counts))
if (!identical(rownames(metadata), colnames(counts))) stop("Metadata mismatch")
if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)) stop("RNA counts must be nonnegative integers")
counts <- round(as.matrix(counts))
keep <- rowSums(counts >= 10) >= 2
if (sum(keep) < 2) stop("Too few expressed genes after filtering")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = metadata, design = ~ condition)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("condition", "HRPS", "HOPS"), alpha = 0.05)
res_df <- as.data.frame(res); res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.csv(res_df, file.path(out, "all_DESeq2_results.csv"), row.names = FALSE)
hrps_up <- subset(res_df, !is.na(padj) & padj < 0.05 & log2FoldChange > 0.5)
hops_up <- subset(res_df, !is.na(padj) & padj < 0.05 & log2FoldChange < -0.5)
write.csv(hrps_up, file.path(out, "HRPS_upregulated.csv"), row.names = FALSE)
write.csv(hops_up, file.path(out, "HOPS_upregulated.csv"), row.names = FALSE)
writeLines(hrps_up$gene, file.path(out, "HRPS_genes.txt"))
writeLines(hops_up$gene, file.path(out, "HOPS_genes.txt"))
write.csv(metadata, file.path(out, "sample_metadata.csv"))
res_df$category <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0.5, "HRPS up",
                          ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.5, "HOPS up", "Not significant"))
plot_data <- res_df[is.finite(res_df$log2FoldChange) & !is.na(res_df$padj), ]
plot_data$mlog10padj <- -log10(pmax(plot_data$padj, .Machine$double.xmin))
top <- head(subset(plot_data, category != "Not significant"), 10)
p <- ggplot2::ggplot(plot_data, ggplot2::aes(log2FoldChange, mlog10padj, color = category)) +
  ggplot2::geom_point(alpha = 0.65, size = 0.7) +
  ggplot2::scale_color_manual(values = c("HRPS up" = "#D55E00", "HOPS up" = "#0072B2", "Not significant" = "grey70")) +
  ggrepel::geom_text_repel(data = top, ggplot2::aes(label = gene), max.overlaps = 10, size = 3) +
  ggplot2::theme_minimal() + ggplot2::labs(x = "log2FC (HRPS vs HOPS)", y = "-log10(BH adjusted p-value)")
ggplot2::ggsave(file.path(out, "volcano.pdf"), p, width = 9, height = 7)
vst <- DESeq2::vst(dds, blind = FALSE)
selected <- head(subset(res_df, !is.na(padj))$gene, 20)
if (length(selected) >= 2) {
  mat <- SummarizedExperiment::assay(vst)[selected, , drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  grDevices::pdf(file.path(out, "top20_heatmap.pdf"), width = 9, height = 8)
  pheatmap::pheatmap(mat, annotation_col = metadata, show_colnames = TRUE)
  grDevices::dev.off()
}
writeLines(c("Contrast: HRPS vs HOPS; positive log2FC indicates HRPS higher.",
             "Sample-level pseudobulk raw counts; no cells treated as independent replicates.",
             "Design: ~ condition; does not correct for donor/batch covariates.",
             "Threshold: BH adjusted p < 0.05 and |log2FC| > 0.5.",
             paste("Independent samples: HRPS", ncol(pa), "; HOPS", ncol(pb)),
             paste("R:", R.version.string)), file.path(out, "analysis_notes.txt"))
message("Finished: ", normalizePath(out))
