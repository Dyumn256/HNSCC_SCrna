#!/usr/bin/env Rscript
# CombinedUMAP.R — compare two Seurat objects on a joint UMAP.
# Usage: Rscript CombinedUMAP.R HRPS.rds HOPS.rds output/umap
# Inputs: two Seurat objects containing raw RNA counts. Cells must represent
# comparable assays; merged UMAP does NOT correct batch effects.
# Output: combined Seurat RDS, UMAP PDF, cluster compositions CSV/PDF, markers CSV.
# Requires: Seurat (v4/v5), ggplot2, patchwork. Optional: MAST for marker test.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript CombinedUMAP.R HRPS.rds HOPS.rds output_dir")
packages <- c("Seurat", "ggplot2", "patchwork")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(Seurat))
set.seed(42)
hrps <- readRDS(args[1]); hops <- readRDS(args[2]); out <- args[3]
if (!inherits(hrps, "Seurat") || !inherits(hops, "Seurat")) stop("Both RDS files must contain Seurat objects")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (!"RNA" %in% Assays(hrps) || !"RNA" %in% Assays(hops)) stop("Both objects need an RNA assay")
DefaultAssay(hrps) <- "RNA"; DefaultAssay(hops) <- "RNA"
hrps$condition <- "HRPS"; hops$condition <- "HOPS"
merged <- merge(hrps, y = hops, add.cell.ids = c("HRPS", "HOPS"), project = "Spheroid_Study")
# Seurat v5 can retain distinct count layers after merging. Join them before
# normalization and downstream analyses (Seurat v4 has no JoinLayers).
if ("JoinLayers" %in% getNamespaceExports("Seurat")) merged <- JoinLayers(merged, assay = "RNA")
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
npcs <- min(20L, ncol(Embeddings(merged, "pca")))
if (npcs < 2L) stop("Insufficient PCA dimensions; inspect data size and quality")
merged <- FindNeighbors(merged, dims = seq_len(npcs), verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)
merged <- RunUMAP(merged, dims = seq_len(npcs), seed.use = 42, verbose = FALSE)
saveRDS(merged, file.path(out, "combined_seurat.rds"))

p1 <- DimPlot(merged, reduction = "umap", group.by = "condition", cols = c(HRPS = "#E41A1C", HOPS = "#377EB8")) + ggplot2::ggtitle("Merged UMAP by phenotype")
p2 <- DimPlot(merged, reduction = "umap", label = TRUE) + ggplot2::ggtitle("Unsupervised clusters")
ggplot2::ggsave(file.path(out, "umap.pdf"), p1 + p2, width = 13, height = 5)
meta <- merged[[]]
tab <- as.data.frame(table(cluster = as.character(meta$seurat_clusters), condition = meta$condition))
names(tab)[3] <- "cells"
tab <- tab[tab$cells > 0, ]
tab$percent_within_cluster <- 100 * tab$cells / ave(tab$cells, tab$cluster, FUN = sum)
tab$percent_within_condition <- 100 * tab$cells / ave(tab$cells, tab$condition, FUN = sum)
write.csv(tab, file.path(out, "cluster_composition.csv"), row.names = FALSE)
p3 <- ggplot2::ggplot(tab, ggplot2::aes(cluster, percent_within_cluster, fill = condition)) +
  ggplot2::geom_col() + ggplot2::labs(x = "Cluster", y = "% of cells within cluster") + ggplot2::theme_classic()
p4 <- ggplot2::ggplot(tab, ggplot2::aes(cluster, percent_within_condition, fill = condition)) +
  ggplot2::geom_col(position = "dodge") + ggplot2::labs(x = "Cluster", y = "% of each condition") + ggplot2::theme_classic()
ggplot2::ggsave(file.path(out, "cluster_composition.pdf"), p3 / p4, width = 10, height = 10)

# MAST is optional; use Wilcoxon if it is absent. These are cell-level marker
# screens and must NOT be interpreted as biological-replicate differential tests.
test <- if (requireNamespace("MAST", quietly = TRUE)) "MAST" else "wilcox"
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = test, verbose = FALSE)
write.csv(markers, file.path(out, "cluster_markers.csv"), row.names = FALSE)
genes <- intersect(c("MKI67", "TOP2A", "CD274", "TGFB1", "HLA-A", "AREG", "EREG", "MIF", "VIM", "SNAI2", "COL1A1"), rownames(merged))
if (length(genes) > 0) {
  p <- DotPlot(merged, features = genes, group.by = "seurat_clusters") + RotatedAxis()
  ggplot2::ggsave(file.path(out, "selected_gene_dotplot.pdf"), p, width = 11, height = 5)
}
writeLines(c(paste("R:", R.version.string), paste("Marker test:", test),
             "Merged visualization does not include batch integration.",
             "Cells are not biological replicates; validate differential results across samples."),
           file.path(out, "analysis_notes.txt"))
message("Finished: ", normalizePath(out))
