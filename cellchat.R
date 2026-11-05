#!/usr/bin/env Rscript
# cellchat.R — ligand–receptor communication inference for one Seurat dataset.
# Usage: Rscript cellchat.R seurat_object.rds output_dir human|mouse [identity_column]
# If identity_column is omitted, the object's active Seurat identities are used.
# Requires: Seurat, CellChat, ggplot2. The chosen organism must match gene IDs.
# CellChat estimates putative interactions from expression, not experimental proof.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:4) stop("Usage: Rscript cellchat.R object.rds output_dir human|mouse [identity_column]")
for (pkg in c("Seurat", "CellChat", "ggplot2")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(CellChat))
set.seed(42)
obj <- readRDS(args[1]); out <- args[2]; species <- tolower(args[3]); dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!inherits(obj, "Seurat")) stop("RDS must contain a Seurat object")
if (!species %in% c("human", "mouse")) stop("Species must be human or mouse")
if (!"RNA" %in% Assays(obj)) stop("Seurat object needs RNA assay")
DefaultAssay(obj) <- "RNA"
if ("JoinLayers" %in% getNamespaceExports("Seurat")) obj <- JoinLayers(obj, assay = "RNA")
if (length(args) == 4L) {
  if (!args[4] %in% names(obj[[]])) stop("Identity column not found: ", args[4])
  groups <- as.character(obj[[]][, args[4]])
} else groups <- as.character(Idents(obj))
if (length(groups) != ncol(obj) || anyNA(groups) || any(!nzchar(groups))) stop("Missing/invalid cell group labels")
meta <- data.frame(group = groups, row.names = colnames(obj))
if (any(table(meta$group) < 10L)) warning("Some groups contain <10 cells and may be filtered")
# A single species-specific database; do NOT overwrite human entries with mouse.
db <- if (species == "human") CellChatDB.human else CellChatDB.mouse
db <- subsetDB(db, search = "Secreted Signaling")
expr <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
if (length(expr@x) == 0L) stop("RNA normalized data are empty. Run NormalizeData before CellChat.")
chat <- createCellChat(object = expr, meta = meta, group.by = "group")
chat@DB <- db
chat <- subsetData(chat)
chat <- identifyOverExpressedGenes(chat)
chat <- identifyOverExpressedInteractions(chat)
# Set raw.use=TRUE to avoid requiring PPI projection (optional and species-specific).
chat <- computeCommunProb(chat, raw.use = TRUE)
chat <- filterCommunication(chat, min.cells = 10)
chat <- computeCommunProbPathway(chat)
chat <- aggregateNet(chat)
saveRDS(chat, file.path(out, "cellchat_results.rds"))
write.csv(as.data.frame(table(group = meta$group)), file.path(out, "cell_group_sizes.csv"), row.names = FALSE)
write.csv(chat@net$count, file.path(out, "interaction_count_matrix.csv"))
write.csv(chat@net$weight, file.path(out, "interaction_weight_matrix.csv"))
communications <- subsetCommunication(chat)
write.csv(communications, file.path(out, "ligand_receptor_interactions.csv"), row.names = FALSE)
writeLines(as.character(chat@netP$pathways), file.path(out, "detected_pathways.txt"))
group_sizes <- as.numeric(table(chat@idents))
grDevices::pdf(file.path(out, "overview_network.pdf"), width = 13, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(chat@net$count, vertex.weight = group_sizes, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction count")
netVisual_circle(chat@net$weight, vertex.weight = group_sizes, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
grDevices::dev.off()
if (length(chat@netP$pathways)) {
  # Only plot a pathway that is actually present (no hard-coded TNF/CCL names).
  pathway <- chat@netP$pathways[1]
  grDevices::pdf(file.path(out, "example_pathway.pdf"), width = 8, height = 8)
  netVisual_aggregate(chat, signaling = pathway, layout = "circle")
  grDevices::dev.off()
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP")
  saveRDS(chat, file.path(out, "cellchat_results.rds"))
}
writeLines(c(paste("Species:", species), "Database subset: Secreted Signaling",
             "CellChat normalized expression input from RNA assay; raw.use=TRUE.",
             "Inspect sample composition/batch differences before biological interpretation.",
             paste("R:", R.version.string)), file.path(out, "analysis_notes.txt"))
message("Finished: ", normalizePath(out))
