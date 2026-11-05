#!/usr/bin/env Rscript
# cibersortx1.R — join CIBERSORTx immune fractions with sample metadata.
# Usage: Rscript cibersortx1.R CIBERSORTx_results.csv clinical.csv output_dir
# Clinical CSV must include SampleID and Clonal_Phenotype. CIBERSORTx must include
# Mixture and numeric immune-fraction columns. Two groups are compared:
# 'High HOPS / Low HRPS' and 'High HRPS / Low HOPS'.
# Output: sample matching audit, immune statistics with BH FDR, and PDF plots.
# Requires: ggplot2, dplyr, tidyr. No packages are installed automatically.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript cibersortx1.R fractions.csv clinical.csv output_dir")
for (p in c("ggplot2", "dplyr", "tidyr")) if (!requireNamespace(p, quietly = TRUE)) stop("Install R package: ", p)
out <- args[3]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
fractions <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
clin <- read.csv(args[2], check.names = FALSE, stringsAsFactors = FALSE)
if (!"Mixture" %in% names(fractions)) stop("CIBERSORTx results lack Mixture")
if (!all(c("SampleID", "Clonal_Phenotype") %in% names(clin))) stop("Clinical CSV requires SampleID and Clonal_Phenotype")
if (anyDuplicated(fractions$Mixture) || anyDuplicated(clin$SampleID)) stop("Sample IDs must be unique within each input file")
# Prefer canonical LM22 cell names; avoid assuming columns 2:23 are immune cells.
lm22 <- c("B cells naive", "B cells memory", "Plasma cells", "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated", "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta", "NK cells resting", "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated", "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils")
cells <- intersect(lm22, names(fractions))
if (length(cells) != 22L) stop("Expected all 22 LM22 columns; found ", length(cells), ". Check your signature matrix and headers.")
for (cell in cells) {
  fractions[[cell]] <- suppressWarnings(as.numeric(fractions[[cell]]))
  if (any(!is.finite(fractions[[cell]]))) stop("Non-numeric/missing immune fraction in ", cell)
}
write.csv(data.frame(Mixture = fractions$Mixture, matched = fractions$Mixture %in% clin$SampleID), file.path(out, "sample_matching_audit.csv"), row.names = FALSE)
merged <- dplyr::inner_join(fractions, clin, by = c("Mixture" = "SampleID"))
valid <- c("High HOPS / Low HRPS", "High HRPS / Low HOPS")ndat <- merged[merged$Clonal_Phenotype %in% valid, , drop = FALSE]
dat$Clonal_Phenotype <- factor(dat$Clonal_Phenotype, levels = valid, labels = c("HOPS (Specialist)", "HRPS (Architect)"))
counts <- table(dat$Clonal_Phenotype)
if (any(counts < 2L)) stop("At least 2 matched samples per condition required; got: ", paste(counts, collapse = ", "))
write.csv(data.frame(group = names(counts), n = as.integer(counts)), file.path(out, "group_sizes.csv"), row.names = FALSE)
stats <- do.call(rbind, lapply(cells, function(cell) {
  x <- dat[[cell]][dat$Clonal_Phenotype == levels(dat$Clonal_Phenotype)[1]]
  y <- dat[[cell]][dat$Clonal_Phenotype == levels(dat$Clonal_Phenotype)[2]]
  p <- tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.frame(CellType = cell, mean_HOPS = mean(x), mean_HRPS = mean(y),
             median_HOPS = median(x), median_HRPS = median(y), p_value = p)
}))
stats$FDR_BH <- p.adjust(stats$p_value, method = "BH")
stats <- stats[order(stats$FDR_BH), ]
write.csv(stats, file.path(out, "immune_statistics.csv"), row.names = FALSE)
long <- tidyr::pivot_longer(dat[, c("Mixture", "Clonal_Phenotype", cells)], cols = dplyr::all_of(cells), names_to = "CellType", values_to = "Fraction")
p <- ggplot2::ggplot(long, ggplot2::aes(Clonal_Phenotype, Fraction, fill = Clonal_Phenotype)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  ggplot2::geom_jitter(width = 0.1, alpha = 0.45, size = 0.7) +
  ggplot2::facet_wrap(~CellType, scales = "free_y", ncol = 4) +
  ggplot2::scale_fill_manual(values = c("HOPS (Specialist)" = "#E64B35", "HRPS (Architect)" = "#4DBBD5")) +
  ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none") +
  ggplot2::labs(x = NULL, y = "Estimated immune-cell fraction", title = "CIBERSORTx LM22 immune fractions (all 22 cell types)")
ggplot2::ggsave(file.path(out, "immune_boxplots.pdf"), p, width = 15, height = 12)
summary <- dplyr::summarise(dplyr::group_by(long, CellType, Clonal_Phenotype), Mean = mean(Fraction), SEM = stats::sd(Fraction)/sqrt(dplyr::n()), .groups = "drop")
p2 <- ggplot2::ggplot(summary, ggplot2::aes(CellType, Mean, fill = Clonal_Phenotype)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, Mean - SEM), ymax = Mean + SEM), position = ggplot2::position_dodge(0.8), width = 0.2) +
  ggplot2::coord_flip() + ggplot2::theme_minimal() + ggplot2::labs(x = NULL, y = "Mean estimated fraction", fill = "Phenotype")
ggplot2::ggsave(file.path(out, "immune_means.pdf"), p2, width = 11, height = 9)
writeLines("Exploratory analysis: Wilcoxon sample-level tests, BH correction across 22 cell types. CIBERSORTx fractions are estimates and compositional; interpret with caution. Ensure each row is an independent biological sample.", file.path(out, "analysis_notes.txt"))
message("Finished: ", normalizePath(out))
