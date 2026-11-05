#!/usr/bin/env Rscript
# drugsensitivity3.R — evaluate an externally defined HOPS gene signature in GEO counts.
# Historical name retained for compatibility; this script scores expression,
# it does NOT calculate drug IC50 or infer drug sensitivity.
# Usage: Rscript drugsensitivity3.R raw_counts_dir sample_metadata.csv HOPS_genes.txt output_dir
# Input counts: one .txt/.tsv optionally .gz per sample, with 2 tab-separated
# columns (Ensembl gene ID, integer read count), no header (HTSeq format).
# Filenames must begin with a GEO GSM ID (e.g., GSM1234567_counts.txt.gz).
# Metadata CSV: ID,Group columns; ID must match GSM IDs in filenames.
# Gene list: one human HGNC gene symbol per line; defined INDEPENDENTLY of test data.
# Requires: ggplot2, org.Hs.eg.db, AnnotationDbi, ggdist (optional).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript drugsensitivity3.R raw_dir metadata.csv genes.txt output_dir")
for (p in c("ggplot2", "org.Hs.eg.db", "AnnotationDbi")) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)
rawdir <- args[1]; out <- args[4]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- list.files(rawdir, pattern = "\\.(txt|tsv)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No .txt/.tsv[.gz] files in ", rawdir)
ids <- sub("^((GSM[0-9]+)).*$", "\\1", basename(files))
if (any(!grepl("^GSM[0-9]+$", ids)) || anyDuplicated(ids)) stop("Files must start with unique GSM IDs")
meta <- read.csv(args[2], stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("ID", "Group") %in% names(meta))) stop("Metadata needs ID and Group columns")
if (anyDuplicated(meta$ID)) stop("Metadata ID must be unique")
genes <- unique(trimws(readLines(args[3], warn = FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes) < 2L) stop("Provide at least 2 pre-specified HGNC genes")
# Parse HTSeq counts; all values must be nonnegative integer raw counts.
read_count <- function(file) {
  df <- read.delim(file, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
  if (ncol(df) != 2L) stop("Expected exactly 2 columns in ", basename(file))
  df <- df[!startsWith(df[[1]], "__"), , drop = FALSE]
  df[[1]] <- sub("\\.[0-9]+$", "", df[[1]])
  counts <- suppressWarnings(as.numeric(df[[2]]))
  if (anyNA(counts) || any(counts < 0) || any(counts != floor(counts))) stop("Invalid raw counts in ", basename(file))
  # Collapse duplicate Ensembl IDs by SUM, not by averaging read counts.
  tapply(counts, df[[1]], sum)
}
counts_list <- lapply(files, read_count); names(counts_list) <- ids
common <- Reduce(intersect, lapply(counts_list, names))
if (!length(common)) stop("No shared Ensembl IDs across input samples")
counts <- vapply(counts_list, function(x) as.numeric(x[common]), numeric(length(common)))
rownames(counts) <- common; colnames(counts) <- ids
# Keep only metadata-matched files; fail if metadata samples are missing.
missing_ids <- setdiff(meta$ID, ids)
if (length(missing_ids)) stop("Metadata samples missing count files: ", paste(missing_ids, collapse = ", "))
meta <- meta[match(ids, meta$ID), , drop = FALSE]
if (anyNA(meta$ID) || anyNA(meta$Group)) stop("Every count file requires metadata")
# Preserve Ensembl IDs until mapping. Many-to-one gene symbols: sum counts.
map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(counts),
                             columns = "SYMBOL", keytype = "ENSEMBL")
map <- map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), c("ENSEMBL", "SYMBOL")]
# Exclude ambiguous Ensembl-to-multiple-symbol mappings.
map <- map[!duplicated(map$ENSEMBL) & !duplicated(map$ENSEMBL, fromLast = TRUE), , drop = FALSE]
map <- map[map$ENSEMBL %in% rownames(counts), , drop = FALSE]
if (!nrow(map)) stop("No human gene symbols mapped; confirm organism and Ensembl IDs")
by_symbol <- rowsum(counts[map$ENSEMBL, , drop = FALSE], group = map$SYMBOL, reorder = FALSE)
# Log2(CPM+1) scoring: simple descriptive signature, not a validated
# pharmacogenomic drug-response predictor or DE test.
lib_sizes <- colSums(counts)
if (any(lib_sizes <= 0)) stop("Some libraries have zero total counts")
logcpm <- log2(sweep(by_symbol, 2, lib_sizes / 1e6, "/") + 1)
available <- intersect(genes, rownames(logcpm))
writeLines(setdiff(genes, available), file.path(out, "unmatched_signature_genes.txt"))
if (length(available) < 2L) stop("Fewer than 2 signature genes matched; see unmatched_signature_genes.txt")
scores <- colMeans(logcpm[available, , drop = FALSE])
result <- data.frame(ID = names(scores), Group = meta$Group, HOPS_score = unname(scores), stringsAsFactors = FALSE)
write.csv(result, file.path(out, "signature_scores.csv"), row.names = FALSE)
writeLines(available, file.path(out, "matched_signature_genes.txt"))
p <- ggplot2::ggplot(result, ggplot2::aes(Group, HOPS_score, fill = Group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  ggplot2::geom_jitter(width = 0.1, height = 0, size = 1.8, alpha = 0.8) +
  ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none") +
  ggplot2::labs(title = "External HOPS signature expression", subtitle = paste(length(available), "matched genes"), x = "Experimental group", y = "Mean log2(CPM+1)")
ggplot2::ggsave(file.path(out, "signature_scores.pdf"), p, width = 10, height = 6)
# Pairwise contrasts are NOT hard-coded: users may examine the group scores and
# define a planned statistical comparison that respects experimental design.
writeLines(c("Exploratory external expression validation only; no drug sensitivity prediction.",
             "Input: per-sample raw integer HTSeq counts; human gene-symbol mapping.",
             "Signature score: unweighted mean of log2(CPM+1) for matched HGNC symbols.",
             "No automatic t-tests: specify design, replication and planned contrasts separately.",
             paste("Signature genes requested:", length(genes)),
             paste("Signature genes matched:", length(available)),
             paste("R:", R.version.string)), file.path(out, "analysis_notes.txt"))
message("Finished: ", normalizePath(out))
