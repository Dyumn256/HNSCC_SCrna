library(Seurat)
library(ggplot2)
library(ggdist)
library(gghalves)
library(dplyr)

# ── CONFIG ────────────────────────────────────────────────────────────────────
data_dir <- "K:/Spatial"
out_dir  <- "K:/Spatial/plots_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_files <- c(
  "sample_1.Robj",  "sample_2.Robj",  "sample_3.Robj",
  "sample_4.Robj",  "sample_5.Robj",  "sample_6.Robj",
  "sample_7.Robj",  "sample_8.Robj",  "sample_9.Robj",
  "sample_10.Robj", "sample_11.Robj", "sample_12.Robj"
)

hrps_genes <- c("TMEM123", "LAPTM4B", "FAM162A", "FOS", "IGFBP7", 
                "LOXL2", "ASB2", "NDRG1", "STC2", "CA9", 
                "IGFBP3", "SPRR1B", "NMB", "CCL20", "FEZ1")

hops_genes <- c("KRT6A", "POLR2J3", "CSTB", "PERP", "TPM3", 
                "S100A8", "DSC2", "FTH1", "S100A6", "S100A16", 
                "CP", "CD24", "F3", "ALCAM", "DSG3")

# ── LOOP: collect metadata from all samples ───────────────────────────────────
all_meta <- list()

for (sfile in sample_files) {
  
  sample_id <- sub("\\.Robj$", "", sfile)
  robj_path <- file.path(data_dir, sfile)
  
  if (!file.exists(robj_path)) {
    warning("Not found, skipping: ", robj_path)
    next
  }
  
  message("Loading: ", sample_id)
  loaded_names <- load(robj_path)
  sobj <- get(loaded_names[1])
  
  # Score if not already present
  if (!"HRPS1" %in% colnames(sobj@meta.data)) {
    sobj <- AddModuleScore(sobj, features = list(hrps_genes), name = "HRPS", seed = 42)
  }
  if (!"HOPS1" %in% colnames(sobj@meta.data)) {
    sobj <- AddModuleScore(sobj, features = list(hops_genes), name = "HOPS", seed = 42)
  }
  
  # Check required columns exist
  has_annot <- "cluster_annotations" %in% colnames(sobj@meta.data)
  if (!has_annot) {
    warning(sample_id, ": 'cluster_annotations' missing — skipping")
    rm(sobj); gc(); next
  }
  
  # Pull relevant columns + tag sample
  tmp <- sobj@meta.data %>%
    select(
      HRPS_score = HRPS1,
      HOPS_score = HOPS1,
      niche      = cluster_annotations
    ) %>%
    mutate(sample = sample_id)
  
  all_meta[[sample_id]] <- tmp
  rm(sobj); gc()
}

# ── POOL & FILTER ─────────────────────────────────────────────────────────────
df_all <- bind_rows(all_meta) %>%
  filter(niche %in% c("core", "transitory", "edge")) %>%
  mutate(niche = factor(niche, levels = c("core", "transitory", "edge")))

n_patients <- length(unique(df_all$sample))
message("Pooled spots: ", nrow(df_all), " | Samples: ", n_patients)

# ── SIGNIFICANCE TESTING ──────────────────────────────────────────────────────
message("\n── HRPS score pairwise Wilcoxon ──")
print(pairwise.wilcox.test(df_all$HRPS_score, df_all$niche, p.adjust.method = "bonferroni"))

message("\n── HOPS score pairwise Wilcoxon ──")
print(pairwise.wilcox.test(df_all$HOPS_score, df_all$niche, p.adjust.method = "bonferroni"))

# ── PLOT FUNCTION ─────────────────────────────────────────────────────────────
region_cols <- c(
  "core"        = "#4DC6CE",
  "transitory" = "#E8D44D",
  "edge"        = "#E05A5A"
)

make_raincloud <- function(df, score_col, y_label, plot_title, out_file) {
  
  df$score <- df[[score_col]]
  y_max    <- max(df$score, na.rm = TRUE)
  y_min    <- min(df$score, na.rm = TRUE)
  y_range  <- y_max - y_min
  
  brk1_y  <- y_max + y_range * 0.08
  brk2_y  <- y_max + y_range * 0.18
  text1_y <- brk1_y + y_range * 0.03
  text2_y <- brk2_y + y_range * 0.03
  y_upper <- y_max + y_range * 0.30
  
  p <- ggplot(df, aes(x = niche, y = score, fill = niche, colour = niche)) +
    
    # Half violin (right side)
    geom_half_violin(
      side   = "r",
      alpha  = 0.85,
      trim   = FALSE,
      scale  = "width",
      colour = NA
    ) +
    
    # Jitter dots (left side)
    geom_jitter(
      width  = 0.06,
      height = 0,
      size   = 0.4,
      alpha  = 0.30
    ) +
    
    # Boxplot
    geom_boxplot(
      width              = 0.08,
      outlier.shape      = NA,
      fill               = "grey20",
      colour             = "grey20",
      alpha              = 0.9,
      median.linewidth   = 1.2
    ) +
    
    # Mean dot
    stat_summary(
      fun    = mean,
      geom   = "point",
      shape  = 21,
      size   = 2,
      fill   = "white",
      colour = "grey30",
      stroke = 0.8
    ) +
    
    scale_fill_manual(values   = region_cols) +
    scale_colour_manual(values = region_cols) +
    coord_cartesian(ylim = c(y_min - y_range * 0.05, y_upper)) +
    
    # Bracket: core vs transitory
    annotate("segment", x = 1,   xend = 2,   y = brk1_y, yend = brk1_y, colour = "grey20") +
    annotate("segment", x = 1,   xend = 1,   y = brk1_y - y_range*0.02, yend = brk1_y, colour = "grey20") +
    annotate("segment", x = 2,   xend = 2,   y = brk1_y - y_range*0.02, yend = brk1_y, colour = "grey20") +
    annotate("text",    x = 1.5, y = text1_y, label = "***", size = 4) +
    
    # Bracket: core vs edge
    annotate("segment", x = 1, xend = 3,   y = brk2_y, yend = brk2_y, colour = "grey20") +
    annotate("segment", x = 1, xend = 1,   y = brk2_y - y_range*0.02, yend = brk2_y, colour = "grey20") +
    annotate("segment", x = 3, xend = 3,   y = brk2_y - y_range*0.02, yend = brk2_y, colour = "grey20") +
    annotate("text",    x = 2, y = text2_y, label = "***", size = 4) +
    
    labs(title = plot_title, x = NULL, y = y_label) +
    
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, size = 11),
      axis.text.x     = element_text(size = 11),
      axis.text.y     = element_text(size = 10),
      axis.title.y    = element_text(size = 11, face = "bold"),
      axis.line       = element_line(colour = "grey40"),
      axis.ticks      = element_line(colour = "grey40")
    )
  
  ggsave(out_file, p, width = 4, height = 5.5, dpi = 300, bg = "white")
  message("Saved → ", out_file)
  return(p)
}

# ── GENERATE PLOTS ────────────────────────────────────────────────────────────
p_hrps <- make_raincloud(
  df         = df_all,
  score_col  = "HRPS_score",
  y_label    = "HRPS score",
  plot_title = paste0("n = ", n_patients, " patients"),
  out_file   = file.path(out_dir, "HRPS_raincloud_pooled.png")
)

p_hops <- make_raincloud(
  df         = df_all,
  score_col  = "HOPS_score",
  y_label    = "HOPS score",
  plot_title = paste0("n = ", n_patients, " patients"),
  out_file   = file.path(out_dir, "HOPS_raincloud_pooled.png")
)

message("\nDone! Plots saved to: ", out_dir)