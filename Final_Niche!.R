# ══════════════════════════════════════════════════════════════════════════════
#  Spatial H&E Overlay — Master Loop
#  Produces 3 plots per sample:
#    1. HRPS score
#    2. HOPS score
#    3. Niche annotation (categorical)
# ══════════════════════════════════════════════════════════════════════════════

library(Seurat)
library(png)
library(scales)

# ── USER CONFIG ───────────────────────────────────────────────────────────────

data_dir   <- "K:/Spatial"                     # folder with sample .Robj files
raw_dir    <- "K:/Spatial/GSE208253_RAW"       # folder with .gz images & matrices
out_dir    <- "K:/Spatial/plots_output"        # where PNGs will be saved

# All sample Robj filenames (add/remove as needed)
sample_files <- c(
  "sample_1.Robj",
  "sample_2.Robj",
  "sample_3.Robj",
  "sample_4.Robj",
  "sample_5.Robj",
  "sample_6.Robj",
  "sample_7.Robj",
  "sample_8.Robj",
  "sample_9.Robj",
  "sample_10.Robj",
  "sample_11.Robj",
  "sample_12.Robj"
)

# Gene sets for module scoring
hrps_genes <- c("TMEM123", "LAPTM4B", "FAM162A", "FOS", "IGFBP7", 
                "LOXL2", "ASB2", "NDRG1", "STC2", "CA9", 
                "IGFBP3", "SPRR1B", "NMB", "CCL20", "FEZ1")

hops_genes <- c("KRT6A", "POLR2J3", "CSTB", "PERP", "TPM3", 
                "S100A8", "DSC2", "FTH1", "S100A6", "S100A16", 
                "CP", "CD24", "F3", "ALCAM", "DSG3")

# Colour palettes
magma_pal_cols  <- c("#0D0887","#6A00A8","#CB4679","#F89441","#F0F921")  # HRPS
ylgnbu_pal_cols <- c("#FFFFD9","#EDF8B1","#7FCDBB","#1D91C0","#0C2C84") # HOPS
region_cols     <- c(
  "core"        = "#E69F00",
  "transitory" = "#56B4E9",
  "edge"        = "#009E73",
  "nc"          = "grey70"
)

# Global plot settings
spot_r_px <- 4.5
out_w      <- 8
res        <- 300

# ── HELPER FUNCTIONS ──────────────────────────────────────────────────────────

#' Draw H&E + continuous score overlay, save PNG
plot_continuous <- function(he, img_h, img_w,
                            px_col, px_row,
                            vals, palette_cols,
                            plot_title, legend_label,
                            out_file,
                            spot_r_px = 4.5, alpha = 0.80,
                            out_w = 8, res = 300) {
  
  q      <- quantile(vals, c(0.02, 0.98), na.rm = TRUE)
  clp    <- pmax(pmin(vals, q[2]), q[1])
  nrm    <- (clp - q[1]) / (q[2] - q[1])
  pal    <- colorRampPalette(palette_cols)(256)
  cols   <- pal[pmax(1L, ceiling(nrm * 255) + 1L)]
  
  asp    <- img_h / img_w
  out_h  <- out_w * asp + 0.8
  
  png(out_file, width = out_w, height = out_h, units = "in", res = res)
  par(mar = c(0.5, 0.5, 2.0, 4.5), bg = "white", xpd = FALSE)
  
  plot.new()
  plot.window(xlim = c(0, img_w), ylim = c(img_h, 0))
  rasterImage(he, 0, img_h, img_w, 0)
  
  symbols(x = px_col, y = px_row,
          circles = rep(spot_r_px, length(px_col)),
          inches = FALSE, add = TRUE, fg = NA,
          bg = adjustcolor(cols, alpha.f = alpha))
  
  title(main = plot_title, cex.main = 1.2, font.main = 2, line = 0.6)
  
  # Colour bar
  n_leg    <- 200
  usr      <- par("usr")
  bar_x0   <- usr[2] + (usr[2] - usr[1]) * 0.02
  bar_x1   <- bar_x0 + (usr[2] - usr[1]) * 0.03
  bar_y0   <- img_h * 0.15
  bar_y1   <- img_h * 0.85
  leg_y    <- seq(bar_y0, bar_y1, length.out = n_leg + 1)
  leg_cols <- colorRampPalette(palette_cols)(n_leg)
  
  for (i in seq_len(n_leg))
    rect(bar_x0, leg_y[i], bar_x1, leg_y[i + 1],
         col = leg_cols[n_leg - i + 1L], border = NA, xpd = TRUE)
  rect(bar_x0, bar_y0, bar_x1, bar_y1,
       col = NA, border = "grey40", lwd = 0.7, xpd = TRUE)
  
  brks <- pretty(c(q[1], q[2]), n = 5)
  brks <- brks[brks >= q[1] & brks <= q[2]]
  ty   <- bar_y0 + (1 - (brks - q[1]) / (q[2] - q[1])) * (bar_y1 - bar_y0)
  text(bar_x1 + (usr[2] - usr[1]) * 0.01, ty,
       labels = formatC(brks, digits = 2, format = "f"),
       cex = 0.7, adj = 0, xpd = TRUE)
  text((bar_x0 + bar_x1) / 2, bar_y1 + img_h * 0.05,
       labels = legend_label, cex = 0.75, font = 2, xpd = TRUE, adj = 0.5)
  
  dev.off()
  message("  Saved → ", out_file)
}

#' Draw H&E + categorical niche overlay, save PNG
plot_categorical <- function(he, img_h, img_w,
                             px_col, px_row,
                             clusters, palette_cols,
                             plot_title, legend_title,
                             out_file,
                             spot_r_px = 4.5, alpha = 0.88,
                             out_w = 8, res = 300) {
  
  spot_cols <- palette_cols[as.character(clusters)]
  
  asp   <- img_h / img_w
  out_h <- out_w * asp + 0.8
  
  png(out_file, width = out_w, height = out_h, units = "in", res = res)
  par(mar = c(0.5, 0.5, 2.0, 4.5), bg = "white", xpd = FALSE)
  
  plot.new()
  plot.window(xlim = c(0, img_w), ylim = c(img_h, 0))
  rasterImage(he, 0, img_h, img_w, 0)
  
  symbols(x = px_col, y = px_row,
          circles = rep(spot_r_px, length(px_col)),
          inches = FALSE, add = TRUE, fg = NA,
          bg = adjustcolor(spot_cols, alpha.f = alpha))
  
  title(main = plot_title, cex.main = 1.2, font.main = 2, line = 0.6)
  
  usr         <- par("usr")
  cats        <- names(palette_cols)
  dot_r       <- (usr[2] - usr[1]) * 0.015
  spacing     <- img_h * 0.08
  leg_x       <- usr[2] + (usr[2] - usr[1]) * 0.02
  leg_y_start <- img_h * 0.25
  
  text(leg_x + dot_r, leg_y_start - spacing * 0.8,
       labels = legend_title, cex = 0.8, font = 2,
       adj = c(0.3, 0.5), xpd = TRUE)
  
  for (i in seq_along(cats)) {
    cy <- leg_y_start + (i - 1) * spacing
    symbols(leg_x + dot_r, cy, circles = dot_r,
            inches = FALSE, add = TRUE, fg = NA,
            bg = adjustcolor(palette_cols[cats[i]], alpha.f = alpha),
            xpd = TRUE)
    text(leg_x + dot_r * 2.8, cy, labels = cats[i],
         cex = 0.75, adj = c(0, 0.5), xpd = TRUE)
  }
  
  dev.off()
  message("  Saved → ", out_file)
}

# ── MAIN LOOP ─────────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (sfile in sample_files) {
  
  sample_id <- sub("\\.Robj$", "", sfile)    # e.g. "sample_1"
  message("\n══ Processing: ", sample_id, " ══")
  
  # ── Load object ────────────────────────────────────────────────────────────
  robj_path <- file.path(data_dir, sfile)
  if (!file.exists(robj_path)) {
    warning("File not found, skipping: ", robj_path)
    next
  }
  
  loaded_names <- load(robj_path)
  sobj <- get(loaded_names[1])
  
  # ── Detect image name (first available image slot) ─────────────────────────
  img_name <- names(sobj@images)[1]
  message("  Image slot: ", img_name)
  
  # ── Module scoring (AddModuleScore adds columns HRPS1, HOPS1, etc.) ──────────
  message("  Scoring gene modules...")
  
  if (!"HRPS1" %in% colnames(sobj@meta.data)) {
    sobj <- AddModuleScore(sobj,
                           features = list(hrps_genes),
                           name     = "HRPS",
                           seed     = 42)
  }
  
  if (!"HOPS1" %in% colnames(sobj@meta.data)) {
    sobj <- AddModuleScore(sobj,
                           features = list(hops_genes),
                           name     = "HOPS",
                           seed     = 42)
  }
  
  # ── Extract shared coords ──────────────────────────────────────────────────
  img_info <- sobj@images[[img_name]]
  he       <- img_info@image
  img_h    <- dim(he)[1]
  img_w    <- dim(he)[2]
  coords   <- GetTissueCoordinates(sobj, image = img_name)
  
  # ── Per-sample output subfolder ────────────────────────────────────────────
  samp_out <- file.path(out_dir, sample_id)
  dir.create(samp_out, showWarnings = FALSE)
  
  # ══ PLOT 1: HRPS score ══════════════════════════════════════════════════════
  message("  Plot 1/3 — HRPS score")
  hrps_vals <- sobj$HRPS1
  shared    <- intersect(rownames(coords), names(hrps_vals))
  
  plot_continuous(
    he           = he,
    img_h        = img_h, img_w = img_w,
    px_col       = coords[shared, "imagecol"],
    px_row       = coords[shared, "imagerow"],
    vals         = hrps_vals[shared],
    palette_cols = magma_pal_cols,
    plot_title   = "HRPS score",
    legend_label = "HRPS\nscore",
    out_file     = file.path(samp_out, paste0(sample_id, "_HRPS_score_HE.png")),
    spot_r_px    = spot_r_px, out_w = out_w, res = res
  )
  
  # ══ PLOT 2: HOPS score ═════════════════════════════════════════════════════
  message("  Plot 2/3 — HOPS score")
  hops_vals <- sobj$HOPS1
  shared    <- intersect(rownames(coords), names(hops_vals))
  
  plot_continuous(
    he           = he,
    img_h        = img_h, img_w = img_w,
    px_col       = coords[shared, "imagecol"],
    px_row       = coords[shared, "imagerow"],
    vals         = hops_vals[shared],
    palette_cols = ylgnbu_pal_cols,
    plot_title   = "HOPS score",
    legend_label = "HOPS\nscore",
    out_file     = file.path(samp_out, paste0(sample_id, "_HOPS_score_HE.png")),
    spot_r_px    = spot_r_px, out_w = out_w, res = res
  )
  
  # ══ PLOT 3: Niche annotation ════════════════════════════════════════════════
  message("  Plot 3/3 — Niche annotation")
  
  if ("cluster_annotations" %in% colnames(sobj@meta.data)) {
    clusters <- sobj$cluster_annotations
    shared   <- intersect(rownames(coords), names(clusters))
    
    plot_categorical(
      he           = he,
      img_h        = img_h, img_w = img_w,
      px_col       = coords[shared, "imagecol"],
      px_row       = coords[shared, "imagerow"],
      clusters     = clusters[shared],
      palette_cols = region_cols,
      plot_title   = "Tumor Niche",
      legend_title = "Region",
      out_file     = file.path(samp_out, paste0(sample_id, "_niche_HE.png")),
      spot_r_px    = spot_r_px, out_w = out_w, res = res
    )
  } else {
    message("  !! 'cluster_annotations' not found — skipping niche plot")
  }
  
  # Clean up to free memory before next sample
  rm(sobj)
  gc()
}

message("\n✓ All done! Plots saved to: ", out_dir)