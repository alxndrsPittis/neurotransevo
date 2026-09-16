# House style for all single-cell expression heatmaps (Figs 4D/E, 5C/D, 7, 8, S4-S14).
#
# - columns: cell types ordered by complete-linkage clustering of 1 - Pearson r on
#   expression of the whole transcriptome (identical for every panel of a
#   species), split into non-neuronal | neuronal blocks; no dendrograms
# - rows: z-scored (scaled) expression, rev(RdYlBu) palette, legend "Scaled expr"
# - panel annotations (family labels, brackets, silhouettes) are added at figure
#   assembly and are not drawn here

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(grid)
})

HM <- list(
  palette   = rev(brewer.pal(7, "RdYlBu")),
  font      = "Helvetica",
  fontsize  = 7,
  cell_w_mm = 3.2,
  cell_h_mm = 3.2,
  gap_mm    = 1
)

NEURON_REGEX <- "neuron|nervous|neural|i_n_|neurosecretory|N\\(|N.early.state"

# Column order shared by all panels of a species -------------------------------
# Complete linkage on 1 - Pearson r of cell-type-averaged expression, reversed
# dendrogram; neuronal cell types are then moved to a right-hand block keeping
# their relative order. Reproduces the column order of the published figures.
cell_type_order <- function(mat) {
  hc <- hclust(as.dist(1 - cor(as.matrix(mat), use = "pairwise.complete.obs", method = "pearson")),
               method = "complete")
  o <- labels(rev(as.dendrogram(hc)))
  c(o[neuron_split(o) == "O"], o[neuron_split(o) == "N"])
}

neuron_split <- function(cell_types) {
  factor(ifelse(grepl(NEURON_REGEX, cell_types, ignore.case = TRUE), "N", "O"), levels = c("O", "N"))
}

scale_rows <- function(m) {
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  z
}

# Generic house-style heatmap ------------------------------------------------------
# mat:          full genes x cell types matrix of the species
# genes/labels: row IDs to show (in display order) and their labels
# row_split:    NULL, a vector of group names (kept in order of appearance), or an
#               integer k to cut the row dendrogram into k blocks
# cluster_rows: order rows by complete-linkage clustering of 1 - Pearson r
# keep_empty:   keep genes without expression as grey rows (default: drop them)
expression_heatmap <- function(mat, genes, labels = genes, row_split = NULL,
                               cluster_rows = FALSE, keep_empty = FALSE,
                               row_titles = TRUE, show_row_names = TRUE,
                               left_annotation = NULL, right_annotation = NULL) {
  cols <- cell_type_order(mat)
  present <- genes %in% rownames(mat)
  m <- matrix(0, nrow = length(genes), ncol = length(cols), dimnames = list(genes, cols))
  m[present, ] <- as.matrix(mat[genes[present], cols, drop = FALSE])
  empty <- rowSums(m) == 0
  keep <- if (keep_empty) rep(TRUE, length(genes)) else !empty
  m <- m[keep, , drop = FALSE]; empty <- empty[keep]
  rownames(m) <- labels[keep]
  if (length(row_split) == length(genes)) {
    lv <- if (is.factor(row_split)) levels(droplevels(row_split[keep])) else unique(row_split[keep])
    row_split <- factor(as.character(row_split[keep]), levels = lv)
  }
  if (!is.null(left_annotation)) left_annotation <- left_annotation[keep, ]
  z <- scale_rows(m)
  z[empty, ] <- NA

  pearson_dist <- function(x) {
    d <- 1 - cor(t(x), use = "pairwise.complete.obs")
    d[!is.finite(d)] <- 1
    as.dist(d)
  }
  row_clust <- FALSE
  if (isTRUE(cluster_rows) && sum(!empty) > 1) {
    # With a grouping vector, rows are clustered within each group (same distance/linkage).
    row_clust <- if (is.factor(row_split)) TRUE else as.dendrogram(hclust(pearson_dist(m), method = "complete"))
  }

  gp <- function(size = HM$fontsize, face = "plain") gpar(fontsize = size, fontfamily = HM$font, fontface = face)
  Heatmap(
    z,
    name = "Scaled expr",
    col = HM$palette,
    na_col = "grey75",
    cluster_columns = FALSE,
    cluster_rows = row_clust,
    clustering_distance_rows = pearson_dist,
    clustering_method_rows = "complete",
    row_dend_reorder = FALSE,
    show_row_dend = FALSE,
    column_split = neuron_split(colnames(z)),
    column_title = NULL,
    column_gap = unit(HM$gap_mm, "mm"),
    row_split = row_split,
    cluster_row_slices = FALSE,
    row_title = if (row_titles) character(0) else NULL,
    show_row_names = show_row_names,
    row_title_rot = 0,
    row_title_gp = gp(8),
    row_gap = unit(HM$gap_mm, "mm"),
    border = FALSE,
    rect_gp = gpar(col = NA),
    row_names_gp = gp(),
    column_names_gp = gp(),
    column_names_rot = 45,
    column_names_side = "bottom",
    width = unit(ncol(z) * HM$cell_w_mm, "mm"),
    height = unit(nrow(z) * HM$cell_h_mm, "mm"),
    left_annotation = left_annotation,
    right_annotation = right_annotation,
    heatmap_legend_param = list(title = "Scaled expr", direction = "horizontal",
                                title_gp = gp(7, "bold"), labels_gp = gp(7),
                                legend_width = unit(2.5, "cm"))
  )
}

# Pearson correlation of panel genes (rows) with partner genes (columns) across
# cell types; non-significant values (p > 0.05) are set to 0, as in Fig. 5C/D.
correlation_heatmap <- function(mat, genes, partners, partner_labels = partners, row_labels = genes,
                                p_max = 0.05) {
  suppressPackageStartupMessages(library(Hmisc))
  partners <- intersect(partners, rownames(mat))
  genes_in <- genes[genes %in% rownames(mat)]
  rc <- rcorr(t(as.matrix(mat[unique(c(genes_in, partners)), , drop = FALSE])))
  r <- matrix(0, nrow = length(genes), ncol = length(partners), dimnames = list(genes, partners))
  sig <- rc$r[genes_in, partners, drop = FALSE]
  sig[rc$P[genes_in, partners, drop = FALSE] > p_max | is.na(rc$P[genes_in, partners, drop = FALSE])] <- 0
  r[genes_in, ] <- sig
  colnames(r) <- partner_labels
  rownames(r) <- row_labels
  gp <- function(size = HM$fontsize, face = "plain") gpar(fontsize = size, fontfamily = HM$font, fontface = face)
  Heatmap(r, name = "Pearson's r",
          col = circlize::colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
          cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gp(),
          column_names_side = "top", column_names_rot = 45, column_names_gp = gp(),
          rect_gp = gpar(col = NA),
          width = unit(ncol(r) * HM$cell_w_mm, "mm"),
          heatmap_legend_param = list(title = "Pearson's r", direction = "horizontal",
                                      title_gp = gp(7, "bold"), labels_gp = gp(7),
                                      legend_width = unit(2.5, "cm")))
}

save_heatmap <- function(ht, stem) {
  dir <- dirname(results_path("figures", paste0(stem, ".pdf")))
  size <- local({
    pdf(NULL); on.exit(dev.off())
    d <- draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE)
    w <- ComplexHeatmap:::width(d); h <- ComplexHeatmap:::height(d)
    c(convertWidth(w, "inches", valueOnly = TRUE), convertHeight(h, "inches", valueOnly = TRUE))
  })
  size <- size + c(0.4, 0.8)
  for (dev in c("pdf", "svg")) {
    path <- file.path(dir, paste0(stem, ".", dev))
    # pdf(): text stays text (standard Helvetica); svglite: real <text> elements.
    # Both are editable in Inkscape / Illustrator.
    if (dev == "pdf") pdf(path, width = size[1], height = size[2], family = HM$font, useDingbats = FALSE)
    else svglite::svglite(path, width = size[1], height = size[2],
                          system_fonts = list(sans = HM$font))
    draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE)
    dev.off()
  }
  log_msg("figure: %s (%.1f x %.1f in)", file.path(dir, stem), size[1], size[2])
}
