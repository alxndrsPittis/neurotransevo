#!/usr/bin/env Rscript
# Gene-family expression heatmaps (Figs 4D/E, 5C/D, 7A/B, S4-S10) in the house style.
# Row sets and orders follow the curated family table (data/families/) and, where the
# published figure used a hand-curated phylogeny order, data/families/row_orders/.
#   Rscript scripts/04_plot_family_heatmaps.R [PANEL ...]

local({
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(file_arg)), "..", "R", "config.R"))
})
for (f in c("matrices.R", "coregulons.R", "heatmaps.R")) source(file.path(CFG$root, "R", f))

FAMS <- as.data.frame(read_families())

ENZYMES_TRANSPORTERS <- c("Acetyltransf", "Biopterin_H", "Pyridoxal-AADC", "Pyridoxal-GAD", "DOMON_clade",
                          "P450-CYP2D6_clade", "PNMT", "SLC18", "VIAAT", "SLC5", "SNF_clade")
CHANNELS <- c("Neur_chan", "Lig_chan", "DLG_clade", "Rapsyn_clade")
TRANSPORTERS_SLC22 <- c("SLC22_clade", "SV2", "SVOP")

SPECIES_MATRIX <- c(MOUSE = "mouse", DROME = "drosophila", HYDRA = "hydra", NEMATOSTELLA = "nematostella")

read_order <- function(file) {
  x <- readLines(file.path(CFG$root, "data/families/row_orders", file))
  x <- x[!grepl("^#", x)]
  block <- cumsum(x == "") + 1
  list(ids = x[x != ""], block = block[x != ""])
}

# Rows of a panel: members of `families` (table order), grouped by family.
family_rows <- function(species, families) {
  f <- FAMS[FAMS$species == species & FAMS$family %in% families, ]
  f <- f[order(match(f$family, families)), ]
  f[!duplicated(f$id), c("id", "name", "family")]
}

PANELS <- list(
  Fig4D_mouse_mAChR_clade   = list(species = "MOUSE", families = "7tm_1_subcladeACh", keep_empty = TRUE, split = FALSE),
  Fig4E_hydra_mAChR_clade   = list(species = "HYDRA", families = "7tm_1_subcladeACh", keep_empty = TRUE, split = FALSE),
  Fig5C_mouse_nAChR         = list(species = "MOUSE", order = "Fig5C_mouse_nachr.txt", keep_empty = TRUE, correlate = "Rapsyn_clade"),
  Fig5D_hydra_nAChR         = list(species = "HYDRA", order = "Fig5D_hydra_nachr.txt", keep_empty = TRUE, correlate = "Rapsyn_clade"),
  Fig7A_mouse_enzymes_transporters = list(species = "MOUSE", families = setdiff(ENZYMES_TRANSPORTERS, "P450-CYP2D6_clade")),
  Fig7B_hydra_enzymes_transporters = list(species = "HYDRA", families = ENZYMES_TRANSPORTERS),
  FigS4_mouse_channels      = list(species = "MOUSE", families = c("Neur_chan_ACh", "Neur_chan_GABA", "Lig_chan", "DLG_clade", "Rapsyn_clade"),
                                   order = "FigS4_mouse_cysloop.txt", keep_empty = TRUE,
                                   groups = c(Neur_chan_ACh = "Cys-loop", Neur_chan_GABA = "Cys-loop")),
  FigS5_drosophila_channels = list(species = "DROME", families = CHANNELS, keep_empty = TRUE),
  FigS6_hydra_channels      = list(species = "HYDRA", families = CHANNELS, keep_empty = TRUE),
  FigS7_nematostella_channels = list(species = "NEMATOSTELLA", families = CHANNELS, keep_empty = TRUE),
  FigS8A_drosophila_enzymes_transporters = list(species = "DROME", families = setdiff(ENZYMES_TRANSPORTERS, c("P450-CYP2D6_clade", "PNMT"))),
  FigS8B_nematostella_enzymes_transporters = list(species = "NEMATOSTELLA", families = setdiff(ENZYMES_TRANSPORTERS, "PNMT")),
  FigS9A_mouse_SLC22        = list(species = "MOUSE", families = TRANSPORTERS_SLC22),
  FigS9B_drosophila_SLC22   = list(species = "DROME", families = TRANSPORTERS_SLC22),
  FigS10A_hydra_SLC22       = list(species = "HYDRA", families = TRANSPORTERS_SLC22),
  FigS10B_nematostella_SLC22 = list(species = "NEMATOSTELLA", families = TRANSPORTERS_SLC22)
)

plot_panel <- function(name, spec) {
  mat <- load_matrix(SPECIES_MATRIX[[spec$species]])
  if (!is.null(spec$families)) {
    rows <- family_rows(spec$species, spec$families)
    group <- rows$family
    if (!is.null(spec$groups)) group[rows$family %in% names(spec$groups)] <- spec$groups[rows$family[rows$family %in% names(spec$groups)]]
    if (!is.null(spec$order)) {           # curated order inside the first group
      o <- read_order(spec$order)$ids
      first <- group == group[1]
      key <- ifelse(rows$id[first] %in% o, rows$id[first], rows$name[first])
      idx <- which(first)[order(match(key, o))]
      rows <- rbind(rows[idx, ], rows[!first, ]); group <- c(group[idx], group[!first])
    }
    split <- if (isFALSE(spec$split)) NULL else group
  } else {
    o <- read_order(spec$order)
    sp_fams <- FAMS[FAMS$species == spec$species, ]
    hit <- ifelse(o$ids %in% sp_fams$id, match(o$ids, sp_fams$id), match(o$ids, sp_fams$name))
    rows <- sp_fams[hit, c("id", "name", "family")]
    split <- o$block
    blocks_untitled <- TRUE
  }
  blocks_untitled <- is.null(spec$families)
  keep_empty <- isTRUE(spec$keep_empty)
  if (!keep_empty) {
    expressed <- rows$id %in% rownames(mat)[rowSums(mat) > 0]
    rows <- rows[expressed, ]; if (length(split) > 1) split <- split[expressed]
  }
  ht <- expression_heatmap(mat, rows$id, labels = rows$name, row_split = split,
                           keep_empty = keep_empty, row_titles = !blocks_untitled,
                           show_row_names = is.null(spec$correlate))
  if (!is.null(spec$correlate)) {
    partners <- family_rows(spec$species, spec$correlate)
    ht <- ht + correlation_heatmap(mat, rows$id, partners$id, partners$name, row_labels = rows$name)
  }
  save_heatmap(ht, name)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

which <- commandArgs(trailingOnly = TRUE)
if (!length(which)) which <- names(PANELS)
for (p in which) plot_panel(p, PANELS[[p]])
