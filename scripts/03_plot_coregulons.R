#!/usr/bin/env Rscript
# Coregulon heatmaps (Fig. 8A/B, Figs S11-S14) in the house expression-heatmap style.
#   Rscript scripts/03_plot_coregulons.R [SPECIES ...]

local({
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(file_arg)), "..", "R", "config.R"))
})
for (f in c("matrices.R", "coregulons.R", "heatmaps.R")) source(file.path(CFG$root, "R", f))

PANELS <- list(
  HYDRA        = "Fig8A_Hydra_coregulon",
  NEMATOSTELLA = "Fig8B_Nematostella_coregulon",
  MNEMIOPSIS   = "FigS11_Mnemiopsis_coregulon",
  SPONGILLA    = "FigS12_Spongilla_coregulon",
  MOUSE        = "FigS13_Mouse_coregulon",
  DROME        = "FigS14_Drosophila_coregulon"
)

short_id <- function(species, ids) switch(species, HYDRA = sub("aep$", "", ids), ids)

gene_labels <- function(species, members) {
  name <- members$name
  name[is.na(name) | name %in% c("", "-")] <- NA
  if (species == "MOUSE") {
    sym <- data.table::fread(mapping_path("mouse_ensembl2symbol.tab"), header = FALSE)
    name <- sym$V2[match(members$id, sym$V1)]
  }
  id <- short_id(species, members$id)
  ifelse(is.na(name) | name == id, id, paste(id, name))
}

seed_family <- function(species, ids) {
  fams <- as.data.frame(read_families())
  fam_order <- c(unlist(CFG$coregulon$seed_families), unlist(CFG$coregulon$seed_fallback))
  fams <- fams[fams$species == species & fams$family %in% fam_order, ]
  fams <- fams[order(match(fams$family, fam_order)), ]
  fams$family[match(ids, fams$id)]
}

FAMILY_COLORS <- setNames(brewer.pal(8, "Dark2"),
                          c("Acetyltransf", "Biopterin_H", "Pyridoxal-AADC", "DOMON_clade",
                            "P450-CYP2D6_clade", "SLC18", "SLC5", "SNF_clade"))
FAMILY_COLORS["SNF"] <- FAMILY_COLORS[["SNF_clade"]]

species <- commandArgs(trailingOnly = TRUE)
if (!length(species)) species <- names(PANELS)

for (sp in species) {
  members <- read_coregulon(sp)
  mat     <- load_matrix(species_config()[[sp]]$matrix)
  fam     <- seed_family(sp, members$id)
  fam[!members$is_seed] <- NA
  present <- intersect(names(FAMILY_COLORS), unique(fam))
  anno <- rowAnnotation(
    `Seed family` = factor(fam, levels = present),
    col = list(`Seed family` = FAMILY_COLORS[present]),
    na_col = "white", show_annotation_name = FALSE,
    simple_anno_size = unit(HM$cell_w_mm, "mm"),
    annotation_legend_param = list(`Seed family` = list(
      title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = HM$font),
      labels_gp = gpar(fontsize = 7, fontfamily = HM$font)))
  )
  ht <- expression_heatmap(mat, members$id, labels = gene_labels(sp, members),
                           cluster_rows = TRUE, left_annotation = anno)
  save_heatmap(ht, PANELS[[sp]])
}
