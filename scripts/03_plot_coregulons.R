#!/usr/bin/env Rscript
# Coregulon heatmaps (Fig. 8A-C, Figs S11-S13) in the house expression-heatmap style.
#   Rscript scripts/03_plot_coregulons.R [SPECIES ...]

local({
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(file_arg)), "..", "R", "config.R"))
})
for (f in c("matrices.R", "coregulons.R", "heatmaps.R")) source(file.path(CFG$root, "R", f))

PANELS <- list(
  HYDRA        = "Fig8A_Hydra_coregulon",
  NEMATOSTELLA = "Fig8B_Nematostella_coregulon",
  MNEMIOPSIS   = "Fig8C_Mnemiopsis_coregulon",
  SPONGILLA    = "FigS11_Spongilla_coregulon",
  MOUSE        = "FigS12_Mouse_coregulon",
  DROME        = "FigS13_Drosophila_coregulon"
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
  if (species %in% c("MOUSE", "DROME")) return(ifelse(is.na(name), id, name))   # bilaterians: gene symbols
  ifelse(is.na(name) | name == id, id, paste(id, name))
}

seed_family <- function(species, ids) {
  fams <- as.data.frame(read_families())
  fam_order <- c(unlist(CFG$coregulon$seed_families), unlist(CFG$coregulon$seed_fallback))
  fams <- fams[fams$species == species & fams$family %in% fam_order, ]
  fams <- fams[order(match(fams$family, fam_order)), ]
  fams$family[match(ids, fams$id)]
}

# Cell-type class used to group coregulon genes into blocks: the class of the cell type
# where the gene's scaled expression peaks. Classes come from the atlas annotation
# (prefix of the Levy et al. / Sebe-Pedros et al. names; Spongilla cell-type families of
# Musser et al. Table 1) or, for mouse and fly, from the lineage groupings below.
MOUSE_LINEAGE <- c(
  "Neural progenitor cells" = "neural", "Neural Tube" = "neural", "Postmitotic premature neurons" = "neural",
  "Sensory neurons" = "neural", "Granule neurons" = "neural", "Excitatory neurons" = "neural",
  "Inhibitory interneurons" = "neural", "Inhibitory neuron progenitors" = "neural", "Inhibitory neurons" = "neural",
  "Cholinergic neurons" = "neural", "Radial glia" = "glia", "Schwann cell precursor" = "glia",
  "Oligodendrocyte Progenitors" = "glia", "Premature oligodendrocyte" = "glia", "Ependymal cell" = "glia",
  "Isthmic organizer cells" = "neural", "Cardiac muscle lineages" = "muscle", "Myocytes" = "muscle",
  "Definitive erythroid lineage" = "blood", "Primitive erythroid lineage" = "blood", "Megakaryocytes" = "blood",
  "White blood cells" = "blood", "Endothelial cells" = "endothelium", "Epithelial cells" = "epithelium",
  "Lens" = "epithelium", "Hepatocytes" = "hepatocytes", "Melanocytes" = "melanocytes",
  "Chondroctye progenitors" = "mesenchyme", "Chondrocytes & osteoblasts" = "mesenchyme",
  "Connective tissue progenitors" = "mesenchyme", "Early mesenchyme" = "mesenchyme",
  "Intermediate Mesoderm" = "mesenchyme", "Jaw and tooth progenitors" = "mesenchyme",
  "Limb mesenchyme" = "mesenchyme", "Osteoblasts" = "mesenchyme", "Stromal cells" = "mesenchyme",
  "Notochord cells" = "notochord")
DROME_LINEAGE <- c(
  "epithelial cell" = "epithelia", "follicle cell" = "reproductive", "adult hindgut" = "epithelia",
  "polar follicle cell" = "reproductive", "enteroendocrine cell" = "epithelia", "eo support cell" = "epithelia",
  "adult fat body" = "fat body", "germline cell" = "reproductive", "female reproductive system" = "reproductive",
  "escort cell" = "reproductive", "follicle cell St. 9+" = "reproductive",
  "prefollicle cell/stalk follicle cell" = "reproductive", "spermatocyte" = "reproductive",
  "male accessory gland" = "reproductive", "cell body glial cell" = "glia", "adult glial cell" = "glia",
  "subperineurial glial cell" = "glia", "adult reticular neuropil associated glial cell" = "glia",
  "perineurial glial sheath" = "glia", "CNS surface associated glial cell" = "glia", "hemocyte" = "hemocytes",
  "indirect flight muscle" = "muscle", "muscle cell" = "muscle", "adult ventral nervous system" = "neurons",
  "adult peripheral nervous system" = "neurons", "multidendritic neuron" = "neurons",
  "leg muscle motor neuron" = "neurons", "adult oenocyte" = "oenocytes", "scolopidial neuron" = "neurons",
  "leg taste bristle chemosensory neuron" = "neurons", "gustatory receptor neuron" = "neurons",
  "adult tracheal cell" = "trachea", "unannotated" = "unannotated")

cell_type_class <- function(species, cell_types) {
  cls <- switch(species,
    HYDRA = , NEMATOSTELLA = sub("_.*", "", cell_types),
    MNEMIOPSIS = { x <- sub("_?cl[0-9]+$", "", cell_types); ifelse(x == "", "unannotated", gsub("_", " ", x)) },
    SPONGILLA = {
      t1 <- read.delim(require_file(sc_path("Musser2020/Table1_celltype_descriptions_final.tsv")), check.names = FALSE)
      unname(setNames(t1[["Cell Type Family"]], t1[["Cell Type Name"]])[cell_types])
    },
    MOUSE = unname(MOUSE_LINEAGE[cell_types]),
    DROME = unname(DROME_LINEAGE[cell_types]))
  if (anyNA(cls)) stop("cell types without class: ", paste(cell_types[is.na(cls)], collapse = ", "))
  cls
}

peak_blocks <- function(species, mat, ids) {
  cols <- cell_type_order(mat)
  z <- scale_rows(as.matrix(mat[ids, cols, drop = FALSE]))
  cls <- cell_type_class(species, cols)
  peak <- cls[apply(z, 1, which.max)]
  factor(peak, levels = unique(cls))
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
  blocks <- peak_blocks(sp, mat, members$id)
  ht <- expression_heatmap(mat, members$id, labels = gene_labels(sp, members),
                           row_split = blocks, cluster_rows = TRUE,
                           left_annotation = anno)
  save_heatmap(ht, PANELS[[sp]])
}
