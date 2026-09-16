# Cell-type-averaged expression matrices (genes x cell types).
#
# For every dataset, UMI counts are summed per annotated cell type and divided
# by the number of cells of that type. Each builder reads the public data from
# CFG$paths$singlecell (layout documented in data/README.md).

suppressPackageStartupMessages(library(Matrix))

# ---- Aggregation ---------------------------------------------------------------
# Locale-independent (C-locale) factor, so column order is identical on every OS.
c_factor <- function(x, levels = sort(unique(x[!is.na(x)]), method = "radix")) factor(x, levels = levels)

aggregate_cols <- function(mat, by) {
  groups <- if (is.factor(by)) by else c_factor(by)
  out <- sapply(levels(groups), function(lvl) Matrix::rowSums(mat[, groups == lvl, drop = FALSE]))
  rownames(out) <- rownames(mat)
  out
}

normalise_by_clusters <- function(mat_sum, clusters) {
  sweep(mat_sum, MARGIN = 2, STATS = table(clusters), FUN = "/")
}

cell_type_means <- function(mat, cell_types) {
  cls <- c_factor(cell_types)
  normalise_by_clusters(aggregate_cols(mat, cls), cls)
}

# ---- Builders --------------------------------------------------------------------

# Hydra vulgaris: Siebert et al. 2019 UMI matrix, cell types of Levy et al. 2021.
build_hydra <- function() {
  dat <- read.delim(require_file(sc_path("Siebert2019/wt_expression_matrix.txt")), row.names = 1)
  colnames(dat) <- gsub("\\.", "-", gsub("^X", "", colnames(dat)))
  cells <- read.delim(require_file(sc_path("Arnau2021/Hydra/Hvul_cell_type_assignments")), row.names = 1)
  cells <- cells[cells$cell_type != "", , drop = FALSE]
  common <- intersect(rownames(cells), colnames(dat))
  m <- cell_type_means(as.matrix(dat[, common]), cells[common, "cell_type"])
  rownames(m) <- vapply(strsplit(rownames(m), "|", fixed = TRUE), `[`, character(1), 1L)
  m
}

# Nematostella vectensis: adult MARS-seq atlas and cell types of Levy et al. 2021.
build_nematostella <- function() {
  dat <- readRDS(require_file(sc_path("Arnau2021/Nematostella/Nvec_adult_sc_UMI_counts.RDS")))
  dat <- dat[!grepl("mOrange", rownames(dat)), ]
  cells <- read.delim(require_file(sc_path("Arnau2021/Nematostella/Nvec_adult_cell_type_assignments")),
                      row.names = 1)
  cells <- cells[cells$cell_type != 0 & !is.na(cells$cell_type), , drop = FALSE]
  m <- cell_type_means(as.matrix(dat[, rownames(cells)]), cells$cell_type)
  rownames(m) <- gsub("^Nvec_", "", rownames(m))
  m
}

# Mnemiopsis leidyi: Sebe-Pedros et al. 2018 (GSM3021563) metacells and annotation.
build_mnemiopsis <- function() {
  dat   <- read.delim(require_file(sc_path("Arnau2018/GSM3021563_Mnemiopsis_leidyi_UMI_table.txt")),
                      row.names = 1)
  cells <- read.delim(require_file(sc_path("Arnau2018/GSM3021563_Mnemiopsis_leidyi_metacell_definition.txt")),
                      row.names = 1, header = FALSE, colClasses = c("character", "character"))
  anno  <- read.delim(require_file(sc_path("Arnau2018/GSM3021563_Mnemiopsis_leidyi_metacell_assignments.txt")),
                      row.names = 1, header = FALSE)
  cell_type <- anno[match(cells[[1]], rownames(anno)), 1]
  cell_type_means(as.matrix(dat[, rownames(cells)]), cell_type)
}

# Spongilla lacustris: Musser et al. 2021 (GSE134912), cell-type names from Table 1.
build_spongilla <- function() {
  dat  <- read.delim(require_file(sc_path("Musser2020/GSE134912_spongilla_10x_count_matrix.txt")),
                     row.names = 1, check.names = FALSE)
  meta <- read.table(require_file(sc_path("Musser2020/clusters_42_34_and_oddcells_reassigned.txt")),
                     row.names = 1, sep = " ")
  dat <- dat[, rownames(meta)]
  rownames(dat) <- vapply(strsplit(rownames(dat), " ", fixed = TRUE), `[`, character(1), 1L)
  m <- cell_type_means(as.matrix(dat), meta[[1]])
  anno <- read.delim(require_file(sc_path("Musser2020/Table1_celltype_descriptions_final.tsv")),
                     check.names = FALSE)
  name_by_cl <- setNames(anno[["Cell Type Name"]], as.character(anno[["Cluster #"]]))
  m <- m[, colnames(m) %in% names(name_by_cl), drop = FALSE]
  colnames(m) <- unname(name_by_cl[colnames(m)])
  m
}

# Mus musculus: Cao et al. 2019 organogenesis atlas, 37 main cell types.
MOUSE_CELL_TYPES <- c(
  "Neural progenitor cells", "Neural Tube", "Postmitotic premature neurons",
  "Sensory neurons", "Granule neurons", "Excitatory neurons",
  "Inhibitory interneurons", "Inhibitory neuron progenitors", "Inhibitory neurons",
  "Cholinergic neurons", "Radial glia", "Schwann cell precursor",
  "Cardiac muscle lineages", "Lens", "Chondroctye progenitors",
  "Chondrocytes & osteoblasts", "Connective tissue progenitors",
  "Definitive erythroid lineage", "Early mesenchyme", "Endothelial cells",
  "Ependymal cell", "Epithelial cells", "Hepatocytes", "Intermediate Mesoderm",
  "Isthmic organizer cells", "Jaw and tooth progenitors", "Limb mesenchyme",
  "Megakaryocytes", "Melanocytes", "Myocytes", "Notochord cells",
  "Oligodendrocyte Progenitors", "Osteoblasts", "Premature oligodendrocyte",
  "Primitive erythroid lineage", "Stromal cells", "White blood cells")
MOUSE_NEURON_TYPES <- MOUSE_CELL_TYPES[1:10]

build_mouse <- function() {
  dat   <- readRDS(require_file(sc_path("Cao2019/gene_count_cleaned.RDS")))
  cells <- read.csv(require_file(sc_path("Cao2019/cell_annotate.csv")), row.names = 1)[
             colnames(dat), "Main_cell_type"]
  keep  <- !is.na(cells) & cells %in% MOUSE_CELL_TYPES
  cls   <- factor(cells[keep], levels = MOUSE_CELL_TYPES)
  m     <- normalise_by_clusters(aggregate_cols(dat[, keep, drop = FALSE], cls), cls)
  rownames(m) <- vapply(strsplit(rownames(m), ".", fixed = TRUE), `[`, character(1), 1L)
  m
}

# Drosophila melanogaster: Fly Cell Atlas body 10x (Li et al. 2022), fine annotation.
# Read directly from the loom (HDF5) file in chunks of cells.
build_drosophila <- function(chunk = 10000L) {
  suppressPackageStartupMessages(library(hdf5r))
  f <- H5File$new(require_file(sc_path("Li2021/s_fca_biohub_body_10x.loom")), mode = "r")
  on.exit(f$close_all(), add = TRUE)
  genes <- f[["row_attrs/Gene"]]$read()
  annot <- f[["col_attrs/annotation"]]$read()
  keep  <- which(annot != "artefact" & !is.na(annot))
  cls   <- c_factor(annot[keep])
  mat   <- f[["matrix"]]            # cells x genes as seen from R
  sums  <- matrix(0, nrow = length(genes), ncol = nlevels(cls),
                  dimnames = list(genes, levels(cls)))
  for (start in seq(1L, length(keep), by = chunk)) {
    idx <- keep[start:min(start + chunk - 1L, length(keep))]
    rs  <- rowsum(mat[idx, , drop = FALSE], group = as.character(cls[match(idx, keep)]))
    sums[, rownames(rs)] <- sums[, rownames(rs)] + t(rs)
  }
  normalise_by_clusters(sums, cls)
}

# ---- Registry --------------------------------------------------------------------
MATRICES <- list(
  hydra        = list(build = build_hydra,        species = "HYDRA"),
  nematostella = list(build = build_nematostella, species = "NEMATOSTELLA"),
  mnemiopsis   = list(build = build_mnemiopsis,   species = "MNEMIOPSIS"),
  spongilla    = list(build = build_spongilla,    species = "SPONGILLA"),
  mouse        = list(build = build_mouse,        species = "MOUSE"),
  drosophila   = list(build = build_drosophila,   species = "DROME")
)

load_matrix <- function(name, rebuild = FALSE) {
  path <- cache_path(name)
  if (!rebuild && file.exists(path)) return(readRDS(path))
  log_msg("building matrix: %s", name)
  m <- MATRICES[[name]]$build()
  saveRDS(m, path)
  log_msg("  %s: %d genes x %d cell types -> %s", name, nrow(m), ncol(m), path)
  m
}
