# Amine-metabolism coregulons.
#
# For each species: Pearson correlation (across cell types) between every seed
# gene (members of the seed families) and every expressed gene; BH correction
# per seed within the candidate pool (GO:0009308 genes plus seeds); a gene
# joins when r >= r_min and BH-adjusted p <= padj_max for at least one active
# seed (a seed with >= 1 non-self partner). The coregulon is the union of the
# GO-annotated partners and the active seeds.

suppressPackageStartupMessages(library(data.table))

EGGNOG_COLS <- c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
                 "max_annot_lvl", "COG_category", "Description", "Preferred_name",
                 "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
                 "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")

# ---- Species configuration: matrix, eggNOG file, eggNOG query -> matrix ID ------
species_config <- function() list(
  HYDRA = list(
    matrix = "hydra",
    eggnog = "Hydra_vulgaris.6087.hydra.emapper.annotations",
    to_matrix_id = function(x) sub("^6087\\.(t[0-9]+aep)_.*", "\\1", sub("^6087\\.(t[0-9]+aep)$", "\\1", x))
  ),
  NEMATOSTELLA = list(
    matrix = "nematostella",
    eggnog = "Nematostella_vectensis.45351.nematostella.emapper.annotations",
    to_matrix_id = function(x) sub("^45351\\.", "", x)
  ),
  MNEMIOPSIS = list(
    matrix = "mnemiopsis",
    eggnog = "Mnemiopsis_leidyi.27923.MneLei_Aug2011.emapper.annotations",
    to_matrix_id = function(x) sub("-PA$", "", sub("^27923\\.", "", x))
  ),
  SPONGILLA = list(
    matrix = "spongilla",
    eggnog = "Spongilla_lacustris.6055.spongilla.emapper.annotations",
    # eggNOG queries are proteins (6055.m.<n>); the matrix uses Trinity genes (c<n>_g<n>).
    to_matrix_id = local({
      g2p <- fread(require_file(mapping_path("spongilla_gene2protein.tab")), header = FALSE,
                   col.names = c("gene", "protein"))
      p2g <- setNames(g2p$gene, g2p$protein)
      function(x) unname(p2g[x])
    })
  ),
  MOUSE = list(
    matrix = "mouse",
    eggnog = "Mus_musculus.10090.emapper.emapper.annotations",
    # eggNOG queries are UniProt accessions; the matrix uses Ensembl gene IDs.
    to_matrix_id = local({
      u2e <- fread(require_file(mapping_path("mouse_uniprot2ensembl.tab")), header = FALSE,
                   col.names = c("uniprot", "ensembl"))
      m <- setNames(u2e$ensembl, u2e$uniprot)
      function(x) unname(m[sub("^10090\\.", "", x)])
    })
  ),
  DROME = list(
    matrix = "drosophila",
    eggnog = "Drosophila_melanogaster.7227.emapper.annotations",
    # eggNOG queries are UniProt accessions; the matrix uses FlyBase gene symbols,
    # taken from the GN= tag of the UniProt reference proteome headers.
    to_matrix_id = local({
      heads <- grep("^>", readLines(require_file(file.path(CFG$paths$proteomes, "UP000000803_7227.fasta"))),
                    value = TRUE)
      acc <- sub("^>(sp|tr)\\|([^|]+)\\|.*", "\\2", heads)
      gn  <- ifelse(grepl(" GN=", heads), sub(".* GN=([^ ]+) .*", "\\1", heads), NA_character_)
      a2g <- setNames(sub("^Dmel\\\\", "", gn), acc)
      function(x) unname(a2g[sub("^7227\\.", "", x)])
    })
  )
)

# ---- Inputs -------------------------------------------------------------------------
read_families <- function() {
  fread(require_file(CFG$paths$families), header = FALSE, sep = "\t",
        col.names = c("family", "module", "name", "id", "species"))
}

seed_genes <- function(species_code, families = read_families()) {
  cc   <- CFG$coregulon
  fams <- families[families$species == species_code]
  seeds <- fams$id[fams$family %in% cc$seed_families]
  for (fam in names(cc$seed_fallback)) {
    if (!any(fams$family == fam)) seeds <- c(seeds, fams$id[fams$family == cc$seed_fallback[[fam]]])
  }
  unique(seeds)
}

read_eggnog <- function(file) {
  fread(require_file(eggnog_path(file)), sep = "\t", header = FALSE, skip = "#query\t",
        fill = TRUE, blank.lines.skip = TRUE, quote = "", col.names = EGGNOG_COLS)
}

# ---- Computation ----------------------------------------------------------------------
seed_correlations <- function(mat, seeds) {
  mat   <- mat[rowSums(mat) > 0, , drop = FALSE]
  seeds <- intersect(seeds, rownames(mat))
  x     <- t(as.matrix(mat))
  n     <- nrow(x)
  r <- matrix(NA_real_, nrow = nrow(mat), ncol = length(seeds), dimnames = list(rownames(mat), seeds))
  p <- r
  for (s in seeds) {
    xs <- x[, s]
    if (sd(xs, na.rm = TRUE) == 0) next
    rs <- apply(x, 2, function(y) {
      if (sd(y, na.rm = TRUE) == 0) return(NA_real_)
      suppressWarnings(cor(xs, y, method = "pearson", use = "pairwise.complete.obs"))
    })
    tstat  <- rs * sqrt(n - 2) / sqrt(1 - rs^2)
    r[, s] <- rs
    p[, s] <- 2 * pt(-abs(tstat), df = n - 2)
  }
  list(r = r, p = p, n_cell_types = n)
}

compute_coregulon <- function(species) {
  cfg <- species_config()[[species]]
  cc  <- CFG$coregulon
  mat <- load_matrix(cfg$matrix)
  ann <- read_eggnog(cfg$eggnog)
  ann[, matrix_id := cfg$to_matrix_id(query)]
  go_ids <- unique(ann[grepl(cc$go_term, GOs, fixed = TRUE) & !is.na(matrix_id) & matrix_id != "", matrix_id])
  seeds  <- seed_genes(species)

  rp   <- seed_correlations(mat, seeds)
  pool <- union(intersect(go_ids, rownames(rp$r)), intersect(colnames(rp$r), rownames(rp$r)))
  r    <- rp$r[pool, , drop = FALSE]
  padj <- apply(rp$p[pool, , drop = FALSE], 2, function(col) {
    out <- rep(NA_real_, length(col)); ok <- !is.na(col)
    if (any(ok)) out[ok] <- p.adjust(col[ok], method = "BH")
    out
  })
  rownames(padj) <- pool

  pass <- (r >= cc$r_min) & (padj <= cc$padj_max)
  pass[is.na(pass)] <- FALSE
  active <- Filter(function(s) sum(pass[setdiff(rownames(pass), s), s]) >= 1,
                   intersect(colnames(pass), rownames(pass)))
  pass <- pass[, active, drop = FALSE]
  partners <- rownames(r)[rowSums(pass) >= 1]
  members  <- union(intersect(partners, go_ids), active)

  rr <- r[members, active, drop = FALSE]; pp <- padj[members, active, drop = FALSE]
  best <- apply(rr, 1, function(v) if (all(is.na(v))) NA_integer_ else which.max(v))
  info <- ann[match(members, ann$matrix_id)]
  table <- data.table(
    id          = members,
    is_seed     = members %in% active,
    max_r       = apply(rr, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)),
    min_padj    = apply(pp, 1, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)),
    best_seed   = colnames(rr)[best],
    name        = info$Preferred_name,
    description = info$Description
  )
  list(species = species, members = table, seeds = seeds, active_seeds = active,
       go_universe = length(go_ids), n_cell_types = rp$n_cell_types,
       r = rp$r, p = rp$p, candidate_pool = pool)
}

write_coregulon <- function(res) {
  out <- results_path("coregulons", paste0(res$species, ".members.tsv"))
  fwrite(res$members, out, sep = "\t")
  saveRDS(res[c("r", "p", "n_cell_types", "candidate_pool", "active_seeds")],
          results_path("coregulons", paste0(res$species, ".correlations.rds")))
  log_msg("%-13s %3d cell types | %3d seeds (%d active) | GO universe %5d | coregulon %3d genes",
          res$species, res$n_cell_types, length(res$seeds), length(res$active_seeds),
          res$go_universe, nrow(res$members))
  invisible(out)
}

read_coregulon <- function(species) {
  fread(require_file(results_path("coregulons", paste0(species, ".members.tsv"), create = FALSE),
                     "run scripts/02_compute_coregulons.R first"), sep = "\t")
}
