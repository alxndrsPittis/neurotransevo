#!/usr/bin/env Rscript
# Compute amine-metabolism coregulons (thresholds in config.yml).
#   Rscript scripts/02_compute_coregulons.R                  # all species
#   Rscript scripts/02_compute_coregulons.R HYDRA NEMATOSTELLA

local({
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(file_arg)), "..", "R", "config.R"))
})
source(file.path(CFG$root, "R", "matrices.R"))
source(file.path(CFG$root, "R", "coregulons.R"))

species <- commandArgs(trailingOnly = TRUE)
if (!length(species)) species <- names(species_config())

for (sp in species) write_coregulon(compute_coregulon(sp))
