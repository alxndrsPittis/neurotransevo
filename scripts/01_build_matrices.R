#!/usr/bin/env Rscript
# Build cell-type-averaged expression matrices.
#   Rscript scripts/01_build_matrices.R                 # all datasets
#   Rscript scripts/01_build_matrices.R hydra mouse     # selected datasets
#   Rscript scripts/01_build_matrices.R --rebuild       # ignore existing results

local({
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(normalizePath(file_arg)), "..", "R", "config.R"))
})
source(file.path(CFG$root, "R", "matrices.R"))

args    <- commandArgs(trailingOnly = TRUE)
rebuild <- "--rebuild" %in% args
which   <- setdiff(args, "--rebuild")
if (!length(which)) which <- names(MATRICES)
unknown <- setdiff(which, names(MATRICES))
if (length(unknown)) stop("Unknown dataset(s): ", paste(unknown, collapse = ", "))

for (name in which) invisible(load_matrix(name, rebuild = rebuild))
