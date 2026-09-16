# Configuration and path resolution shared by all R scripts.
#
# Precedence (highest first): NEUROTRANSEVO_<KEY> environment variables,
# config.local.yml, config.yml. Relative paths resolve against the repo root.

suppressPackageStartupMessages(library(yaml))

find_repo_root <- function() {
  env <- Sys.getenv("NEUROTRANSEVO_ROOT")
  if (nzchar(env)) return(normalizePath(env))
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  starts <- c(if (length(file_arg)) dirname(normalizePath(file_arg[1])), getwd())
  for (dir in starts) {
    repeat {
      if (file.exists(file.path(dir, "config.yml")) &&
          file.exists(file.path(dir, "environment.yml"))) return(dir)
      parent <- dirname(dir)
      if (identical(parent, dir)) break
      dir <- parent
    }
  }
  stop("Cannot locate the repository root (config.yml). Run from inside the ",
       "repository or set NEUROTRANSEVO_ROOT.")
}

merge_config <- function(base, override) {
  for (k in names(override)) {
    base[[k]] <- if (is.list(base[[k]]) && is.list(override[[k]]) && !is.null(names(override[[k]])))
      merge_config(base[[k]], override[[k]]) else override[[k]]
  }
  base
}

load_config <- function(root = find_repo_root()) {
  cfg <- read_yaml(file.path(root, "config.yml"))
  local <- file.path(root, "config.local.yml")
  if (file.exists(local)) cfg <- merge_config(cfg, read_yaml(local))
  for (k in names(cfg$paths)) {
    env <- Sys.getenv(paste0("NEUROTRANSEVO_", toupper(k)))
    if (nzchar(env)) cfg$paths[[k]] <- env
  }
  cfg$paths <- lapply(cfg$paths, function(p)
    if (grepl("^(/|~|[A-Za-z]:)", p)) path.expand(p) else file.path(root, p))
  cfg$root <- root
  cfg
}

CFG <- load_config()

# Path helpers ----------------------------------------------------------------
sc_path      <- function(...) file.path(CFG$paths$singlecell, ...)
eggnog_path  <- function(...) file.path(CFG$paths$eggnog, ...)
mapping_path <- function(...) file.path(CFG$paths$mappings, ...)
results_path <- function(..., create = TRUE) {
  p <- file.path(CFG$paths$results, ...)
  if (create) dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  p
}
cache_path <- function(name) results_path("matrices", paste0(name, ".rds"))

require_file <- function(path, hint = "see data/README.md") {
  if (!file.exists(path)) stop(sprintf("Missing input: %s (%s)", path, hint), call. = FALSE)
  path
}

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
