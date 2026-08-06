#!/usr/bin/env Rscript
# Fallback exporter for EWAS Data Hub sample-info when only .RData is present.
# Prefer the Python path (mbs.registry.sample_info) which reads the .txt member.
#
# Usage:
#   Rscript scripts/export_ewas_sample_info.R FAMILY INPUT.zip OUTPUT.parquet
# Requires: R, arrow (or use write.csv then convert).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: export_ewas_sample_info.R FAMILY INPUT.zip OUTPUT.parquet")
}
family <- args[[1]]
zip_path <- args[[2]]
out_path <- args[[3]]

tmpdir <- tempfile(paste0("ewas_sample_", family, "_"))
dir.create(tmpdir, recursive = TRUE)
utils::unzip(zip_path, exdir = tmpdir)

rdata <- list.files(tmpdir, pattern = "\\.RData$", full.names = TRUE, ignore.case = TRUE)
txt <- list.files(tmpdir, pattern = "\\.txt$", full.names = TRUE, ignore.case = TRUE)

if (length(txt) >= 1) {
  message("Found .txt member; prefer Python export. Writing CSV sidecar for inspection.")
  tbl <- utils::read.table(txt[[1]], header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
} else if (length(rdata) >= 1) {
  env <- new.env(parent = emptyenv())
  load(rdata[[1]], envir = env)
  objs <- ls(envir = env)
  if (length(objs) < 1) stop("RData contained no objects")
  tbl <- get(objs[[1]], envir = env)
  if (!is.data.frame(tbl)) {
    tbl <- as.data.frame(tbl)
  }
} else {
  stop("No .txt or .RData member in zip")
}

csv_path <- sub("\\.parquet$", ".csv", out_path)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tbl, csv_path, row.names = FALSE)
message("Wrote ", csv_path)
message("Convert to parquet with: uv run python -c \"import pandas as pd; pd.read_csv(r'", csv_path, "').to_parquet(r'", out_path, "')\"")
