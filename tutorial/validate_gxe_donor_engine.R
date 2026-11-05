#!/usr/bin/env Rscript

## Compare Step 2 output from the donor-aggregate and cell-reference engines.
##
## Usage:
##   Rscript validate_gxe_donor_engine.R CELL_OUTPUT DONOR_OUTPUT

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: validate_gxe_donor_engine.R CELL_OUTPUT DONOR_OUTPUT")
}

read_step2 <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

split_numeric <- function(x) {
  do.call(
    rbind,
    lapply(strsplit(as.character(x), ",", fixed = TRUE), as.numeric)
  )
}

cell <- read_step2(args[[1]])
donor <- read_step2(args[[2]])

if (!identical(cell$MarkerID, donor$MarkerID)) {
  stop("MarkerID columns are not identical")
}

fields <- c("Beta_ge", "seBeta_ge", "pval_ge", "pval_noSPA_ge")
tolerance <- 5e-6

for (field in fields) {
  cell_value <- split_numeric(cell[[field]])
  donor_value <- split_numeric(donor[[field]])
  difference <- abs(cell_value - donor_value)
  maximum <- max(difference, na.rm = TRUE)
  cat(sprintf("%-16s max absolute difference: %.3e\n", field, maximum))
  if (!is.finite(maximum) || maximum > tolerance) {
    stop(field, " differs by more than ", tolerance)
  }
}

if (all(c("pval_ge_CCT") %in% names(cell)) &&
    all(c("pval_ge_CCT") %in% names(donor))) {
  maximum <- max(abs(cell$pval_ge_CCT - donor$pval_ge_CCT), na.rm = TRUE)
  cat(sprintf("%-16s max absolute difference: %.3e\n", "pval_ge_CCT", maximum))
  if (!is.finite(maximum) || maximum > tolerance) {
    stop("pval_ge_CCT differs by more than ", tolerance)
  }
}

cat("PASS: donor and cell GxE conditioning engines agree.\n")
