#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---------------- Command-line options ----------------
option_list <- list(
  make_option("--input", type="character", default="",
              help="Path to the main_results_filtered.txt file (input) [required]"),
  make_option("--outdir", type="character", default=".",
              help="Output directory [default='.']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt  <- args$options

if (opt$input == "")
  stop("Error: --input is required.\n")

# ---------------- Paths ----------------
infile  <- opt$input
outdir  <- opt$outdir
outfile <- file.path(outdir, "step3_longformat.txt")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Input file:", infile, "\n")
cat("Output directory:", outdir, "\n")

# ---------------- Define ACAT function (self-contained) ----------------
ACAT <- function(p) {
  p <- p[is.finite(p) & p > 0 & p < 1]
  if (!length(p)) return(NA_real_)
  t <- tan((0.5 - p) * pi)
  pmin(pmax(0.5 - atan(mean(t)) / pi, 0), 1)
}

# ---------------- Read input ----------------
data <- fread(infile, header = TRUE)

# Define columns
exclude_cols <- c("MarkerID", "AF_Allele2", "Gene")
pval_cols <- setdiff(names(data), exclude_cols)
cat("Processing p-value columns:", paste(pval_cols, collapse = ", "), "\n")

# Sanity check
if (!all(c("Gene","MarkerID") %in% names(data))) {
  stop("Need columns 'Gene' and 'MarkerID'. Got: ", paste(names(data), collapse=", "))
}

# ---------------- Convert p-columns to numeric ----------------
data[, (pval_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = pval_cols]

# ---------------- Clean p-values ----------------
EPS <- 1e-300
data[, (pval_cols) := lapply(.SD, function(p) {
  p[!is.finite(p)] <- NA_real_
  p[p < 0]  <- NA_real_
  p[p > 1]  <- NA_real_
  p[p == 0] <- EPS
  p
}), .SDcols = pval_cols]

# ---------------- Long format ----------------
long <- melt(
  data,
  id.vars = c("Gene","MarkerID"),
  measure.vars = pval_cols,
  variable.name = "pval_column",
  value.name   = "p"
)

# ---------------- ACAT per (Gene, pval_column) + top SNP ----------------
res <- long[, {
  p_clean <- p[is.finite(p) & p > 0 & p < 1]
  acat    <- ACAT(p_clean)
  idx     <- which.min(replace(p, is.na(p), Inf))
  list(
    ACAT_p       = acat,
    top_MarkerID = MarkerID[idx],
    top_pval     = p[idx]
  )
}, by = .(Gene, pval_column)][order(Gene, pval_column)]

# ---------------- Write per-gene files (with Gene column) ----------------
res[, {
  tmp <- copy(.SD)
  tmp[, Gene := .BY$Gene]
  fwrite(tmp, file.path(outdir, paste0(.BY$Gene, "_ACAT_results.txt")),
         sep = "\t", quote = FALSE)
  NULL
}, by = Gene]

# ---------------- Combine all per-gene files ----------------
cat("Combining per-gene files into one long-format file...\n")

gene_files <- list.files(outdir, pattern = "_ACAT_results\\.txt$", full.names = TRUE)

read_with_gene <- function(fp) {
  dt <- fread(fp)
  if (!"Gene" %in% names(dt)) {
    g <- sub("_ACAT_results\\.txt$", "", basename(fp))
    dt[, Gene := g]
  }
  dt
}

combined <- rbindlist(lapply(gene_files, read_with_gene), use.names = TRUE, fill = TRUE)
setcolorder(combined, c("Gene", setdiff(names(combined), "Gene")))

fwrite(combined, outfile, sep = "\t", quote = FALSE)
cat("Done! Combined long-format file saved to:", outfile, "\n")

