#!/usr/bin/env Rscript
#
# Convert a plain-text phenotype file (TSV/CSV) to HDF5 format for CASTIE.
#
# The HDF5 file stores:
#   /meta/<colname>    - one dataset per metadata column (IDs, covariates)
#   /counts/data       - non-zero gene expression values (int32)
#   /counts/indices    - row indices of non-zeros (int32, 0-based)
#   /counts/indptr     - column pointers (int64, length = ngenes + 1)
#   /counts/gene_names - gene name strings
#
# Usage:
#   Rscript convert_phenoFile_to_h5.R \
#     --phenoFile input.tsv \
#     --output output.h5 \
#     --metaCols "IID,cellID,age,sex,pc1,pc2,pc3,log_cell_read_counts"
#
# All columns NOT listed in --metaCols are treated as gene expression columns
# and stored in the sparse counts matrix.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(hdf5r)
})

option_list <- list(
  make_option("--phenoFile", type = "character", default = "",
    help = "Required. Path to the plain-text phenotype file (TSV/CSV, optionally gzipped)."),
  make_option("--output", type = "character", default = "",
    help = "Required. Path for the output HDF5 file."),
  make_option("--metaCols", type = "character", default = "",
    help = "Required. Comma-separated list of metadata column names (IDs, covariates, offsets). All other columns are treated as gene expression."),
  make_option("--sampleIDCol", type = "character", default = "IID",
    help = "Sample ID column name. Read as character to preserve leading zeros / non-numeric IDs [default: IID]."),
  make_option("--cellIDCol", type = "character", default = "",
    help = "Optional. Cell ID / barcode column name. When set, this column is read as character so barcodes like 'AAACCTGAGAAACCAT-1' aren't coerced to numeric / NA."),
  make_option("--chunkSize", type = "integer", default = 5000L,
    help = "Number of gene columns to process at a time [default: 5000]. Larger = faster but more memory.")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

if (opt$phenoFile == "" || opt$output == "" || opt$metaCols == "") {
  stop("--phenoFile, --output, and --metaCols are all required.\n")
}

metaColNames <- strsplit(opt$metaCols, ",")[[1]]

cat("=== Convert phenotype file to HDF5 ===\n")
cat("Input :", opt$phenoFile, "\n")
cat("Output:", opt$output, "\n")
cat("Meta columns:", paste(metaColNames, collapse = ", "), "\n")

## Step 1: read header to identify gene columns
isCompressed <- grepl("\\.(gz|bgz)$", opt$phenoFile)
if (isCompressed) {
  headerLine <- system(paste0("gunzip -c ", opt$phenoFile, " | head -n 1"), intern = TRUE)
} else {
  headerLine <- readLines(opt$phenoFile, n = 1)
}

if (grepl("\t", headerLine)) {
  allCols <- strsplit(headerLine, "\t", fixed = TRUE)[[1]]
} else if (grepl(",", headerLine)) {
  allCols <- strsplit(headerLine, ",", fixed = TRUE)[[1]]
} else {
  allCols <- strsplit(trimws(headerLine), "\\s+")[[1]]
}

geneCols <- setdiff(allCols, metaColNames)
cat("Total columns:", length(allCols), "\n")
cat("Metadata columns:", length(metaColNames), "\n")
cat("Gene columns:", length(geneCols), "\n")

## Step 2: read metadata columns
## force sample ID (and cell ID, if specified) as character to preserve
## barcodes / leading zeros / mixed-format IDs.
charCols <- opt$sampleIDCol
if (nzchar(opt$cellIDCol)) charCols <- c(charCols, opt$cellIDCol)

cat("Reading metadata columns...\n")
metaData <- fread(opt$phenoFile, header = TRUE, select = metaColNames,
  colClasses = list(character = charCols),
  stringsAsFactors = FALSE, data.table = FALSE)
nrows <- nrow(metaData)
cat("Rows:", nrows, "\n")

## Step 3: create HDF5 file and write metadata
if (file.exists(opt$output)) file.remove(opt$output)
h5 <- H5File$new(opt$output, mode = "w")

metaGrp <- h5$create_group("meta")
for (col in colnames(metaData)) {
  vals <- metaData[[col]]
  if (is.character(vals) || is.factor(vals)) {
    metaGrp[[col]] <- as.character(vals)
  } else {
    metaGrp[[col]] <- vals
  }
  cat("  Written meta/", col, " (", class(metaData[[col]])[1], ")\n", sep = "")
}
rm(metaData); gc()

## Step 4: build sparse CSC matrix in chunks
cat("Building sparse gene expression matrix...\n")
countsGrp <- h5$create_group("counts")

ngenes <- length(geneCols)
chunkSize <- opt$chunkSize

## pre-allocate lists for CSC components
## indptr / nnz_total stored as double so total nnz can exceed 2^31
## (R doubles give 53-bit integers, ample for any realistic dataset).
all_data    <- vector("list", ceiling(ngenes / chunkSize))
all_indices <- vector("list", ceiling(ngenes / chunkSize))
indptr      <- numeric(ngenes + 1L)
indptr[1]   <- 0
nnz_total   <- 0

chunk_i <- 0L
for (gStart in seq(1, ngenes, by = chunkSize)) {
  gEnd <- min(gStart + chunkSize - 1L, ngenes)
  chunk_i <- chunk_i + 1L
  chunkGenes <- geneCols[gStart:gEnd]

  cat(sprintf("  Chunk %d: genes %d-%d of %d...\n", chunk_i, gStart, gEnd, ngenes))
  chunkData <- fread(opt$phenoFile, header = TRUE, select = chunkGenes,
    stringsAsFactors = FALSE, data.table = FALSE)

  ## collect per-gene non-zero values/indices into lists, then concatenate
  ## once per chunk -- avoids quadratic c() growth across ~chunkSize iterations
  vals_list <- vector("list", length(chunkGenes))
  idx_list  <- vector("list", length(chunkGenes))

  for (j in seq_along(chunkGenes)) {
    vec_raw <- chunkData[[j]]
    ## reject non-integer counts (e.g., normalized / log-transformed input)
    if (!is.integer(vec_raw)) {
      bad <- which(!is.na(vec_raw) & vec_raw != as.integer(vec_raw))
      if (length(bad) > 0L) {
        stop("Non-integer values found in gene column '", chunkGenes[j],
             "' (e.g., row ", bad[1], " = ", vec_raw[bad[1]],
             "). The converter expects raw integer counts.")
      }
    }
    vec <- as.integer(vec_raw)
    nz <- which(vec != 0L)
    vals_list[[j]] <- vec[nz]
    idx_list[[j]]  <- nz - 1L  # 0-based row indices
    nnz_total <- nnz_total + length(nz)
    indptr[gStart + j] <- nnz_total
  }

  all_data[[chunk_i]]    <- unlist(vals_list, use.names = FALSE)
  all_indices[[chunk_i]] <- unlist(idx_list,  use.names = FALSE)

  rm(chunkData, vals_list, idx_list); gc()
}

cat(sprintf("Total non-zeros: %d (%.1f%% sparse)\n",
  nnz_total, 100 * (1 - nnz_total / (as.numeric(nrows) * ngenes))))

## Step 5: write CSC arrays to HDF5
cat("Writing sparse matrix to HDF5...\n")
all_data_vec    <- unlist(all_data, use.names = FALSE)
all_indices_vec <- unlist(all_indices, use.names = FALSE)
rm(all_data, all_indices); gc()

countsGrp[["data"]]       <- all_data_vec
countsGrp[["indices"]]    <- all_indices_vec
countsGrp[["indptr"]]     <- indptr   # double, supports nnz > 2^31
countsGrp[["gene_names"]] <- geneCols

rm(all_data_vec, all_indices_vec); gc()

h5$close()

outSize <- file.info(opt$output)$size / 1024 / 1024
cat(sprintf("\nDone! Output: %s (%.1f MB)\n", opt$output, outSize))
cat("Use --phenoFile", opt$output, "in step1_fitNULLGLMM_qtl.R\n")
