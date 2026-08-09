#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(qvalue)
  library(optparse)
})

option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Path to Step 3 long-format input file [required]"),
  make_option("--outdir", type = "character", default = ".",
              help = "Output directory [default: current directory]"),
  make_option("--prefix", type = "character", default = "",
              help = "Optional output file prefix [default: none]"),
  make_option("--fdr", type = "numeric", default = 0.05,
              help = "FDR threshold for eGene calling [default: 0.05]"),
  make_option("--main_col", type = "character", default = "pval_main",
              help = "Column name for main eQTL p-values [default: pval_main]"),
  make_option("--exclude_cols", type = "character", default = "pval_main,pval_ge_CCT",
              help = "Comma-separated columns to exclude from context union [default: pval_main,pval_ge_CCT]")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input)) stop("--input is required")
if (!file.exists(opt$input)) stop("Input file does not exist: ", opt$input)
if (!is.finite(opt$fdr) || opt$fdr <= 0 || opt$fdr > 1) {
  stop("--fdr must be greater than 0 and no greater than 1")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
exclude_cols <- trimws(strsplit(opt$exclude_cols, ",")[[1]])

df <- fread(opt$input)
required_cols <- c("Gene", "pval_column", "ACAT_p")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols)) {
  stop("Input is missing required columns: ", paste(missing_cols, collapse = ", "))
}

df[, Gene := as.factor(Gene)]
df[, pval_column := as.factor(pval_column)]

cat("===== CASTIE Step 4 eGene calling start =====\n")

# Per-context eGene calling
egene_list <- list()

for (col in levels(df$pval_column)) {
  sub <- df[pval_column == col]
  valid <- is.finite(sub$ACAT_p) & sub$ACAT_p > 0 & sub$ACAT_p <= 1
  sub[, qv := NA_real_]
  if (any(valid)) {
    sub[valid, qv := qvalue(ACAT_p, pi0 = 1, lfdr.out = FALSE)$qvalue]
  }

  message(col, " qvalue summary:")
  print(summary(sub$qv))

  egenes <- sub[!is.na(qv) & qv <= opt$fdr][order(qv)]
  message(col, " eGenes at FDR ", opt$fdr, ": ", nrow(egenes))
  egene_list[[col]] <- as.character(egenes$Gene)

  outfile <- file.path(opt$outdir, paste0(opt$prefix, col, "_egene.tsv"))
  fwrite(egenes, outfile, sep = "\t")
  message("Written: ", outfile)

  txtfile <- file.path(opt$outdir, paste0(opt$prefix, col, "_egene_genes.txt"))
  writeLines(as.character(egenes$Gene), txtfile)
  message("Written: ", txtfile)
}

# Context-only eGenes
main_genes <- egene_list[[opt$main_col]]
if (is.null(main_genes)) {
  main_genes <- character()
  warning("Main p-value column not found: ", opt$main_col)
}
context_cols <- setdiff(names(egene_list), exclude_cols)

message("\nContext columns used for union: ", paste(context_cols, collapse = ", "))

all_ctx_genes <- unique(as.character(unlist(egene_list[context_cols], use.names = FALSE)))
ctx_only_genes <- setdiff(all_ctx_genes, main_genes)

ctx_union_file <- file.path(opt$outdir, paste0(opt$prefix, "all_contexts_egenes.txt"))
writeLines(all_ctx_genes, ctx_union_file)
message("\nAll-context union eGenes (", length(all_ctx_genes), "): ", ctx_union_file)

ctx_only_file <- file.path(opt$outdir, paste0(opt$prefix, "context_only_egenes.txt"))
writeLines(ctx_only_genes, ctx_only_file)
message("Context-only eGenes (not in ", opt$main_col, ") (", length(ctx_only_genes), "): ", ctx_only_file)

# Shared eGenes across contexts
message("\n-- Shared eGenes across contexts --")
ctx_egenes <- egene_list[context_cols]

if (length(ctx_egenes) >= 2) {
  for (i in seq_along(ctx_egenes)) {
    for (j in seq_along(ctx_egenes)) {
      if (j <= i) next
      shared <- intersect(ctx_egenes[[i]], ctx_egenes[[j]])
      if (length(shared) > 0) {
        message(names(ctx_egenes)[i], " & ", names(ctx_egenes)[j], ": ",
                length(shared), " shared -- ", paste(shared, collapse = ", "))
      }
    }
  }
}

gene_ctx_count <- sort(table(as.character(unlist(ctx_egenes, use.names = FALSE))), decreasing = TRUE)
multi_ctx <- gene_ctx_count[gene_ctx_count > 1]

if (length(multi_ctx) > 0) {
  message("\nGenes in >1 context:")
  print(multi_ctx)
  multi_ctx_dt <- data.table(Gene = names(multi_ctx), n_contexts = as.integer(multi_ctx))
  multi_ctx_dt[, contexts := vapply(Gene, function(g) {
    paste(names(ctx_egenes)[vapply(ctx_egenes, function(x) g %in% x, logical(1))], collapse = ",")
  }, character(1))]
  multi_ctx_file <- file.path(opt$outdir, paste0(opt$prefix, "shared_context_egenes.tsv"))
  fwrite(multi_ctx_dt, multi_ctx_file, sep = "\t")
  message("Written: ", multi_ctx_file)
} else {
  message("No genes shared across contexts")
}

message("\nCASTIE Step 4 complete.")
