SAIGE-QTL is an R package developed with Rcpp for scalable and accurate expression quantitative trait locus mapping for single-cell studies 

The method
- Model repeat and complex data structure, due to multiple cells per individual and relatedness between individuals 
- Model discrete read counts
- Fast and scalable for large data, test 20k genes, 200 cell types, millions of cells, millions of variants
- Test rare variations. Single-variant test is underpowered

Please see [https://weizhou0.github.io/SAIGE-QTL-doc/](https://weizhou0.github.io/SAIGE-QTL-doc/) for how to run SAIGE-QTL.

## Version Updates

### Version 0.2.5.1

**New Features:**
- **Custom G×E Permutation Support**: Added `--permute_ginge_fam_file` parameter to Step2 for using custom genotype permutation orders from FAM files
  - When used with `--is_permute_ginge=TRUE`, allows specifying a pre-generated permuted FAM file instead of random permutation
  - Ensures consistent genotype permutation order between simulation and analysis workflows
  - Maintains backward compatibility - without the parameter, original random permutation behavior is preserved

**Usage:**
```bash
# Original random permutation (unchanged)
step2_tests_qtl.R --is_permute_ginge=TRUE

# New custom permutation using FAM file
step2_tests_qtl.R --is_permute_ginge=TRUE --permute_ginge_fam_file=path/to/permuted.fam
```
## Version Updates

### Version 0.2.5.2

- Updated to pixi installation

### Version 0.2.5.7

**Step 1 (`step1_fitNULLGLMM_qtl.R`) — new flags:**

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--smwCacheMemLimitMB` | numeric | NA | Skip the SMW cache and fall back to per-vector Sigma solves when estimated cache size exceeds this many MB. Useful to bound peak memory on large datasets |
| `--isWriteReport` | logical | FALSE | Save a fitting report (solver used, convergence, offset flag) to a `report/` subdirectory next to `outputPrefix` |

**Step 2 (`step2_tests_qtl.R`) — new flags:**

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--output_format` | character | parquet | Output format for association results: `parquet` or `txt` |

**Usage examples:**
```bash
# Step 1 with memory-bounded SMW cache and fitting report
step1_fitNULLGLMM_qtl.R --smwCacheMemLimitMB=8000 --isWriteReport=TRUE ...

# Step 2 writing plain text output instead of parquet
step2_tests_qtl.R --output_format=txt ...
```