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