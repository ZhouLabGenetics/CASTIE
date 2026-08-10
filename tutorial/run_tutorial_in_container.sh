#!/usr/bin/env bash
set -euo pipefail

TUTORIAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${TUTORIAL_DIR}"

# A source installation sets CASTIE_SOURCE_ROOT so the scripts under extdata/
# resolve exactly like the commands installed in the container image.
if [[ -n "${CASTIE_SOURCE_ROOT:-}" ]]; then
  export PATH="${CASTIE_SOURCE_ROOT}/extdata:${PATH}"
fi

mkdir -p output

step1_fitNULLGLMM_qtl.R \
  --useSparseGRMtoFitNULL=FALSE \
  --useGRMtoFitNULL=FALSE \
  --phenoFile=data/phenotypes.tsv \
  --phenoCol=gene_1 \
  --covarColList=X1,X2,pf1,pf2 \
  --sampleCovarColList=X1,X2 \
  --dynamicCovarColList=pf1,pf2 \
  --sampleIDColinphenoFile=IND_ID \
  --traitType=count \
  --outputPrefix=output/gene_1 \
  --skipVarianceRatioEstimation=FALSE \
  --isRemoveZerosinPheno=FALSE \
  --isCovariateOffset=FALSE \
  --isCovariateTransform=TRUE \
  --skipModelFitting=FALSE \
  --tol=0.00001 \
  --plinkFile=data/grm_variants \
  --IsOverwriteVarianceRatioFile=TRUE

step2_tests_qtl.R \
  --bedFile=data/genotypes.bed \
  --bimFile=data/genotypes.bim \
  --famFile=data/genotypes.fam \
  --SAIGEOutputFile=output/gene_1_cis \
  --chrom=1 \
  --minMAF=0 \
  --minMAC=1 \
  --LOCO=FALSE \
  --GMMATmodelFile=output/gene_1.rda \
  --SPAcutoff=2 \
  --varianceRatioFile=output/gene_1.varianceRatio.txt \
  --rangestoIncludeFile=data/gene_1_region.tsv \
  --markers_per_chunk=1000 \
  --output_format=txt

# Prepare the combined Step 2 table consumed by Step 3. For a multi-gene run,
# use the Polars concatenation workflow described in the documentation.
Rscript -e '
  library(data.table)
  result <- fread("output/gene_1_cis")
  dynamic_p <- tstrsplit(result$pval_ge, ",", fixed = TRUE)
  result[, pf1 := dynamic_p[[1]]]
  result[, pf2 := dynamic_p[[2]]]
  result <- result[AF_Allele2 > 0.05 & AF_Allele2 < 0.95,
                   .(Gene = "gene_1", MarkerID, AF_Allele2,
                     pval_main = p.value, pf1, pf2)]
  fwrite(result, "output/step3_input.txt", sep = "\t")
'

step3_gene_pvalue.R \
  --input=output/step3_input.txt \
  --outdir=output/step3

step4_get_egenes.R \
  --input=output/step3/step3_longformat.txt \
  --outdir=output/step4 \
  --fdr=0.05

echo "Tutorial complete. Results are in ${TUTORIAL_DIR}/output"
