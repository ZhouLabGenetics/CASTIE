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
  --pval_cutoff_for_gxe=1 \
  --varianceRatioFile=output/gene_1.varianceRatio.txt \
  --rangestoIncludeFile=data/gene_1_region.tsv \
  --markers_per_chunk=1000 \
  --output_format=txt

concat_step2_results.py \
  --input-dir=output \
  --output=output/step3_input.txt \
  --contexts=pf1,pf2 \
  --file-pattern='*_cis' \
  --gene-regex='^(?P<gene>.+)_cis$' \
  --maf-min=0.05 \
  --maf-max=0.95

step3_gene_pvalue.R \
  --input=output/step3_input.txt \
  --outdir=output/step3

step4_get_egenes.R \
  --input=output/step3/step3_longformat.txt \
  --outdir=output/step4 \
  --fdr=0.05

echo "Tutorial complete. Results are in ${TUTORIAL_DIR}/output"
