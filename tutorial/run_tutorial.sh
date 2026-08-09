#!/usr/bin/env bash
set -euo pipefail

TUTORIAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${CASTIE_IMAGE:-yijia0802/castie:Latest}"

python3 "${TUTORIAL_DIR}/generate_simulated_data.py"
mkdir -p "${TUTORIAL_DIR}/output"

docker run --rm --platform linux/amd64 \
  -v "${TUTORIAL_DIR}:/tutorial" \
  -w /tutorial \
  "${IMAGE}" \
  bash -lc '
    set -euo pipefail

    step1_fitNULLGLMM_qtl.R \
      --useSparseGRMtoFitNULL=FALSE \
      --useGRMtoFitNULL=FALSE \
      --phenoFile=data/phenotypes.tsv \
      --phenoCol=gene_1 \
      --covarColList=X1,X2,pf1,pf2 \
      --sampleCovarColList=X1,X2 \
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

    step3_gene_pvalue_qtl.R \
      --assocFile=output/gene_1_cis \
      --geneName=gene_1 \
      --genePval_outputFile=output/gene_1_gene_pvalue.tsv
  '

echo "Tutorial complete. Results are in ${TUTORIAL_DIR}/output"
