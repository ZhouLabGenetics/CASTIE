# CASTIE tutorial with simulated data

This tutorial runs the four-stage CASTIE workflow on a small deterministic
dataset. It creates 40 simulated donors, five cells per donor, 200 variants,
cell- and donor-level covariates, and one count phenotype (`gene_1`). The data
are for software demonstration only and have no biological interpretation.

## Requirements

- Python 3
- Docker Desktop, or another Docker-compatible engine
- The `yijia0802/castie:Latest` image available locally

## Run the tutorial

From the repository root:

```bash
bash tutorial/run_tutorial.sh
```

The script regenerates `tutorial/data/`, mounts the tutorial directory into the
container, and writes results to `tutorial/output/`.

To use another compatible image:

```bash
CASTIE_IMAGE=your-account/castie:tag bash tutorial/run_tutorial.sh
```

## Generate only the input data

```bash
python3 tutorial/generate_simulated_data.py
```

Generated inputs:

- `phenotypes.tsv`: cell-level phenotype and covariates
- `genotypes.{bed,bim,fam}`: PLINK 1 genotype dataset for association testing
- `grm_variants.{bed,bim,fam}`: PLINK 1 dataset used in Step 1
- `gene_1_region.tsv`: chromosome interval tested in Step 2

## Workflow

1. Step 1 fits the null generalized linear mixed model and estimates the
   variance ratio.
2. Step 2 tests variants in the simulated cis region.
3. Step 3 combines variant-level p-values into a gene-level p-value.
4. Step 4 applies FDR control to call eGenes. The tutorial converts its
   one-gene Step 3 result to the same long format used in a multi-gene run.

The complete commands are kept in `run_tutorial.sh` so they can be copied and
adapted to real datasets.

## Step 4 output

The tutorial runs Step 4 automatically. Its inputs and outputs are:

- `output/step3_longformat.txt`: one-gene Step 3 result in the required schema
- `output/step4/pval_main_egene.tsv`: significant main-effect eGenes
- `output/step4/pval_main_egene_genes.txt`: significant gene names
- `output/step4/all_contexts_egenes.txt`: union of context-dependent eGenes
- `output/step4/context_only_egenes.txt`: context-dependent eGenes absent from
  the main-effect list

For a real multi-gene dynamic analysis, create `step3_longformat.txt` with
`step3_0.2.5.7.R`, then use the same Step 4 command:

```bash
step4_get_egenes.R \
  --input output/step3_longformat.txt \
  --outdir output/step4 \
  --fdr 0.05
```

Step 4 expects the columns `Gene`, `pval_column`, and `ACAT_p`. It writes an
eGene table and gene list for every p-value context, plus context-union,
context-only, and shared-context summaries.
