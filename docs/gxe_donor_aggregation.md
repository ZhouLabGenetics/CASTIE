# Donor-aggregate exact GxC conditioning

CASTIE's exact model-based GxC test conditions the context-specific score on
the marginal genotype score. The original implementation performs the exact
calculation with cell-length vectors for every variant and every context. The
donor engine implements the same statistic from sufficient statistics cached
once per null model.

No Step 1 change or new Step 1 model is required. The main-effect test and its
variance ratio are unchanged.

## Calculation

For donor \(d\), let \(g_d\) be the donor genotype, \(c_k\) context \(k\),
\(X\) the cell-level nuisance design, and \(P\) the exact null-model
projection under \(\Sigma\). CASTIE tests the score conditioned on the
marginal genotype score:

\[
I_{11}=G^T P G,\qquad
I_{21,k}=(G C_k)^T P G,\qquad
I_{22,k}=(G C_k)^T P(G C_k),
\]

\[
T_{c,k}=T_{GC,k}-\frac{I_{21,k}}{I_{11}}T_G,
\qquad
A_{c,k}=I_{22,k}-\frac{I_{21,k}^2}{I_{11}}.
\]

Because genotype is constant within donor, each quadratic form is a sum over
donors. CASTIE caches the donor contributions to \(I_{11}\), \(I_{21,k}\),
\(I_{22,k}\), and the nuisance projections \(X^T\Sigma^{-1}G\) and
\(X^T\Sigma^{-1}(GC_k)\). A marker then requires donor-length rather than
cell-length operations.

When Step 1 selected sandwich variance for a context, the donor engine uses
the same corrected score and donor-clustered meat as the reference engine:

\[
Q_{d,k}=\sum_{i\in d}r_i
\left(\widetilde{GC}_{ik}-\frac{I_{21,k}}{I_{11}}\widetilde G_i\right),
\]

\[
V_{SW,k}=\sum_d(Q_{d,k}-\bar Q_k)^2\frac{D}{D-p}.
\]

The reported effect and standard error remain

\[
\widehat\beta_k=T_{c,k}/A_{c,k},\qquad
SE_k=\sqrt{V_k}/A_{c,k},
\]

where \(V_k=A_{c,k}\) for model variance and \(V_k=V_{SW,k}\) for sandwich
variance.

## Usage

The donor engine is the default for exact model conditioning:

```bash
step2_tests_qtl.R \
  --gxe_conditioning=model \
  --gxe_modelcond_engine=donor \
  ...
```

The previous cell-level implementation remains available as a reference:

```bash
step2_tests_qtl.R \
  --gxe_conditioning=model \
  --gxe_modelcond_engine=cell \
  ...
```

Permutation modes currently select the cell engine automatically. Models fit
with a GRM continue to fall back to legacy variance-ratio conditioning because
the donor-block closed form is not applicable.

## Validation

Run Step 2 once with each engine and compare the outputs:

```bash
Rscript tutorial/validate_gxe_donor_engine.R \
  output.cell.txt \
  output.donor.txt
```

Small differences in the final displayed digit of a combined CCT p-value can
occur because the mathematically equivalent engines accumulate floating-point
sums in a different order.
