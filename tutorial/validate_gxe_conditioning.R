#!/usr/bin/env Rscript
## Independent check of the step-2 dynamic (GxE) test conditioned on the marginal G score.
##
## Recomputes, in plain R, the conditional GxE chi-square for every tested variant from the
## step-1 model file and the genotypes, using
##   model    : Var(T_G), Cov(T_GE, T_G), Var(T_GE) = u' P v with the null-model covariance
##              Sigma_d = tau0 W_d^-1 + tau1 11' + tau2 E_d E_d' (each donor block inverted directly)
##   varratio : the legacy approximation  VR_g * <g~,g~>_W,  sqrt(VR_eg VR_g) * <ge~,g~>_W,
##              VR_eg * <ge~,ge~>_W  (variance ratios from step 1)
## and compares them with step-2 output run with --gxe_conditioning=model / =varratio.
## Intended for small data (dense per-donor inverses), e.g. the tutorial.
##
## Usage: Rscript validate_gxe_conditioning.R <model.rda> <plink prefix> <varianceRatio.txt> \
##          <step2 output, model> <step2 output, varratio>
args <- commandArgs(TRUE)
stopifnot(length(args) == 5)
load(args[1])                     # modglmm
m <- modglmm

## ---- null model pieces ----
X    <- t(m$obj.noK$XV) / as.vector(m$obj.noK$V)        # cell-level design (as used in step 1)
mu   <- as.vector(m$fitted.values)
vw   <- as.vector(m$varWeights)
W    <- mu * vw                                          # count trait: mu.eta^2/variance * var_weights
res  <- as.vector(m$residuals)
tau  <- as.vector(m$theta); tau <- c(tau, 0, 0)[1:3]
E    <- as.matrix(m$eMat)                                # standardized contexts (GxE and random slopes)
did  <- match(m$sampleID, unique(m$sampleID))
stopifnot(m$traitType == "count", nrow(X) == length(mu))

blocks <- split(seq_along(did), did)
SigmaInv <- function(V) {                                # dense per-donor inverse, no Woodbury
  V <- as.matrix(V); out <- V
  for (ix in blocks) {
    Sd <- tau[1] * diag(1 / W[ix], length(ix)) + tau[2] + tau[3] * tcrossprod(E[ix, , drop = FALSE])
    out[ix, ] <- solve(Sd, V[ix, , drop = FALSE])
  }
  out
}
SiX <- SigmaInv(X); XSXinv <- solve(crossprod(X, SiX))
uPv <- function(u, v, Siv) sum(u * Siv) - drop(crossprod(crossprod(SiX, u), XSXinv %*% crossprod(SiX, v)))
XVXinv_XV <- solve(crossprod(X, as.vector(m$obj.noK$V) * X), t(X * as.vector(m$obj.noK$V)))
Xadj <- function(v) v - drop(X %*% (XVXinv_XV %*% v))

## ---- variance ratios (as parsed in readInGLMM) ----
vr <- read.table(args[3], header = FALSE, stringsAsFactors = FALSE)
VR_g  <- as.numeric(vr[vr[, 2] == "null" & vr[, 3] > 0, 1])[1]
VR_eg <- as.numeric(vr[vr[, 2] == "null" & vr[, 3] == 0, 1])

## ---- genotypes (PLINK 1 bed, SNP-major) ----
pfx <- args[2]
fam <- read.table(paste0(pfx, ".fam"), stringsAsFactors = FALSE)
bim <- read.table(paste0(pfx, ".bim"), stringsAsFactors = FALSE)
nS <- nrow(fam); nb <- ceiling(nS / 4)
raw <- readBin(paste0(pfx, ".bed"), "raw", 3 + nb * nrow(bim))[-(1:3)]
decode <- function(j) {
  by <- as.integer(raw[(j - 1) * nb + seq_len(nb)])
  code <- as.vector(sapply(by, function(b) bitwAnd(bitwShiftR(b, c(0, 2, 4, 6)), 3L)))[1:nS]
  c(2, NA, 1, 0)[code + 1]                               # 00=hom A1, 01=missing, 10=het, 11=hom A2 -> A1 dosage
}
cell2fam <- match(m$sampleID, fam$V2)

## ---- recompute conditional GxE chi-square per variant and context ----
use_sw <- if (is.null(m$use_sandwich)) rep(FALSE, ncol(E)) else as.logical(m$use_sandwich)
nD <- length(blocks); pX <- ncol(X)
eff_res <- res * vw / tau[1]
crve <- function(score_cell) {                           # donor-clustered variance, as in step 2
  Q <- as.vector(rowsum(score_cell, did)); sum((Q - mean(Q))^2) * nD / (nD - pX)
}
pv <- function(Tc, Vc) if (Vc <= 1e-5) 1 else pchisq(Tc^2 / Vc, 1, lower.tail = FALSE)
rows <- list()
for (j in seq_len(nrow(bim))) {
  gd <- decode(j); if (anyNA(gd)) next
  g <- gd[cell2fam]
  if (mean(g) / 2 > 0.5) g <- 2 - g                     # chi-square is invariant; mirrors CASTIE's flip
  gt <- Xadj(g); if (sum(gt^2) < 1e-8) next
  TG <- sum(gt * res * vw) / tau[1]
  SiG <- SigmaInv(gt); I11 <- uPv(gt, gt, SiG)
  dG <- sum(m$obj.noK$V * gt^2)
  for (k in seq_len(ncol(E))) {
    get <- Xadj(E[, k] * g)
    TGE <- sum(get * res * vw) / tau[1]
    I21 <- uPv(get, gt, SiG); I22 <- uPv(get, get, SigmaInv(get))
    cm <- I21 / I11
    Tc_m <- TGE - cm * TG
    Vc_m <- if (use_sw[k]) crve(eff_res * (get - cm * gt)) else I22 - cm * I21
    cr <- sqrt(VR_eg[k] * VR_g) * sum(m$obj.noK$V * get * gt)
    VGE <- VR_eg[k] * sum(m$obj.noK$V * get^2); VG <- VR_g * dG
    Tc_v <- TGE - cr / VG * TG
    Vc_v <- if (use_sw[k]) {
      crve(get * (eff_res - sqrt(VR_eg[k]) * sqrt(VR_g) * m$obj.noK$V * tau[1] * gt * TG / VG))
    } else VGE - cr^2 / VG
    rows[[length(rows) + 1]] <- data.frame(MarkerID = bim$V2[j], context = colnames(E)[k],
                                           p_model_R = pv(Tc_m, Vc_m), p_varratio_R = pv(Tc_v, Vc_v))
  }
}
ref <- do.call(rbind, rows)

## ---- compare with CASTIE step-2 output ----
read_step2 <- function(f) {                             # pval_ge holds one comma-separated value per context
  d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  pm <- do.call(rbind, strsplit(as.character(d$pval_ge), ","))
  do.call(rbind, lapply(seq_len(ncol(E)), function(k)
    data.frame(MarkerID = d$MarkerID, context = colnames(E)[k], p = as.numeric(pm[, k]))))
}
cmp <- merge(ref, setNames(read_step2(args[4]), c("MarkerID", "context", "p_model_CASTIE")))
cmp <- merge(cmp, setNames(read_step2(args[5]), c("MarkerID", "context", "p_varratio_CASTIE")))
lp <- function(p) -log10(pmax(p, 1e-300))
cat(sprintf("compared %d variant x context tests (%d variants, contexts: %s)\n",
            nrow(cmp), length(unique(cmp$MarkerID)), paste(unique(cmp$context), collapse = ", ")))
cat(sprintf("model-based : max |log10 p (CASTIE) - log10 p (R)| = %.2e\n",
            max(abs(lp(cmp$p_model_CASTIE) - lp(cmp$p_model_R)))))
cat(sprintf("varratio    : max |log10 p (CASTIE) - log10 p (R)| = %.2e\n",
            max(abs(lp(cmp$p_varratio_CASTIE) - lp(cmp$p_varratio_R)))))
cat(sprintf("model vs varratio (CASTIE): median |delta log10 p| = %.3f, max = %.3f\n",
            median(abs(lp(cmp$p_model_CASTIE) - lp(cmp$p_varratio_CASTIE))),
            max(abs(lp(cmp$p_model_CASTIE) - lp(cmp$p_varratio_CASTIE)))))
