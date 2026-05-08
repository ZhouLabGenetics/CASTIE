#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios
#' @param phenoFile character. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarCol. All categorical covariates must also be included in covarColList.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary" or "quantitative". By default, "binary"
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. e.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the null model. All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex").
#' @param eCovarCol vector of characters. Covariates of environmental factors/cell context to be used in the null model. All covariates listed in eCovarCol must be also in covarColList,  e,g c("cellType").
#' @param sampleIDColinphenoFile character. Column name for the sample IDs in the phenotype file e.g. "IID".
#' @param tol numeric. The tolerance for fitting the null model to converge. By default, 0.02.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param tolPCG numeric. The tolerance for PCG to converge. By default, 1e-5.
#' @param maxiterPCG integer. The maximum number of iterations for PCG. By default, 500.
#' @param nThreads integer. Number of threads to be used. By default, 1
#' @param SPAcutoff numeric. The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2.
#' @param numMarkersForVarRatio integer (>0). Minimum number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param tauInit vector of numbers. e.g. c(1,1), Initial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is 0,0, the second element will be 0.5 for binary traits and the initial tau vector for quantitative traits is 1,0
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. By default, TRUE
#' @param traceCVcutoff numeric. The threshold for coefficient of variation (CV) for the trace estimator to increase nrun. By default, 0.0025
#' @param ratioCVcutoff numeric. The threshold for coefficient of variation (CV) for the variance ratio estimate. If ratioCV > ratioCVcutoff. numMarkersForVarRatio will be increased by 10. By default, 0.001
#' @param outputPrefix character. Path to the output files with prefix.
#' @param outputPrefix_varRatio character. Path to the output variance ratio file with prefix. variace ratios will be output to outputPrefix_varRatio.varianceRatio.txt. If outputPrefix_varRatio is not specified, outputPrefix_varRatio will be the same as the outputPrefix. By default, ""
#' @param IsOverwriteVarianceRatioFile logical. Whether to overwrite the variance ratio file if the file exists. By default, FALSE
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. By default, ""
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. By default, ""
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples that are included in the sparse GRM. By default, 2000
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated in the sparse GRM.
#' @param cateVarRatioIndexVec vector of integer 0 or 1. The length of cateVarRatioIndexVec is the number of MAC categories for variance ratio estimation. 1 indicates variance ratio in the MAC category is to be estimated, otherwise 0. By default, NULL. If NULL, variance ratios corresponding to all specified MAC categories will be estimated. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param isDiagofKinSetAsOne logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param useSparseGRMtoFitNULL logical. Whether to use sparse GRM to fit the null model. By default, FALSE
#' @param useSparseGRMforVarRatio logical. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE
#' @param minCovariateCount integer. If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues. By default, -1 (no covariates will be excluded)
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @param includeNonautoMarkersforVarRatio logical. Whether to allow for non-autosomal markers for variance ratio. By default, FALSE
#' @param FemaleOnly logical. Whether to run Step 1 for females only. If TRUE, sexCol and FemaleCode need to be specified. By default, FALSE
#' @param MaleOnly logical. Whether to run Step 1 for males only. If TRUE, sexCol and MaleCode need to be specified. By default, FALSE
#' @param FemaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for females. By default, '1'
#' @param MaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for males. By default, '0'
#' @param sexCol character. Coloumn name for sex in the phenotype file, e.g Sex. By default, ''
#' @param isCovariateOffset logical. Whether to estimate fixed effect coeffciets. By default, FALSE.
#' @param isStoreSigma logical. Whether to store sigma matrix. By default, FALSE. If number of individuals is greater than 10,000, this option may use large memory
#' @param isShrinkModelOutput logical. remove unnecessary objects for step2 from the model output. By default, FALSE.
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM_multiV <- function(plinkFile = "",
                               bedFile = "",
                               bimFile = "",
                               famFile = "",
                               phenoFile = "",
                               phenoCol = "",
                               isRemoveZerosinPheno = FALSE,
                               traitType = "binary",
                               invNormalize = FALSE,
                               covarColList = NULL,
                               qCovarCol = NULL,
                               eCovarCol = NULL,
                               sampleCovarCol = NULL,
                               offsetCol = NULL,
                               varWeightsCol = NULL,
                               longlCol = "",
                               sampleIDColinphenoFile = "",
                               tol = 0.02,
                               maxiter = 20,
                               tolPCG = 1e-5,
                               maxiterPCG = 500,
                               nThreads = 1,
                               SPAcutoff = 2,
                               numMarkersForVarRatio = 30,
                               skipModelFitting = FALSE,
                               memoryChunk = 2,
                               tauInit = c(0, 0),
                               LOCO = TRUE,
                               isLowMemLOCO = FALSE,
                               traceCVcutoff = 0.0025,
                               ratioCVcutoff = 0.001,
                               outputPrefix = "",
                               outputPrefix_varRatio = "",
                               IsOverwriteVarianceRatioFile = FALSE,
                               sparseGRMFile = "",
                               sparseGRMSampleIDFile = "",
                               numRandomMarkerforSparseKin = 1000,
                               relatednessCutoff = 0.125,
                               isCateVarianceRatio = FALSE,
                               cateVarRatioIndexVec = NULL,
                               cateVarRatioMinMACVecExclude = c(10, 20.5),
                               cateVarRatioMaxMACVecInclude = c(20.5),
                               isCovariateTransform = TRUE,
                               isDiagofKinSetAsOne = FALSE,
                               minCovariateCount = -1,
                               minMAFforGRM = 0.01,
                               maxMissingRateforGRM = 0.15,
                               useSparseGRMtoFitNULL = FALSE,
                               useSparseGRMforVarRatio = FALSE,
                               includeNonautoMarkersforVarRatio = FALSE,
                               sexCol = "",
                               FemaleCode = 1,
                               FemaleOnly = FALSE,
                               MaleCode = 0,
                               MaleOnly = FALSE,
                               SampleIDIncludeFile = "",
                               isCovariateOffset = FALSE,
                               skipVarianceRatioEstimation = FALSE,
                               nrun = 30,
                               VmatFilelist = "",
                               VmatSampleFilelist = "",
                               VcellmatFilelist = "",
                               VcellmatSampleFilelist = "",
                               useGRMtoFitNULL = TRUE,
                               isStoreSigma = FALSE,
                               isShrinkModelOutput = FALSE,
                               smwCacheMemLimitMB = NULL,
                               verbose = FALSE) {
  ## set up output files
  modelOut <- paste0(outputPrefix, ".rda")

  if (skipModelFitting) {
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
  } else {
    if (LOCO & isLowMemLOCO) {
      modelOut <- paste(c(outputPrefix, "_noLOCO.rda"), collapse = "")
    }
    file.create(modelOut, showWarnings = TRUE)
  }
  if (plinkFile != "") {
    bimFile <- paste0(plinkFile, ".bim")
    bedFile <- paste0(plinkFile, ".bed")
    famFile <- paste0(plinkFile, ".fam")
  }
  setgenoNULL()

  if (!useGRMtoFitNULL) {
    useSparseGRMtoFitNULL <- FALSE
    useSparseGRMforVarRatio <- FALSE
    LOCO <- FALSE
    #nThreads <- 1
    #cat("No GRM will be used to fit the NULL model and nThreads is set to 1\n")
  }


  if (useSparseGRMtoFitNULL & bedFile == "") {
    cat("Sparse GRM is used to fit the null model and plink file is not specified, so variance ratios won't be estimated\n")
    skipVarianceRatioEstimation <- TRUE
  }

  if (!skipVarianceRatioEstimation) {
    SPAGMMATOut <- paste0(outputPrefix, "_", numMarkersForVarRatio, "markers.SAIGE.results.txt")
    # Check_OutputFile_Create(SPAGMMATOut)

    if (outputPrefix_varRatio == "") {
      outputPrefix_varRatio <- outputPrefix
    }
    varRatioFile <- paste0(outputPrefix_varRatio, ".varianceRatio.txt")

    if (!file.exists(varRatioFile)) {
      file.create(varRatioFile, showWarnings = TRUE)
    } else {
      if (!IsOverwriteVarianceRatioFile) {
        stop(
          "WARNING: The variance ratio file ", varRatioFile,
          " already exists. The new variance ratios will be output to ",
          varRatioFile, ". In order to avoid overwriting the file, please remove the ",
          varRatioFile, " or use the argument outputPrefix_varRatio to specify a different prefix to output the variance ratio(s). Otherwise, specify --IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten with new variance ratio(s)\n"
        )
      } else {
        cat("The variance ratio file ", varRatioFile, " already exists. IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten\n")
      }
    }
  } else {
    cat("Variance ratio estimation will be skipped\n.")
    useSparseGRMforVarRatio <- FALSE
  }


  if (useSparseGRMtoFitNULL) {
    # useSparseGRMforVarRatio = FALSE
    LOCO <- FALSE
    nThreads <- 1
    if (bedFile != "") {
      cat("sparse GRM will be used to fit the NULL model and nThreads is set to 1\n")
    }
    cat("Leave-one-chromosome-out is not applied\n")
  }



  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    if (!file.exists(sparseGRMFile)) {
      stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
    }
    if (!file.exists(sparseGRMSampleIDFile)) {
      stop(
        "sparseGRMSampleIDFile ", sparseGRMSampleIDFile,
        " does not exist!"
      )
    }
  }


  RcppParallel:::setThreadOptions(numThreads = nThreads)
  set_g_omp_num_threads(nThreads)
  data.table::setDTthreads(nThreads)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(nThreads)
    RhpcBLASctl::omp_set_num_threads(nThreads)
  }
  cat(nThreads, " thread(s) will be used\n")

  ## Memory-checkpoint helper.  Called at each [TIMING] phase boundary.
  ## Reports both R-side gc() (used / running peak) and OS-level RSS
  ## (which captures C++/Armadillo allocations that R's gc doesn't see).
  ## R peak is cumulative since session start; OS RSS is current at call time.
  ## Column count of gc() varies (6 vs 7) by R version / memory-limit setting,
  ## so use `ncol` for the trailing "max used (Mb)" column and `2` for "used (Mb)"
  ## (always the second column).
  .report_mem <- function(label) {
    if (!verbose) return(invisible(NULL))
    .gc <- gc(verbose = FALSE)
    rss_mb <- tryCatch(
      as.numeric(system2("ps", c("-o", "rss=", "-p", Sys.getpid()),
                         stdout = TRUE, stderr = FALSE)) / 1024,
      error = function(e) NA_real_
    )
    cat(sprintf("[MEM] %s: R-used=%.0f MB, R-peak=%.0f MB, OS-RSS=%.0f MB\n",
                label, sum(.gc[, 2]), sum(.gc[, ncol(.gc)]),
                if (length(rss_mb) == 0 || is.na(rss_mb)) NA_real_ else rss_mb))
  }

  ## Estimate the SMW cache size in MB and decide whether to populate it.
  ## Cache is per-donor and dominated by 4 N x k floats (U_ind + A_inv_U_ind),
  ## plus 2 N floats (w_ind + w_inv_tau0) and D x k^2 floats (middle_inv).
  ## Returns list(estMB, allowed) where `allowed=FALSE` means caller should skip
  ## prepareSigmaInvSMW_multiV and use the per-vector fallback.
  .smw_cache_decision <- function(eMat, label) {
    if (is.null(eMat)) return(list(estMB = 0, allowed = FALSE))
    N <- nrow(eMat); k <- ncol(eMat)
    n_donors <- if (exists("modglmm", inherits = FALSE) &&
                    !is.null(modglmm$sampleID))
                  length(unique(modglmm$sampleID)) else NA_real_
    bytes <- 4 * (N * (2 + 2 * k) + (if (is.na(n_donors)) 0 else n_donors) * k^2)
    estMB <- bytes / (1024 * 1024)
    over  <- !is.null(smwCacheMemLimitMB) && estMB > smwCacheMemLimitMB
    cat(sprintf("[SMW] %s: estimated cache size = %.0f MB (N=%d, k=%d%s)%s\n",
                label, estMB, N, k,
                if (is.na(n_donors)) "" else sprintf(", D=%d", as.integer(n_donors)),
                if (over) sprintf(" -> EXCEEDS smwCacheMemLimitMB=%.0f, skipping",
                                  smwCacheMemLimitMB) else ""))
    list(estMB = estMB, allowed = !over)
  }

  if (FemaleOnly & MaleOnly) {
    stop("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job\n")
  }

  if (FemaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_FemaleOnly")
    cat(
      "Female-specific model will be fitted. Samples coded as ",
      FemaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  } else if (MaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_MaleOnly")
    cat(
      "Male-specific model will be fitted. Samples coded as ",
      MaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  }


  sampleListwithGeno <- NULL
  if ((!useSparseGRMtoFitNULL & useGRMtoFitNULL) | !skipVarianceRatioEstimation) {
    if (!file.exists(bedFile)) {
      stop("ERROR! bed file does not exsit\n")
    }
    if (!file.exists(bimFile)) {
      stop("ERROR! bim file does not exsit\n")
    } else {
      if (LOCO) {
        chrVec <- data.table:::fread(bimFile, header = F, data.table = F, select = 1)
        updatechrList <- updateChrStartEndIndexVec(chrVec)
        LOCO <- updatechrList$LOCO
        chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
        chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
      }
      if (!LOCO) {
        chromosomeStartIndexVec <- rep(NA, 22)
        chromosomeEndIndexVec <- rep(NA, 22)
      }
    }


    if (!file.exists(famFile)) {
      stop("ERROR! fam file does not exsit\n")
    } else {
      sampleListwithGenov0 <- data.table:::fread(famFile, header = F, , colClasses = list(character = 1:4))
      sampleListwithGenov0 <- data.frame(sampleListwithGenov0)
      colnames(sampleListwithGenov0) <- c(
        "FIDgeno", "IIDgeno",
        "father", "mother", "sex", "phe"
      )
      sampleListwithGeno <- NULL
      sampleListwithGeno$IIDgeno <- sampleListwithGenov0$IIDgeno
      sampleListwithGeno <- data.frame(sampleListwithGeno)
      sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
        by = 1
      )
      cat(nrow(sampleListwithGeno), " samples have genotypes\n")
    }
  } else {
    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
      sampleListwithGenov0 <- data.table:::fread(sparseGRMSampleIDFile,
        header = F, , colClasses = c("character"), data.table = F
      )
      colnames(sampleListwithGenov0) <- c("IIDgeno")
      sampleListwithGeno <- NULL
      sampleListwithGeno$IIDgeno <- sampleListwithGenov0$IIDgeno
      sampleListwithGeno <- data.frame(sampleListwithGeno)
      sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
        by = 1
      )
      cat(nrow(sampleListwithGeno), " samples are in the sparse GRM\n")
    }
  }
  .t0_total <- proc.time()

  if (!file.exists(phenoFile)) {
    stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
  } else {

    if (length(eCovarCol) > 0) {
      cat(eCovarCol, "are environmental covariates\n")
      if (!all(eCovarCol %in% covarColList)) {
        stop("ERROR! all covariates in eCovarCol must be in covarColList\n")
      }
      longlCol = ""
    }

    if (longlCol == "") {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile)
    } else {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile, longlCol)
    }

    if (length(offsetCol) > 0) {
      cat(offsetCol, "is the offset term\n")
      checkColList <- c(checkColList, offsetCol)
    }

    if (length(varWeightsCol) > 0) {
      cat(varWeightsCol, " is the weights for variance\n")
      checkColList <- c(checkColList, varWeightsCol)
    }


    .t_pheno <- proc.time()
    isPhenoHDF5 <- grepl("\\.(h5|hdf5)$", phenoFile, ignore.case = TRUE)

    if (isPhenoHDF5) {
      ## ---- HDF5 format: read metadata + one gene column from sparse matrix ----
      data <- read_pheno_h5(phenoFile, phenoCol, sampleIDColinphenoFile, checkColList)

    } else {
      ## ---- Plain text format (tsv/csv/gz/bgz) ----
      isCompressed <- grepl("\\.(gz|bgz)$", phenoFile)

      ## detect file size to decide strategy
      if (isCompressed) {
        fileSize_KB <- as.numeric(system(paste0("du -k ", phenoFile, " | cut -f1"), intern = TRUE))
        isphenoFileLarge <- (fileSize_KB > 200000)
      } else {
        fileSize_KB <- file.info(phenoFile)$size / 1024
        isphenoFileLarge <- (fileSize_KB > 500000)
      }

      ## pick decompression command: pigz (parallel) > gunzip
      if (isCompressed) {
        if (system("command -v pigz > /dev/null 2>&1") == 0) {
          dcatCmd <- "pigz -dc"
        } else {
          dcatCmd <- "gunzip -c"
        }
      }

      if (isphenoFileLarge) {
        ## large file: use cut/awk to extract only needed columns to temp file
        if (isCompressed) {
          firstLine <- system(paste0(dcatCmd, " ", phenoFile, " | head -n 1"), intern = TRUE)
        } else {
          firstLine <- readLines(phenoFile, n = 1)
        }
        ## detect delimiter and parse column names
        if (grepl("\t", firstLine)) {
          allCols <- strsplit(firstLine, "\t", fixed = TRUE)[[1]]
          cutDelim <- ""
          useAwk <- FALSE
        } else if (grepl(",", firstLine)) {
          allCols <- strsplit(firstLine, ",", fixed = TRUE)[[1]]
          cutDelim <- " -d','"
          useAwk <- FALSE
        } else {
          ## space-delimited: cut cannot handle multiple consecutive spaces,
          ## use awk for correctness (works in BSD awk, gawk, and mawk)
          allCols <- strsplit(trimws(firstLine), "\\s+")[[1]]
          useAwk <- TRUE
        }
        colIndices <- which(allCols %in% checkColList)

        ## extract columns to temp file, then fread from file (fastest path)
        ## use a stable cache name so repeated calls (fallback stages) skip re-extraction
        phenoFiletemp <- file.path(dirname(outputPrefix),
          paste0(".saige_subcols_", phenoCol, "_cache"))
        if (!file.exists(phenoFiletemp)) {
          if (useAwk) {
            awkFields <- paste0("$", colIndices, collapse = ",")
            awkCmd <- paste0("awk 'BEGIN{OFS=\"\\t\"}{print ", awkFields, "}'")
            if (isCompressed) {
              system(paste0("LC_ALL=C ", dcatCmd, " ", phenoFile, " | ", awkCmd, " > ", phenoFiletemp))
            } else {
              system(paste0("LC_ALL=C ", awkCmd, " ", phenoFile, " > ", phenoFiletemp))
            }
          } else {
            colIndicesStr <- paste(colIndices, collapse = ",")
            if (isCompressed) {
              system(paste0("LC_ALL=C ", dcatCmd, " ", phenoFile, " | cut", cutDelim, " -f ", colIndicesStr, " > ", phenoFiletemp))
            } else {
              system(paste0("LC_ALL=C cut", cutDelim, " -f ", colIndicesStr, " ", phenoFile, " > ", phenoFiletemp))
            }
          }
          cat("Phenotype column cache created:", phenoFiletemp, "\n")
        } else {
          cat("Reusing cached phenotype columns:", phenoFiletemp, "\n")
        }

        data <- data.table:::fread(phenoFiletemp,
          header = T, stringsAsFactors = FALSE,
          colClasses = list(character = sampleIDColinphenoFile), data.table = F
        )
      } else {
        ## small file: fread with select is fast enough
        if (isCompressed) {
          data <- data.table:::fread(
            cmd = paste0(dcatCmd, " ", phenoFile),
            header = T, stringsAsFactors = FALSE,
            colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
          )
        } else {
          data <- data.table:::fread(phenoFile,
            header = T, stringsAsFactors = FALSE,
            colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
          )
        }
      }
    } ## end plain text vs HDF5

    cat(sprintf("[TIMING] Phenotype file reading: %.1fs\n", (proc.time() - .t_pheno)[3]))
    .report_mem("after Phenotype file reading")

    if (isRemoveZerosinPheno) {
      data <- data[which(data[, which(colnames(data) == phenoCol)] > 0), ]
      cat("Removing all zeros in the phenotype\n")
      if (nrow(data) == 0) {
        stop("ERROR: no samples are left after removing zeros in the phenotype\n")
      }
    }



    if (SampleIDIncludeFile != "") {
      if (!file.exists(SampleIDIncludeFile)) {
        stop("ERROR! SampleIDIncludeFile ", SampleIDIncludeFile, " does not exsit\n")
      } else {
        sampleIDInclude <- data.table:::fread(SampleIDIncludeFile, header = F, stringsAsFactors = FALSE, colClasses = c("character"), data.table = F)
        sampleIDInclude <- as.vector(sampleIDInclude[!duplicated(sampleIDInclude), ])
        cat(length(sampleIDInclude), " non-duplicated sample IDs were found in SampleIDIncludeFile\n")
        data <- data[which(as.vector(data[, which(colnames(data) == sampleIDColinphenoFile)]) %in% sampleIDInclude), , drop = F]
        cat(nrow(data), " samples in sampleIDInclude have non-missing phenotypes and covariates\n")
      }
    }


    if (length(qCovarCol) > 0) {
      cat(qCovarCol, "are categorical covariates\n")
      if (!all(qCovarCol %in% covarColList)) {
        stop("ERROR! all covariates in qCovarCol must be in covarColList\n")
      } else {
        for (q in qCovarCol) {
          data[, q] <- as.factor(data[, q])
        }
      }
    }

    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")
      if (!all(sampleCovarCol %in% covarColList)) {
        stop("ERROR! all covariates in sampleCovarCol must be in covarColList\n")
      }
    }

    if (FemaleOnly | MaleOnly) {
      if (!sexCol %in% colnames(data)) {
        stop("ERROR! column for sex ", sexCol, " does not exist in the phenoFile \n")
      } else {
        if (FemaleOnly) {
          data <- data[which(data[, which(colnames(data) ==
            sexCol)] == FemaleCode), ]
          if (nrow(data) == 0) {
            stop(
              "ERROR! no samples in the phenotype are coded as ",
              FemaleCode, " in the column ", sexCol,
              "\n"
            )
          }
        } else if (MaleOnly) {
          data <- data[which(data[, which(colnames(data) ==
            sexCol)] == MaleCode), ]
          if (nrow(data) == 0) {
            stop(
              "ERROR! no samples in the phenotype are coded as ",
              MaleCode, " in the column ", sexCol, "\n"
            )
          }
        }
      }
    }


    .t_dataprep <- proc.time()
    print("HERERE2")

    if (length(covarColList) > 0) {
      formula <- paste0(phenoCol, "~", paste0(covarColList,
        collapse = "+"
      ))
      hasCovariate <- TRUE
    } else {
      formula <- paste0(phenoCol, "~ 1")
      hasCovariate <- FALSE
    }

    cat("formula is ", formula, "\n")
    formula.null <- as.formula(formula)
    mmat <- model.matrix(formula.null, data, na.action = NULL)
    mmat <- data.frame(mmat)
    mmat <- cbind(mmat, data[, which(colnames(data) == phenoCol), drop = F])
    colnames(mmat)[ncol(mmat)] <- phenoCol

    coln <- 1
    if (length(offsetCol) > 0) {
      mmat <- cbind(mmat, data[, which(colnames(data) == offsetCol), drop = F])
      colnames(mmat)[ncol(mmat)] <- offsetCol
      coln <- coln + 1
    }

    if (length(varWeightsCol) > 0) {
      mmat <- cbind(mmat, data[, which(colnames(data) == varWeightsCol), drop = F])
      colnames(mmat)[ncol(mmat)] <- varWeightsCol

      coln <- coln + 1
    }



    if (length(covarColList) > 0) {
      if (length(qCovarCol) > 0) {
        covarColList <- colnames(mmat)[2:(ncol(mmat) - coln)]
        formula <- paste0(phenoCol, "~", paste0(covarColList, collapse = "+"))
        formula.null <- as.formula(formula)
      }
    }

    mmat$IID <- data[, which(sampleIDColinphenoFile == colnames(data))]
    if (longlCol != "") {
      mmat$longlVar <- data[, which(longlCol == colnames(data))]
    }

    mmat_nomissing <- mmat[complete.cases(mmat), ]
    mmat_nomissing$IndexPheno <- seq(1, nrow(mmat_nomissing),
      by = 1
    )
    cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")

    if (length(varWeightsCol) > 0) {
      varWeights <- mmat_nomissing[, which(colnames(mmat_nomissing) == varWeightsCol)]
    } else {
      varWeights <- NULL
    }
    if (sparseGRMSampleIDFile != "") {
      # if((useSparseGRMtoFitNULL & !skipVarianceRatioEstimation) | useSparseGRMforVarRatio){
      sampleListwithGenov0 <- data.table:::fread(sparseGRMSampleIDFile,
        header = F, , colClasses = c("character"), data.table = F
      )
      colnames(sampleListwithGenov0) <- c("IIDgeno")
      cat(length(sampleListwithGenov0$IIDgeno), " samples are in the sparse GRM\n")
      mmat_nomissing <- mmat_nomissing[which(mmat_nomissing$IID %in% sampleListwithGenov0$IIDgeno), ]
      cat(nrow(mmat_nomissing), " samples who have non-missing phenotypes are also in the sparse GRM\n")
    }


    if (longlCol == "") {
      if (any(duplicated(mmat_nomissing$IID))) {
        cat("Duplicated sample IDs are detected in the phenotype file. Assuming repeated measurements\n")
      }
    } else {
      cat("Longitudinal variable ", longlCol, " is specified\n")
      if (!any(duplicated(mmat_nomissing$IID))) {
        stop("No duplicated sample IDs are detected in the phenotype file\n")
      }
    }


    if (!is.null(sampleListwithGeno)) {
      dataMerge <- merge(mmat_nomissing, sampleListwithGeno,
        by.x = "IID", by.y = "IIDgeno"
      )
      dataMerge_sort <- dataMerge[with(dataMerge, order(IndexGeno)), ]
      # dataMerge_sort = dataMerge[with(dataMerge, order(IndexPheno)),]
    } else {
      dataMerge_sort <- mmat_nomissing
      dataMerge_sort$IIDgeno <- dataMerge_sort$IID
    }

    print("Test")
    print(head(dataMerge_sort))
    # dataMerge_sort = dataMerge

    if (length(eCovarCol) > 0) {
      cat(eCovarCol, "are environmental covariates\n")
      eMat <- dataMerge_sort[, which(colnames(dataMerge_sort) %in% eCovarCol), drop = F]
      for (em in 1:ncol(eMat)) {
        #eMat[, em] <- (eMat[, em] - mean(eMat[, em])) / (sqrt(ncol(eMat))*sd(eMat[, em]))
        eMat[, em] <- (eMat[, em] - mean(eMat[, em])) / (sd(eMat[, em]))
      }
    eMat = as.matrix(eMat)
    #print("HERE set_EMat")
    #eMat = preprocess_E(eMat)
    set_EMat(eMat) 
    
    eMatcolNames <- checkColList[checkColList %in% eCovarCol]
    eMatIsSample <- eMatcolNames %in% sampleCovarCol
    
    
    ##remove e cov from the covariates list
    #print("Here covarColList 0 ")
    #print(covarColList)
    #covarColList =  covarColList[!(covarColList %in% eCovarCol)]
    #formula <- paste0(phenoCol, "~", paste0(covarColList,
    #  collapse = "+"
    #))
    #formula.null <- as.formula(formula)
    #print("eCovarCol")
    #print(eCovarCol)


    }else{

      eMat = NULL
    }

    #print("Here covarColList")
    #print(covarColList)

    rm(mmat)
    rm(mmat_nomissing)
    gc()
    isSparseGRMIdentity <- FALSE
    if (useGRMtoFitNULL) {
      indicatorGenoSamplesWithPheno <- (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)

      if (length(unique(dataMerge_sort$IIDgeno)) < length(unique(sampleListwithGeno$IIDgeno))) {
        cat(
          length(unique(sampleListwithGeno$IIDgeno)) - length(unique(dataMerge_sort$IIDgeno)),
          " samples in geno file do not have phenotypes\n"
        )
      }
      cat(length(unique(dataMerge_sort$IIDgeno)), " samples will be used for analysis\n")
    } else {
      indicatorGenoSamplesWithPheno <- rep(TRUE, nrow(dataMerge_sort))
    }

    if (any(duplicated(dataMerge_sort$IID))) {
      cat(nrow(dataMerge_sort), " observations will be used for analysis\n")
      # if(longlCol != ""){
      # dataMerge_sort = dataMerge_sort[with(dataMerge_sort, order(IndexGeno, longlVar)),]
      #  dataMerge_sort = dataMerge_sort[with(dataMerge_sort, order(longlVar)),]
      # }
      set_I_mat_inR(dataMerge_sort$IID)
      if (longlCol != "") {
        set_T_mat_inR(dataMerge_sort$IID, dataMerge_sort$longlVar)
      }
    } else {
      if (!useGRMtoFitNULL) {
        # stop("No duplicated IDs are observed in the phenotype file, so GRM must be used to fit the null model. Please set useGRMtoFitNULL=TRUE\n")
        cat("No duplicated IDs are observed in the phenotype file, so the identity matrix will be used as a sparse GRM will be used to fit the null model\n")
        isSparseGRMIdentity <- TRUE
        useSparseGRMtoFitNULL <- TRUE
        useGRMtoFitNULL <- TRUE
      }
    }
    set_useGRMtoFitNULL(useGRMtoFitNULL)
  }



  print("Test3")
  print(head(dataMerge_sort))

  if (traitType == "quantitative" & invNormalize) {
    cat(
      "Perform the inverse nomalization for ", phenoCol,
      "\n"
    )
    invPheno <- qnorm((rank(dataMerge_sort[, which(colnames(dataMerge_sort) ==
      phenoCol)], na.last = "keep") - 0.5) / sum(!is.na(dataMerge_sort[
      ,
      which(colnames(dataMerge_sort) == phenoCol)
    ])))
    dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] <- invPheno
  }

  #if(length(eCovarCol) > 0){
  #      dataMerge_sort <- dataMerge_sort[, !(names(dataMerge_sort) %in%
  #	      eCovarCol)]
  #}
  print("Test4")
  print(head(dataMerge_sort))
  if (traitType == "binary" & (length(covarColList) > 0)) {
    out_checksep <- checkPerfectSep(formula.null,
      data = dataMerge_sort,
      minCovariateCount
    )
    covarColList <- covarColList[!(covarColList %in% out_checksep)]
    formula <- paste0(phenoCol, "~", paste0(covarColList,
      collapse = "+"
    ))
    formula.null <- as.formula(formula)
    if (length(covarColList) == 1) {
      hasCovariate <- FALSE
    } else {
      hasCovariate <- TRUE
    }
    dataMerge_sort <- dataMerge_sort[, !(names(dataMerge_sort) %in%
      out_checksep)]
  }
  if (!hasCovariate) {
    print("No covariate is includes so isCovariateOffset = FALSE")
    isCovariateOffset <- FALSE
  }

  if (isCovariateTransform & hasCovariate) {
    cat("qr transformation has been performed on covariates\n")
    out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
    formulaNewList <- c(phenoCol, " ~ ", out.transform$Param.transform$X_name[1])
    if (length(out.transform$Param.transform$X_name) > 1) {
      for (i in c(2:length(out.transform$Param.transform$X_name))) {
        formulaNewList <- c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
      }
    }
    formulaNewList <- paste0(formulaNewList, collapse = "")
    formulaNewList <- paste0(formulaNewList, "-1")
    formula.new <- as.formula(paste0(formulaNewList, collapse = ""))
    data.new <- data.frame(cbind(out.transform$Y, out.transform$X1))
    colnames(data.new) <- c(phenoCol, out.transform$Param.transform$X_name)
    cat("colnames(data.new) is ", colnames(data.new), "\n")
    cat(
      "out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr),
      "\n"
    )

    if (length(offsetCol) > 0) {
      data.new <- cbind(data.new, dataMerge_sort[,which(colnames(dataMerge_sort) == offsetCol)])
      colnames(data.new)[ncol(data.new)] = offsetCol
    }
  } else {
    formula.new <- formula.null
    data.new <- dataMerge_sort
    out.transform <- NULL
  }

  if (traitType == "binary") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = binomial, weights = varWeights
      )
    } else {
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = binomial, weights = varWeights
      )
    }
  } else if (traitType == "quantitative") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = gaussian(link = "identity"), weights = varWeights
      )
    } else {
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = gaussian(link = "identity"), weights = varWeights
      )
    }
  } else if (traitType == "count") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = "poisson", weights = varWeights
      )
    } else {
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = "poisson", weights = varWeights
      )
    }
  } else if (traitType == "count_nb") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = NegBin(), weights = varWeights
      )
    } else {
      print(head(data.new))
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = NegBin(), weights = varWeights
      )
    }
  }
  mmat <- model.matrix(formula.new, data = data.new, na.action = NULL)

  if (isCovariateOffset) {
    covoffset <- mmat[, -1, drop = F] %*% modwitcov$coefficients[-1]
    print("isCovariateOffset=TRUE, so fixed effects coefficnets won't be estimated.")
    formula.new.withCov <- formula.new
    formula_nocov <- paste0(phenoCol, "~ 1")
    formula.new <- as.formula(formula_nocov)
    hasCovariate <- FALSE
  } else {
    covoffset <- rep(0, nrow(data.new))
  }

  data.new$covoffset <- covoffset


  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    if (!isSparseGRMIdentity) {
      getsubGRM_orig(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID)
    } else {
      print(length(dataMerge_sort$IIDgeno))
      sparseGRM <- Matrix:::sparseMatrix(i = as.vector(1:nrow(data.new)), j = as.vector(1:nrow(data.new)), x = rep(1, nrow(data.new)), symmetric = TRUE)
      setupSparseGRM_new(sparseGRM)
    }
    # getsubGRM(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID, dataMerge_sort$longlVar)
    # m4 = gen_sp_v2(sparseGRMtest)
    # cat("Setting up sparse GRM using ", sparseGRMFile, " and ", sparseGRMSampleIDFile, "\n")
    # cat("Dimension of the sparse GRM is ", dim(m4), "\n")
    # A = summary(m4)
    # locationMatinR = rbind(A$i - 1, A$j - 1)
    # valueVecinR = A$x
    # setupSparseGRM(dim(m4)[1], locationMatinR, valueVecinR)
    # setupSparseGRM_new(sparseGRMtest)
    # rm(sparseGRMtest)
    gc()
  }

  # allow for multiple variance components
  # set_Vmat_vec(VmatFilelist, VmatSampleFilelist, dataMerge_sort$IID, dataMerge_sort$longlVar)
  set_Vmat_vec_orig(VmatFilelist, VmatSampleFilelist, dataMerge_sort$IID)

  numofV <- get_numofV()

  print(dataMerge_sort$IID[1:200])
  print(any(duplicated(dataMerge_sort$IID)))

  cat("numofV ", numofV, "\n")

if(is.null(eMat)){  
  if (any(duplicated(dataMerge_sort$IID))) {
    print("HERE")
    if (longlCol == "") {
      print("HERE1")
      print(useGRMtoFitNULL)
      if (useGRMtoFitNULL) {
        print("HERE2")
        num_Kmat <- numofV + 3
        cat("num_Kmat ", num_Kmat, "\n")
      } else {
        num_Kmat <- numofV + 2
      }
      # k = num_Kmat + 3
    } else {
      if (useGRMtoFitNULL) {
        num_Kmat <- 7 + numofV * 3
      } else {
        num_Kmat <- 4 + numofV * 3
      }
      # k = num_Kmat + 2
    }
  } else {
    if (useGRMtoFitNULL) {
      num_Kmat <- numofV + 2
    } else {
      num_Kmat <- numofV + 1
    }
    # k = 2
  }
}else{ #if(nrow(eMat) == 0){
  if (any(duplicated(dataMerge_sort$IID))) {
    num_Kmat <- numofV + 3
  }else{
    if (useGRMtoFitNULL) {
      num_Kmat <- numofV + 3
    }else{
      num_Kmat <- numofV + 2
    }
  }
}


  k <- num_Kmat

  
  set_num_Kmat(num_Kmat)
  cat("num_Kmat ", num_Kmat, "\n")


  if(!is.null(eMat)){
        print("dim(eMat)")
  	print(dim(eMat))
  #	k = 3
  }


  if (longlCol != "") {
    covarianceIdxMat <- set_covarianceidx_Mat()
  } else {
    covarianceIdxMat <- NULL
  }

  if (!skipVarianceRatioEstimation) {
    isVarianceRatioinGeno <- TRUE
    if (isCateVarianceRatio) {
      minMAC_varRatio <- min(cateVarRatioMinMACVecExclude)
      maxMAC_varRatio <- max(cateVarRatioMaxMACVecInclude)
      cat("Categorical variance ratios will be estimated. Please make sure there are at least 200 markers in each MAC category.\n")
    } else {
      minMAC_varRatio <- 20
      maxMAC_varRatio <- -1 # will randomly select markers from the plink file and leave them out when constructing GRM
    }
    setminMAC_VarianceRatio(minMAC_varRatio, maxMAC_varRatio, isVarianceRatioinGeno)
  }

  # set up parameters
  if (minMAFforGRM > 0) {
    cat(
      "Markers in the Plink file with MAF < ", minMAFforGRM,
      " will be removed before constructing GRM\n"
    )
  }
  if (maxMissingRateforGRM > 0) {
    cat("Markers in the Plink file with missing rate > ", maxMissingRateforGRM, " will be removed before constructing GRM\n")
  }

  setminMAFforGRM(minMAFforGRM)
  setmaxMissingRateforGRM(maxMissingRateforGRM)

  cat(sprintf("[TIMING] Data prep (merge, model matrix, eMat, GRM): %.1fs\n", (proc.time() - .t_dataprep)[3]))
  .report_mem("after Data prep")

  .t_glm <- proc.time()
  if (traitType == "binary") {
    cat(phenoCol, " is a binary trait\n")
    uniqPheno <- sort(unique(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)]))
    if (uniqPheno[1] != 0 | uniqPheno[2] != 1) {
      stop("ERROR! phenotype value needs to be 0 or 1 \n")
    }
    print("formula.new")
    print(formula.new)
    print("head(data.new)")
    print(head(data.new))
    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = binomial, weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = binomial, weights = varWeights)
      }
      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = binomial, weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = binomial, weights = varWeights
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = binomial, weights = varWeights)
      }
    }
  } else if (traitType == "quantitative") {
    cat(phenoCol, " is a quantitative trait\n")
    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = gaussian(link = "identity"), weights = varWeights)
      }
      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = gaussian(link = "identity"), weights = varWeights
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = gaussian(link = "identity"), weights = varWeights)
      }
    }
  } else if (traitType == "count") {
    cat(phenoCol, " is a count trait\n")
    # print("before remove zeros")
    # print(dim(dataMerge_sort))
    # if(isRemoveZerosinPheno){
    #    dataMerge_sort = dataMerge_sort[which(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] > 0),]
    #    cat("Removing all zeros in the phenotype\n")
    #    if(nrow(dataMerge_sort) == 0){
    #        stop("ERROR: no samples are left after removing zeros in the phenotype\n")

    #    }
    # }
    # print("after remove zeros")
    # print(dim(dataMerge_sort))
    miny <- min(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)])
    if (miny < 0) {
      stop("ERROR! phenotype value needs to be non-negative \n")
    }

    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = "poisson", weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = "poisson", weights = varWeights)
      }
      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = "poisson", weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = "poisson", weights = varWeights
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = "poisson", weights = varWeights)
      }
    }
  } else if (traitType == "count_nb") {
    cat(phenoCol, " is a count_nb trait\n")
    if (isRemoveZerosinPheno) {
      dataMerge_sort <- dataMerge_sort[which(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] > 0), ]
      cat("Removing all zeros in the phenotype\n")
      if (nrow(dataMerge_sort) == 0) {
        stop("ERROR: no samples are left after removing zeros in the phenotype\n")
      }
    }
    miny <- min(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)])
    if (miny < 0) {
      stop("ERROR! phenotype value needs to be non-negative \n")
    }
    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = NegBin(), weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = NegBin(), weights = varWeights)
      }

      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = NegBin(), weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = NegBin()
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = NegBin(), weights = varWeights)
      }
    }
  }


  cat("glm:\n")
  print(fit0)
  obj.noK <- NULL


  # if(length(fit0$y) > 200000){
  #  isStoreSigma = FALSE
  # }else{
  #  isStoreSigma = TRUE
  # }
  print("isStoreSigma")
  print(isStoreSigma)
  #set_store_sigma(isStoreSigma)

  cat(sprintf("[TIMING] Initial GLM fit: %.1fs\n", (proc.time() - .t_glm)[3]))
  .report_mem("after Initial GLM fit")

  if (!skipModelFitting) {
    # setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
    cat("Start fitting the NULL GLMM\n")
    t_begin <- proc.time()
    print(t_begin)


    tau <- rep(0, k)
    fixtau <- rep(0, k)
    # tauInit = tau

    set_isSparseGRM(useSparseGRMtoFitNULL)
    set_useGRMtoFitNULL(useGRMtoFitNULL)

    if (traitType != "count_nb") {
    if(length(eCovarCol) == 0){
      system.time(modglmm <- glmmkin.ai_PCG_Rcpp_multiV(bedFile, bimFile, famFile, Xorig, isCovariateOffset,
        fit0,
        tau = tau, fixtau = fixtau, maxiter = maxiter,
        tol = tol, verbose = TRUE, nrun = nrun, tolPCG = tolPCG,
        maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno,
        obj.noK = obj.noK, out.transform = out.transform,
        tauInit = tauInit, memoryChunk = memoryChunk,
        LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec,
        chromosomeEndIndexVec = chromosomeEndIndexVec,
        traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        isLowMemLOCO = isLowMemLOCO, covarianceIdxMat = covarianceIdxMat, isStoreSigma = isStoreSigma, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL, useGRMtoFitNULL = useGRMtoFitNULL, isSparseGRMIdentity = isSparseGRMIdentity
      ))
    }else{ #if(length(eCovarCol) == 0){

     print("glmmkin.ai_PCG_Rcpp_multiV_eMat")
     print(fit0)

      system.time(modglmm <- glmmkin.ai_PCG_Rcpp_multiV_eMat(bedFile, bimFile, famFile, Xorig, isCovariateOffset,
        fit0,
        tau = tau, fixtau = fixtau, maxiter = maxiter,
        tol = tol, verbose = TRUE, nrun = nrun, tolPCG = tolPCG,
        maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno,
        obj.noK = obj.noK, out.transform = out.transform,
        tauInit = tauInit, memoryChunk = memoryChunk,
        LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec,
        chromosomeEndIndexVec = chromosomeEndIndexVec,
        traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        isLowMemLOCO = isLowMemLOCO, covarianceIdxMat = covarianceIdxMat, isStoreSigma = isStoreSigma, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL, useGRMtoFitNULL = useGRMtoFitNULL, isSparseGRMIdentity = isSparseGRMIdentity
      ))


    }


      modglmm$obj.glm.null$model <- data.frame(modglmm$obj.glm.null$model)
      cat(sprintf("[TIMING] AIREML null model fitting: %.1fs\n", (proc.time() - t_begin)[3]))
      .report_mem("after AIREML null model fitting")
    } else {
      system.time(modglmm <- glmmkin.ai_PCG_Rcpp_multiV_NB(bedFile, bimFile, famFile, Xorig, isCovariateOffset,
        fit0,
        tau = tau, fixtau = fixtau, maxiter = maxiter,
        tol = tol, verbose = TRUE, nrun = nrun, tolPCG = tolPCG,
        maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno,
        obj.noK = obj.noK, out.transform = out.transform,
        tauInit = tauInit, memoryChunk = memoryChunk,
        LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec,
        chromosomeEndIndexVec = chromosomeEndIndexVec,
        traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        isLowMemLOCO = isLowMemLOCO, covarianceIdxMat = covarianceIdxMat, isStoreSigma = isStoreSigma, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL, useGRMtoFitNULL = useGRMtoFitNULL, isSparseGRMIdentity = isSparseGRMIdentity
      ))

      data.new$y <- modglmm$y
      varWeights <- modglmm$varWeights
      if (!isCovariateOffset) {
        if (length(offsetCol) == 0) {
          fit0 <- glm(formula.new, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
        } else {
          offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
          fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = gaussian(link = "identity"), weights = varWeights)
        }
        Xorig <- NULL
      } else {
        fit0orig <- glm(formula.new.withCov, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
        Xorig <- model.matrix(fit0orig)
        rm(fit0orig)
        gc()
        if (length(offsetCol) == 0) {
          fit0 <- glm(formula.new,
            data = data.new, offset = covoffset,
            family = gaussian(link = "identity"), weights = varWeights
          )
        } else {
          offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
          fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = gaussian(link = "identity"), weights = varWeights)
        }
      }

      modglmm$obj.glm.null <- fit0
    }


    # spSigma_final = getSparseSigma_new()
    # modglmm$spSigma = spSigma_final
    # rm(spSigma_final)
    if (traitType != "count_nb") {
      for (x in names(modglmm$obj.glm.null)) {
        attr(modglmm$obj.glm.null[[x]], ".Environment") <- c()
      }
    }
    # modglmm$offset = covoffset
    if (isCovariateOffset) {
      modglmm$offset <- covoffset
    } else {
      if (hasCovariate) {
        data.new.X <- model.matrix(fit0)[, -1, drop = F]
        print(head(data.new))
        print(head(data.new.X))
        print(head(modglmm$coefficients[-1]))
        modglmm$offset <- data.new.X %*% (as.vector(modglmm$coefficients[-1]))
        if (LOCO) {
          for (j in 1:22) {
            if (modglmm$LOCOResult[[j]]$isLOCO) {
              modglmm$LOCOResult[[j]]$offset <- data.new.X %*% (as.vector(modglmm$LOCOResult[[j]]$coefficients[-1]))
            }
          }
        }
      } else {
        modglmm$offset <- covoffset
      }
    }

####to simplify

    if (length(eCovarCol) > 0) {
      cat(eCovarCol, "are environmental covariates\n")
      modglmm$eMat = eMat
      eMatcolNames = checkColList[checkColList %in% eCovarCol]
      eMatIsSample = eMatcolNames %in% sampleCovarCol
      modglmm$eMatcolNames = eMatcolNames
      modglmm$eMatIsSample = eMatIsSample
      
      
      #modglmm$eMat <- data.new[, which(colnames(data.new) %in% eCovarCol), drop = F]
      #for (em in 1:ncol(modglmm$eMat)) {
      #  modglmm$eMat[, em] <- (modglmm$eMat[, em] - mean(modglmm$eMat[, em])) / (sd(modglmm$eMat[, em]))
      #}
    }
#####



    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")
      modglmm$sampleXMat <- modglmm$X[, which(colnames(modglmm$X) %in% sampleCovarCol), drop = F]
      modglmm$sampleXMat <- cbind(modglmm$X[, 1], modglmm$sampleXMat)

      uniqsampleind <- which(!duplicated(modglmm$sampleID))
      modglmm$sampleXMat <- modglmm$sampleXMat[uniqsampleind, ]
    }








    .t_post <- proc.time()
    # if((skipVarianceRatioEstimation & useSparseGRMtoFitNULL)){
    family <- fit0$family
    eta <- modglmm$linear.predictors
    mu <- modglmm$fitted.values
    mu.eta <- family$mu.eta(eta)
    sqrtW <- mu.eta / sqrt(family$variance(mu))
    W <- sqrtW^2
    W <- W * modglmm$varWeights
    tauVecNew <- modglmm$theta
    .t_pp_sigmaiX <- proc.time()
    Sigma_iX <- getSigma_X_multiV(W, tauVecNew, modglmm$X, maxiterPCG, tolPCG, LOCO = FALSE)
    if (verbose) cat(sprintf("[TIMING:PP] getSigma_X_multiV: %.2fs\n", (proc.time() - .t_pp_sigmaiX)[3]))
    if (!isShrinkModelOutput) {
      Sigma_iXXSigma_iX <- Sigma_iX %*% (solve(t(modglmm$X) %*% Sigma_iX))
      modglmm$Sigma_iXXSigma_iX <- Sigma_iXXSigma_iX
    }
    # }

    modglmm$useSparseGRMforVarRatio <- useSparseGRMforVarRatio


    # if(any(duplicated(modglmm$sampleID))){
    if (useGRMtoFitNULL) {
      modglmm$tauVal_sp <- modglmm$theta[3]
    } else {
      modglmm$tauVal_sp <- modglmm$theta[2]
    }
    # }

    # if(FALSE){
    # if(length(fit0$y) <= 10000){

    cat("isStoreSigma is ", isStoreSigma, "\n")
    .t_pp_spsigma <- proc.time()
    if (isStoreSigma) {
      ## Try cached SMW Case 2 path: replaces n_donors full Case 2 setups with
      ## one setup + n_donors cheap applies. Falls back to gettI_Sigma_I_multiV
      ## if preconditions don't hold.
      pp_smw_cached <- FALSE
      if (!is.null(modglmm$eMat) && length(tauVecNew) >= 3 && tauVecNew[3] != 0) {
        .pp_dec <- .smw_cache_decision(modglmm$eMat, "post-processing")
        if (.pp_dec$allowed) {
          pp_smw_cached <- isTRUE(prepareSigmaInvSMW_multiV(W, tauVecNew))
        }
      }
      if (pp_smw_cached) {
        modglmm$spSigma <- gettI_Sigma_I_multiV_cached()
        clearSigmaInvSMWcache()
      } else {
        modglmm$spSigma <- gettI_Sigma_I_multiV(W, tauVecNew, maxiterPCG, tolPCG, LOCO = FALSE)
      }
      if (verbose) cat(sprintf("[TIMING:PP] gettI_Sigma_I_multiV: %.2fs (cached=%s)\n",
                  (proc.time() - .t_pp_spsigma)[3], pp_smw_cached))
      # }else{
      # 	gen_sp_Sigma_multiV(W, tauVecNew)
      # 	spSigma = get_sp_Sigma_to_R()
      # 	SigmaMat_sp = chol2inv(chol(spSigma))
      # 	modglmm$spSigma = SigmaMat_sp
      # }
    }
    # }
    # save(modglmm, file = modelOut)
    tau <- modglmm$theta
    alpha0 <- modglmm$coefficients

    #Sigma_iX
 if(FALSE){   
    if (!is.null(modglmm$eMat)) {
      etileMat = NULL
      P2Mat_E = NULL
      Scorevec_E = NULL
      VarMat_E = matrix(0, nrow=ncol(eMat), ncol=ncol(eMat))
      #if(isShrinkModelOutput){Sigma_iXXSigma_iX <- Sigma_iX %*% (solve(t(modglmm$X) %*% Sigma_iX))}
	for (ne in 1:ncol(eMat)) {
		evec = eMat[, ne]
		evec_tilde <- evec - modglmm$obj.noK$XXVX_inv %*% (modglmm$obj.noK$XV %*% evec)
		etileMat = cbind(etileMat, evec_tilde)
	}



      for (ne in 1:ncol(eMat)) {
      	evec_tilde = etileMat[, ne]
	S_E <- innerProduct(evec_tilde, modglmm$residuals * modglmm$varWeights)
	Scorevec_E = c(Scorevec_E, S_E)
	Sigma_iE <- getSigma_G_multiV(W, tauVecNew, evec_tilde, maxiterPCG, tolPCG, LOCO = FALSE)
	for(ne2 in 1:ncol(eMat)){
		evec_tilde2 = etileMat[, ne2]
		varE = t(evec_tilde2) %*% Sigma_iE  - t(evec_tilde2)%*% Sigma_iX %*% (solve(t(modglmm$X) %*% Sigma_iX)) %*% t(modglmm$X) %*% Sigma_iE	 
		VarMat_E[ne, ne2] = varE
	}
      }
      modglmm$etileMat = etileMat
      modglmm$Scorevec_E = Scorevec_E
      modglmm$VarMat_E = VarMat_E
    }  

}# if(FALSE){


    if (!is.null(out.transform) & is.null(fit0$offset)) {
      coef.alpha <- Covariate_Transform_Back(alpha0, out.transform$Param.transform)
      modglmm$coefficients <- coef.alpha
    }


    if (LOCO & isLowMemLOCO) {
      modglmm$LOCOResult <- NULL
      modglmm$LOCO <- FALSE
      chromosomeStartIndexVec <- modglmm$chromosomeStartIndexVec
      chromosomeEndIndexVec <- modglmm$chromosomeEndIndexVec
      modglmm$chromosomeStartIndexVec <- NULL
      modglmm$chromosomeEndIndexVec <- NULL
      modelOut <- paste(c(outputPrefix, "_noLOCO.rda"), collapse = "")
      save(modglmm, file = modelOut)
      modglmm$LOCO <- TRUE
      modglmm$Y <- NULL
      eta0 <- modglmm$linear.predictors
      modglmm$linear.predictors <- NULL
      modglmm$coefficients <- NULL
      modglmm$cov <- NULL
      modglmm$fitted.values <- NULL
      modglmm$residuals <- NULL
      modglmm$obj.noK <- NULL
      offset0 <- modglmm$offset
      modglmm$offset <- NULL
      y <- fit0$y
      gc()
      # save(modglmm, file = modelOut)
      set_Diagof_StdGeno_LOCO()
      modglmm$LOCOResult <- list()
      for (j in 1:22) {
        startIndex <- chromosomeStartIndexVec[j]
        endIndex <- chromosomeEndIndexVec[j]
        if (!is.na(startIndex) && !is.na(endIndex)) {
          cat("leave chromosome ", j, " out\n")
          setStartEndIndex(startIndex, endIndex, j - 1)

          re.coef_LOCO <- Get_Coef_multiV(y, X = model.matrix(fit0), tau, family = fit0$family, alpha = alpha0, eta = eta0, offset = offset0, verbose = TRUE, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = TRUE, var_weights = var_weights)
          cov <- re.coef_LOCO$cov
          alpha <- re.coef_LOCO$alpha
          eta <- re.coef_LOCO$eta
          Y <- re.coef_LOCO$Y
          mu <- re.coef_LOCO$mu
          if (family$family == "binomial") {
            mu2 <- mu * (1 - mu)
          } else if (family$family == "poisson") {
            mu2 <- mu
          } else if (family$family == "gaussian") {
            mu2 <- rep((1 / (tau[1])), length(res))
          } else if (traitType == "count_nb") {
            # mu2 = fit0$family$variance(mu)
            mu2 <- rep((1 / (tau[1])), length(res))
          }
          res <- y - mu

          if (!is.null(out.transform) & is.null(fit0$offset)) {
            coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
          } else {
            coef.alpha <- alpha
          }

          mu2_rescaled <- mu2 * var_weights
          mu_rescaled <- mu * var_weights


          if (!isCovariateOffset) {
            obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
          } else {
            obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
          }
          modglmm$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK)
          if (!isCovariateOffset & hasCovariate) {
            data.new.X <- model.matrix(fit0)[, -1, drop = F]
            modglmm$LOCOResult[[j]]$offset <- data.new.X %*% (as.vector(modglmm$LOCOResult[[j]]$coefficients[-1]))
          }
          modelOutbychr <- paste(c(outputPrefix, "_chr", j, ".rda"), collapse = "")
          if (j != 22) {
            for (j1 in (j + 1):22) {
              modglmm$LOCOResult[[j1]] <- list(NULL)
            }
          }
          save(modglmm, file = modelOutbychr)
          modglmm$LOCOResult[[j]] <- list(NULL)
          gc()
        } else {
          modglmm$LOCOResult[[j]] <- list(isLOCO = FALSE)
        }
      }
      gc()
      # modelOut_nonauto = paste(c(outputPrefix,"_noLOCO.rda"), collapse="")
    } else {
      # b = as.numeric(factor(dataMerge_sort$IID, levels =  unique(dataMerge_sort$IID)))
      # I_mat = Matrix::sparseMatrix(i = 1:length(b), j = b, x = rep(1, length(b)))
      # I_mat = 1.0 * I_mat
      # modglmm$I_longl_mat = I_mat
      # modglmm$I_longl_vec = b - 1
      # modglmm$T_longl_mat = I_mat * (dataMerge_sort$longlVar)
      modglmm$T_longl_vec <- dataMerge_sort$longlVar
      .t_pp_save <- proc.time()
      save(modglmm, file = modelOut)
      if (verbose) cat(sprintf("[TIMING:PP] save modglmm: %.2fs\n", (proc.time() - .t_pp_save)[3]))
    }

  if(sum(modglmm$theta[2:length(modglmm$theta)]) < 0 || sum(modglmm$theta[2:length(modglmm$theta)]) > 10 || (!modglmm$isCovariateOffset && any(modglmm$theta[2:length(modglmm$theta)] == 0))){

	cat(modglmm$theta)
	stop("Tau estimates out of bound: possible model divergence")
  } 
    


    t_end <- proc.time()
    print(t_end)
    cat("t_end - t_begin, fitting the NULL model took\n")
    print(t_end - t_begin)


    # if(bedFile != "" & !useGRMtoFitNULL){
    if (bedFile != "") {
      subSampleInGeno <- dataMerge_sort$IndexGeno
      if (is.null(dataMerge_sort$IndexGeno)) {
        subSampleInGeno <- dataMerge_sort$IndexPheno
      }

      print(subSampleInGeno[1:1000])
      print(head(dataMerge_sort))
      print("HEREHRE")

      subSampleInGeno_unique <- subSampleInGeno[!duplicated(subSampleInGeno)]

      .t_pp_setgeno <- proc.time()
      # setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
      setgeno(bedFile, bimFile, famFile, subSampleInGeno_unique, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
      if (verbose) cat(sprintf("[TIMING:PP] setgeno: %.2fs\n", (proc.time() - .t_pp_setgeno)[3]))
    }
  } else {
    cat("Skip fitting the NULL GLMM\n")
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
    load(modelOut)
    if (is.null(modglmm$LOCO)) {
      modglmm$LOCO <- FALSE
    }

    # need check
    subSampleInGeno <- dataMerge_sort$IndexGeno
    if (is.null(dataMerge_sort$IndexGeno)) {
      subSampleInGeno <- dataMerge_sort$IndexPheno
    }

    print(subSampleInGeno[1:1000])
    print(head(dataMerge_sort))
    print("HEREHRE")

    subSampleInGeno_unique <- subSampleInGeno[!duplicated(subSampleInGeno)]

    # setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    setgeno(bedFile, bimFile, famFile, subSampleInGeno_unique, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)


    # setgeno(bedFile, bimFile, famFile, dataMerge_sort$IndexGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    tau <- modglmm$theta

    if (any(duplicated(modglmm$sampleID))) {
      set_I_mat_inR(modglmm$sampleID)
      # if(longlCol != ""){
      #        set_T_mat_inR(dataMerge_sort$IID, dataMerge_sort$longlVar)
      # }
    }
    set_dup_sample_index(as.numeric(factor(modglmm$sampleID, levels = unique(modglmm$sampleID))))
    # setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
  }

  cat(sprintf("[TIMING] Post-processing (Sigma_iX, spSigma, save): %.1fs\n", (proc.time() - .t_post)[3]))
  .report_mem("after Post-processing")

  .t_varratio <- proc.time()
  if (!skipVarianceRatioEstimation) {
    if (LOCO) {
      MsubIndVec <- getQCdMarkerIndex()
      print(length(MsubIndVec))
      chrVec <- data.table:::fread(bimFile, header = F)[, 1]
      print(length(chrVec))
      chrVec <- chrVec[which(MsubIndVec == TRUE)]
      updatechrList <- updateChrStartEndIndexVec(chrVec)
      LOCO <- updatechrList$LOCO
      chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
      chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
      set_Diagof_StdGeno_LOCO()
    }
    cat("Start estimating variance ratios\n")
    load(modelOut)
    use_sandwich = extractVarianceRatio_multiV(
      obj.glmm.null = modglmm,
      obj.glm.null = fit0, maxiterPCG = maxiterPCG,
      tolPCG = tolPCG, numMarkers = numMarkersForVarRatio, varRatioOutFile = varRatioFile,
      ratioCVcutoff = ratioCVcutoff, testOut = SPAGMMATOut,
      bedFile = bedFile, bimFile = bimFile, famFile = famFile, chromosomeStartIndexVec = chromosomeStartIndexVec,
      chromosomeEndIndexVec = chromosomeEndIndexVec,
      isCateVarianceRatio = isCateVarianceRatio, cateVarRatioIndexVec = cateVarRatioIndexVec,
      useSparseGRMforVarRatio = useSparseGRMforVarRatio, sparseGRMFile = sparseGRMFile,
      sparseGRMSampleIDFile = sparseGRMSampleIDFile,
      numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
      relatednessCutoff = relatednessCutoff, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL,
      nThreads = nThreads, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
      minMAFforGRM = minMAFforGRM, isDiagofKinSetAsOne = isDiagofKinSetAsOne,
      includeNonautoMarkersforVarRatio = includeNonautoMarkersforVarRatio, isStoreSigma = isStoreSigma, useGRMtoFitNULL = useGRMtoFitNULL,
      smwCacheMemLimitMB = smwCacheMemLimitMB,
      verbose = verbose
    )
  } else {
    use_sandwich = NULL
    cat("Skip estimating variance ratios\n")
  }
  closeGenoFile_plink()
  # clean up saved model (as in ReadModel)
  if (isShrinkModelOutput) {
    load(modelOut)
    modglmm$use_sandwich = use_sandwich
    modglmm$Y <- NULL
    modglmm$linear.predictors <- NULL
    modglmm$coefficients <- NULL
    modglmm$cov <- NULL
    if (sum(duplicated(modglmm$sampleID)) > 0) {
      modglmm$obj.noK$Sigma_iXXSigma_iX <- matrix(1)
      modglmm$X <- NULL
      if (is.null(modglmm$eMat)) {
        modglmm$obj.noK$XV <- NULL
        modglmm$obj.noK$XVX <- NULL
        modglmm$obj.noK$XXVX_inv <- NULL
        modglmm$obj.noK$XVX_inv <- NULL
        modglmm$obj.noK$XVX_inv_XV <- NULL
      }
    }
  } else {
    modglmm$use_sandwich = use_sandwich
  }
  save(modglmm, file = modelOut)

  cat(sprintf("[TIMING] Variance ratio estimation + final save: %.1fs\n", (proc.time() - .t_varratio)[3]))
  .report_mem("after Variance ratio estimation")
  cat(sprintf("[TIMING] ===== Total fitNULLGLMM_multiV: %.1fs =====\n", (proc.time() - .t0_total)[3]))
}



extractVarianceRatio_multiV <- function(obj.glmm.null,
                                        obj.glm.null,
                                        maxiterPCG = 500,
                                        tolPCG = 0.01,
                                        numMarkers,
                                        varRatioOutFile,
                                        ratioCVcutoff,
                                        testOut,
                                        bedFile,
                                        bimFile,
                                        famFile,
                                        chromosomeStartIndexVec,
                                        chromosomeEndIndexVec,
                                        isCateVarianceRatio,
                                        cateVarRatioIndexVec,
                                        useSparseGRMforVarRatio,
                                        sparseGRMFile,
                                        sparseGRMSampleIDFile,
                                        numRandomMarkerforSparseKin,
                                        relatednessCutoff,
                                        useSparseGRMtoFitNULL,
                                        nThreads,
                                        cateVarRatioMinMACVecExclude,
                                        cateVarRatioMaxMACVecInclude,
                                        minMAFforGRM,
                                        isDiagofKinSetAsOne,
                                        includeNonautoMarkersforVarRatio,
                                        isStoreSigma = FALSE,
                                        useGRMtoFitNULL = TRUE,
                                        smwCacheMemLimitMB = NULL,
                                        verbose = FALSE) {
  obj.noK <- obj.glmm.null$obj.noK
  if (file.exists(testOut)) {
    file.remove(testOut)
  }
  bimPlink <- data.frame(data.table:::fread(bimFile, header = F))
  if (sum(sapply(bimPlink[, 1], is.numeric)) != nrow(bimPlink)) {
    stop("ERROR: chromosome column in plink bim file is no numeric!\n")
  }

  family <- obj.glm.null$family
  print(family)
  eta <- obj.glmm.null$linear.predictors
  mu <- obj.glmm.null$fitted.values
  mu.eta <- family$mu.eta(eta)

  # var_weights = weights(obj.glm.null)
  var_weights <- obj.glmm.null$varWeights
  # if(!is.null(var_weights)){
  sqrtW <- mu.eta / sqrt(family$variance(mu))
  # }else{
  #  sqrtW = mu.eta/sqrt(family$variance(mu))
  # }


  print("mu[1:20]")
  print(mu[1:20])
  W <- sqrtW^2 ## (mu*(1-mu) for binary)
  W <- W * var_weights
  print("mu[1:20]")
  print(W[1:20])
  tauVecNew <- obj.glmm.null$theta

  # isStoreSigma=FALSE
  # if(isStoreSigma){
  #       gen_sp_Sigma_multiV(W, tauVecNew)
  # }
  X <- obj.glmm.null$X


  set_isSparseGRM(useSparseGRMtoFitNULL)
  set_useGRMtoFitNULL(useGRMtoFitNULL)

  # if(!useGRMtoFitNULL){
  # if(useSparseGRMforVarRatio){
  # }
  # useSparseGRMforVarRatio = FALSE
  # }

  ## Try to enable cached SMW (Case 2) solver for ALL Sigma^-1 calls inside this
  ## function (Sigma_iX_noLOCO + per-marker Sigma_iG / Sigma_iGE).
  ## Falls back to getSigma_X_multiV / getSigma_G_multiV (which still route to
  ## their own Case 2 / 3) if preconditions don't hold.
  smw_cached <- FALSE
  if (!is.null(obj.glmm.null$eMat) && length(tauVecNew) >= 3 && tauVecNew[3] != 0) {
    .vr_N <- nrow(obj.glmm.null$eMat); .vr_k <- ncol(obj.glmm.null$eMat)
    .vr_D <- length(unique(obj.glmm.null$sampleID))
    .vr_estMB <- 4 * (.vr_N * (2 + 2 * .vr_k) + .vr_D * .vr_k^2) / (1024 * 1024)
    .vr_over  <- !is.null(smwCacheMemLimitMB) && .vr_estMB > smwCacheMemLimitMB
    cat(sprintf("[SMW] variance-ratio: estimated cache size = %.0f MB (N=%d, k=%d, D=%d)%s\n",
                .vr_estMB, .vr_N, .vr_k, .vr_D,
                if (.vr_over) sprintf(" -> EXCEEDS smwCacheMemLimitMB=%.0f, skipping",
                                      smwCacheMemLimitMB) else ""))
    if (!.vr_over) {
      smw_cached <- isTRUE(prepareSigmaInvSMW_multiV(W, tauVecNew))
    }
  }
  on.exit(if (smw_cached) clearSigmaInvSMWcache(), add = TRUE)
  cat(sprintf("SMW cache for variance-ratio Sigma^-1 calls: %s\n",
              if (smw_cached) "ENABLED" else "disabled (falling back to getSigma_X/G_multiV)"))

  Sigma_iX_noLOCO <- if (smw_cached) {
    applySigmaInvSMW_multiV_mat(X)
  } else {
    getSigma_X_multiV(W, tauVecNew, X, maxiterPCG, tolPCG, LOCO = FALSE)
  }

  ## Hoist the constant pieces of the var1 / var1GE / I_21 quadratic form
  ## (X and Sigma_iX_noLOCO don't change across markers or eMat columns) so
  ## solve(t(X) %*% Sigma_iX) is computed once instead of per iteration.
  XtSigma_iX     <- t(X) %*% Sigma_iX_noLOCO   # p x p
  XtSigma_iX_inv <- solve(XtSigma_iX)          # p x p

  y <- obj.glmm.null$y

  if (any(duplicated(obj.glmm.null$sampleID))) {
    dupSampleIndex <- as.numeric(factor(obj.glmm.null$sampleID, levels = unique(obj.glmm.null$sampleID)))
  }


  ## randomize the marker orders to be tested
  if (FALSE) {
    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
      sparseSigma <- getSparseSigma(
        bedFile = bedFile, bimFile = bimFile, famFile = famFile,
        outputPrefix = varRatioOutFile,
        sparseGRMFile = sparseGRMFile,
        sparseGRMSampleIDFile = sparseGRMSampleIDFile,
        numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
        relatednessCutoff = relatednessCutoff,
        minMAFforGRM = minMAFforGRM,
        nThreads = nThreads,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        obj.glmm.null = obj.glmm.null,
        W = W, tauVecNew = tauVecNew
      )
      if (length(tauVecNew) > 2) {
        sparseSigma <- sparseSigma + getProdTauKmat(tauVecNew[3:length(tauVecNew)])
      }
    }
  }


  mMarkers <- gettotalMarker()
  listOfMarkersForVarRatio <- list()
  MACvector <- getMACVec()
  isVarianceRatioinGeno <- getIsVarRatioGeno()

  if (isVarianceRatioinGeno) {
    MACvector_forVarRatio <- getMACVec_forVarRatio()
    Indexvector_forVarRatio <- getIndexVec_forVarRatio()
    cat("length(MACvector): ", length(MACvector), "\n")
    cat("length(MACvector_forVarRatio): ", length(MACvector_forVarRatio), "\n")

    if (length(MACvector_forVarRatio) > 0) {
      MACdata <- data.frame(MACvector = MACvector_forVarRatio, geno_ind = rep(1, length(MACvector_forVarRatio)), indexInGeno = seq(1, length(MACvector_forVarRatio)))
    } else {
      stop("No markers were found for variance ratio estimation. Please make sure there are at least 200 markers in each MAC category\n")
    }
  } else {
    MACdata <- data.frame(MACvector = MACvector, geno_ind = rep(0, length(MACvector)), indexInGeno = seq(1, length(MACvector)))
  }

  if (!isCateVarianceRatio) {
    cat("Only one variance ratio will be estimated using randomly selected markers with MAC >= 20\n")
    MACindex <- 1:nrow(MACdata)
    listOfMarkersForVarRatio[[1]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    cateVarRatioIndexVec <- c(1)
  } else {
    cat("Categorical variance ratios will be estimated.\n")

    if (is.null(cateVarRatioIndexVec)) {
      cateVarRatioIndexVec <- rep(1, length(cateVarRatioMinMACVecExclude))
    }
    numCate <- length(cateVarRatioIndexVec)
    for (i in 1:(numCate - 1)) {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[i] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[i])
      listOfMarkersForVarRatio[[i]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    }

    if (length(cateVarRatioMaxMACVecInclude) == (numCate - 1)) {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate])
    } else {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[numCate])
    }

    listOfMarkersForVarRatio[[numCate]] <- sample(MACindex, size = length(MACindex), replace = FALSE)

    for (k in 1:length(cateVarRatioIndexVec)) {
      if (k <= length(cateVarRatioIndexVec) - 1) {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], "\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], " is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      } else {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC  is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }
    }
  } # if(!isCateVarianceRatio){



  b <- as.numeric(factor(obj.glmm.null$sampleID, levels = unique(obj.glmm.null$sampleID)))
  I_mat <- Matrix::sparseMatrix(i = 1:length(b), j = b, x = rep(1, length(b)))
  I_mat <- 1.0 * I_mat

  freqVec <- getAlleleFreqVec()
  Nnomissing <- length(mu)
  varRatioTable <- NULL


  # uniqsampleind = which(!duplicated(obj.glmm.null$sampleID))
  # var_weights_sample = var_weights[uniqsampleind]


  Vsample0 <- as.vector(t(obj.noK$V) %*% I_mat)
  Xsample0 <- obj.glmm.null$sampleXMat
  XVsample0 <- t(Xsample0 * Vsample0)
  XVXsample0 <- t(Xsample0) %*% (t(XVsample0))
  XVXsample_inv0 <- solve(XVXsample0)
  XXVXsample_inv0 <- Xsample0 %*% XVXsample_inv0
  XVX_inv_XVsample0 <- XXVXsample_inv0 * Vsample0
  # XVsample0_e_list = list()
  # for(ne in 1:ncol(obj.glmm.null$eMat)){
  #  evec = obj.glmm.null$eMat[,ne]
  #  Vsample0_e = as.vector(t(obj.noK$V * evec) %*% I_mat)
  #
  #    XVsample0_e =  t(Xsample0 * (Vsample0_e))
  #    XVsample0_e_list[[ne]] = XVsample0_e
  #  }

  for (k in 1:length(listOfMarkersForVarRatio)) {
    if (cateVarRatioIndexVec[k] == 1) {
      numMarkers0 <- numMarkers
      ## Preallocate accumulators to the marker pool size and track n_acc
      ## accepted markers; trim before summarizing. Avoids quadratic c()/rbind
      ## growth across the marker loop.
      N_max <- length(listOfMarkersForVarRatio[[k]])
      varRatio_sparseGRM_vec   <- numeric(N_max)
      varRatio_NULL_vec        <- numeric(N_max)
      varRatio_NULL_sample_vec <- numeric(N_max)
      varRatio_NULL_noXadj_vec <- numeric(N_max)
      n_acc <- 0L

      if (!is.null(obj.glmm.null$eMat)) {
        nE <- ncol(obj.glmm.null$eMat)
        varRatio_NULL_eg_mat     <- matrix(0, nrow = N_max, ncol = nE)
        varRatio_NULL_eg_vec     <- NULL
        varRatio_sparse_eg_mat   <- matrix(0, nrow = N_max, ncol = nE)
        varRatio_sparse_eg_vec   <- NULL
        varModel_egcondg_mat     <- matrix(0, nrow = N_max, ncol = nE)
        varSW_egcondg_mat        <- matrix(0, nrow = N_max, ncol = nE)
        varSWtoModel_egcondg_mat <- NULL
      }

      indexInMarkerList <- 1
      numTestedMarker <- 0
      ratioCV <- ratioCVcutoff + 0.1

      ## Per-category timing accumulators (seconds). Printed at end of category.
      tcat_loop   <- 0  # total wall time of the marker loop body
      tcat_geno   <- 0  # Get_OneSNP_Geno reads
      tcat_sig_G  <- 0  # Sigma^-1 G  (one call per accepted marker)
      tcat_ctx    <- 0  # entire per-marker context for(ne) loop
      tcat_sig_GE <- 0  # Sigma^-1 GE_tilde (one call per context per marker), summed
      tcat_t0     <- Sys.time()

      while (ratioCV > ratioCVcutoff) {
        while (numTestedMarker < numMarkers0) {
          macdata_i <- listOfMarkersForVarRatio[[k]][indexInMarkerList]
          i <- (MACdata$indexInGeno)[macdata_i]
          genoInd <- (MACdata$geno_ind)[macdata_i]
          cat(i, "th marker in geno ", genoInd, "\n")
          cat("MAC: ", (MACdata$MACvector)[macdata_i], "\n")
          .tg0 <- Sys.time()
          if (genoInd == 0) {
            G0 <- Get_OneSNP_Geno(i - 1)
          } else if (genoInd == 1) {
            G0 <- Get_OneSNP_Geno_forVarRatio(i - 1)
          }
          tcat_geno <- tcat_geno + as.numeric(difftime(Sys.time(), .tg0, units = "secs"))

          # if(sum(duplicated(obj.glmm.null$sampleID)) > 0){
          if (sum(G0) / (2 * length(G0)) > 0.5) {
            G0 <- 2 - G0
          }
          G0sample <- G0
          G0 <- as.numeric(I_mat %*% G0sample)

          CHR <- bimPlink[Indexvector_forVarRatio[i] + 1, 1]
          cat("CHR ", CHR, "\n")
          if (sum(G0) / (2 * length(G0)) > 0.5) {
            G0 <- 2 - G0
          }
          NAset <- which(G0 == 0)
          AC <- sum(G0)
          indexInMarkerList <- indexInMarkerList + 1
          if ((CHR >= 1 & CHR <= 22 & AC > 0 & AC < length(G0)) | includeNonautoMarkersforVarRatio) {
            AF <- AC / (2 * Nnomissing)
            n_acc <- n_acc + 1L   # row index into preallocated accumulators
            if (CHR >= 1 & CHR <= 22) {
              autoMarker <- TRUE
            } else {
              autoMarker <- FALSE
            }

            G <- G0 - obj.noK$XXVX_inv %*% (obj.noK$XV %*% G0) # G1 is X adjusted
            # g = G/sqrt(AC)
            # q = innerProduct(g * sqrt(var_weights),y)
            # q = innerProduct(G,y)
            # eta = obj.glmm.null$linear.predictors
            # mu = obj.glmm.null$fitted.values
            # mu.eta = family$mu.eta(eta)
            # sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
            # W = sqrtW^2
            # W = W * var_weights
            # W = W * var_weights
            # print("W[1:10]")
            # print(W[1:10])
            # print("mu.eta[1:10]")
            # print(mu.eta[1:10])
            # print("obj.glm.null$family$variance(mu)[1:10]")
            # print(obj.glm.null$family$variance(mu)[1:10])
            # print("mu[1:100]")
            # print(mu[1:100])
            # print("y[1:100]")
            # print(y[1:100])


            set_isSparseGRM(useSparseGRMtoFitNULL)
            set_useGRMtoFitNULL(useGRMtoFitNULL)


            .tsg0 <- Sys.time()
            Sigma_iG <- if (smw_cached) {
              applySigmaInvSMW_multiV(G)
            } else {
              getSigma_G_multiV(W, tauVecNew, G, maxiterPCG, tolPCG, LOCO = FALSE)
            }
            tcat_sig_G <- tcat_sig_G + as.numeric(difftime(Sys.time(), .tsg0, units = "secs"))
            Sigma_iX <- Sigma_iX_noLOCO

            ## var1 = G' Sigma^-1 G - G' Sigma^-1 X (X' Sigma^-1 X)^-1 X' Sigma^-1 G
            ## Reorder so every intermediate is length-p (avoids the 1 x n
            ## vector that the chained '%*%' form produces).  XtSigma_iX_inv
            ## is hoisted to function scope above the marker loop.
            GtSig_iG  <- as.numeric(crossprod(G, Sigma_iG))
            GtSig_iX  <- as.vector(crossprod(G, Sigma_iX))     # length p
            tX_Sig_iG <- as.vector(crossprod(X, Sigma_iG))     # length p (reused below for I_21)
            var1      <- as.numeric(GtSig_iG - GtSig_iX %*% XtSigma_iX_inv %*% tX_Sig_iG)
            S <- innerProduct(G, obj.glmm.null$residuals * var_weights)
            p_exact <- pchisq(S^2 / var1, df = 1, lower.tail = F)
            cat("AC ", AC, "  S ", S, "  p_exact ", p_exact, "\n")
            # res_sample = as.vector(obj.glmm.null$residuals %*% I_mat)



            if (!is.null(obj.glmm.null$eMat)) {
              ## Vectorized context loop: build all GE / GE_tilde / Sigma_iGE columns
              ## once, then derive var1GE / I_21 / sandwich quantities via vector
              ## and matrix ops.  Mathematically equivalent to the prior per-ne loop;
              ## collapses nE small R-level matmuls and Sigma^-1 calls into one each.
              getilde_sample0_Mat <- NULL  # preserved for count_nb branch compatibility
              n_donors    <- ncol(I_mat)
              ## Coerce to plain length-N vectors: obj.glmm.null$residuals can be
              ## stored as an N x 1 matrix, which would broadcast incompatibly
              ## against GE_tilde_mat (N x nE) when nE > 1.
              residWeight <- as.numeric(obj.glmm.null$residuals) * as.numeric(var_weights)
              S_cell      <- as.numeric(G) * residWeight
              S_donor     <- as.numeric(t(I_mat) %*% S_cell)
              .tctx0 <- Sys.time()

              eMat_loc <- obj.glmm.null$eMat
              nE       <- ncol(eMat_loc)

              ## GE_mat: column ne is G0 * eMat[,ne].  R recycles G0 (length N) down columns.
              GE_mat       <- eMat_loc * as.numeric(G0)
              ## One projection matmul replaces nE scalar projections.
              GE_tilde_mat <- GE_mat - obj.noK$XXVX_inv %*% (obj.noK$XV %*% GE_mat)
              rm(GE_mat)   # consumed; release N x nE
              getildeMat   <- GE_tilde_mat   # used outside the loop for var2nullGE

              ## One batched Sigma^-1 call when SMW cache is available; otherwise
              ## fall back to per-column getSigma_G_multiV (no batched C++ API).
              .tsge0 <- Sys.time()
              if (smw_cached) {
                Sigma_iGE_mat <- applySigmaInvSMW_multiV_mat(GE_tilde_mat)
              } else {
                Sigma_iGE_mat <- matrix(0, nrow = nrow(GE_tilde_mat), ncol = nE)
                for (ne in seq_len(nE)) {
                  Sigma_iGE_mat[, ne] <- getSigma_G_multiV(
                    W, tauVecNew, GE_tilde_mat[, ne], maxiterPCG, tolPCG, LOCO = FALSE)
                }
              }
              tcat_sig_GE <- tcat_sig_GE + as.numeric(difftime(Sys.time(), .tsge0, units = "secs"))

              ## Quadratic forms — all length-nE vectors / nE x p matrices.
              GEtSig_iGE_vec <- colSums(GE_tilde_mat * Sigma_iGE_mat)              # nE
              GEtSig_iX_mat  <- crossprod(GE_tilde_mat, Sigma_iX)                  # nE x p
              tX_Sig_iGE_mat <- crossprod(X, Sigma_iGE_mat)                        # p x nE
              rm(Sigma_iGE_mat)  # consumed; release N x nE
              GEtSig_iX_Ainv <- GEtSig_iX_mat %*% XtSigma_iX_inv                   # nE x p
              ## var1GE[ne] = GEtSig_iGE[ne] - GEtSig_iX[ne,] %*% A^-1 %*% tX_Sig_iGE[,ne]
              var1GE_vec     <- as.numeric(GEtSig_iGE_vec -
                                           rowSums(GEtSig_iX_Ainv * t(tX_Sig_iGE_mat)))

              ## I_21[ne] = GEtSig_iG[ne] - GEtSig_iX[ne,] %*% A^-1 %*% tX_Sig_iG (constant)
              GEtSig_iG_vec  <- as.vector(crossprod(GE_tilde_mat, Sigma_iG))       # nE
              I_21_vec       <- as.numeric(GEtSig_iG_vec -
                                           GEtSig_iX_Ainv %*% tX_Sig_iG)

              I_11               <- var1
              c_coeff_vec        <- I_21_vec / I_11
              varModelGEcond_vec <- var1GE_vec - I_21_vec^2 / I_11

              ## Donor sandwich (HC1) — column-broadcast residual weights, then per-donor
              ## contributions, then centered sum-of-squares per context.
              S_GE_cell_mat <- GE_tilde_mat * residWeight                          # N x nE
              R_donor_mat   <- as.matrix(t(I_mat) %*% S_GE_cell_mat)               # D x nE
              rm(S_GE_cell_mat)  # consumed; release N x nE
              Q_donor_mat   <- R_donor_mat - outer(S_donor, c_coeff_vec)           # D x nE
              T_cond_vec    <- colSums(Q_donor_mat)                                # nE
              Q_bar_vec     <- T_cond_vec / n_donors
              centered_mat  <- sweep(Q_donor_mat, 2, Q_bar_vec, `-`)
              varSWGEcond_vec <- colSums(centered_mat^2) * n_donors /
                                 (n_donors - ncol(Sigma_iX))
              rm(centered_mat, Q_donor_mat, R_donor_mat)  # D x nE each

              ## Sparse-GRM variance ratio (per-context: no batched C++ API).
              if (useSparseGRMforVarRatio) {
                set_isSparseGRM(useSparseGRMforVarRatio)
                var2sparseGE_vec <- numeric(nE)
                for (ne in seq_len(nE)) {
                  GEt_ne           <- GE_tilde_mat[, ne]
                  Sigma_iGE_sparse <- getSigma_G_noV(W, tauVecNew, GEt_ne,
                                                    maxiterPCG, tolPCG, LOCO = FALSE)
                  var2sparseGE_vec[ne] <- as.numeric(crossprod(GEt_ne, Sigma_iGE_sparse))
                }
              } else if (any(duplicated(obj.glmm.null$sampleID))) {
                tauVal <- if (useGRMtoFitNULL) tauVecNew[3] else tauVecNew[2]
                var2sparseGE_vec <- numeric(nE)
                for (ne in seq_len(nE)) {
                  GEt_ne           <- GE_tilde_mat[, ne]
                  Sigma_iGE_sparse <- getSigma_G_V(W, tauVal, tauVecNew[1], GEt_ne,
                                                  maxiterPCG, tolPCG)
                  var2sparseGE_vec[ne] <- as.numeric(crossprod(GEt_ne, Sigma_iGE_sparse))
                }
              } else {
                var2sparseGE_vec <- var1GE_vec
              }

              tcat_ctx <- tcat_ctx + as.numeric(difftime(Sys.time(), .tctx0, units = "secs"))
            }


            G_noXadj <- as.vector(G0sample - mean(G0sample))
            G0_sample_tilde <- G0sample - XXVXsample_inv0 %*% (XVsample0 %*% G0sample)

            # if(useSparseGRMforVarRatio){
            # 	set_isSparseGRM(useSparseGRMforVarRatio)
            # 	Sigma_iG = getSigma_G_noV(W, tauVecNew, G, maxiterPCG, tolPCG, LOCO=FALSE)
            # 	var2_a = t(G) %*% Sigma_iG
            # 	var2sparseGRM = var2_a[1,1]
            # 	cat("var2sparseGRM ", var2sparseGRM, "\n")
            # 	varRatio_sparseGRM_vec = c(varRatio_sparseGRM_vec, var1/var2sparseGRM)

            # }else{
            if (any(duplicated(obj.glmm.null$sampleID))) {
              if (useGRMtoFitNULL) {
                tauVal <- tauVecNew[3]
              } else {
                tauVal <- tauVecNew[2]
              }
              # Sigma_iG = getSigma_G_V(W, tauVal, tauVecNew[1], G, maxiterPCG, tolPCG)
              # Sigma_iG = getSigma_G_multiV(W, tauVal, tauVecNew[1], G, maxiterPCG, tolPCG)
              # var2_a = t(G) %*% Sigma_iG
              if (isStoreSigma) {
                Sigma_iG <- (obj.glmm.null$spSigma) %*% G0_sample_tilde
                var2_a <- t(G0_sample_tilde) %*% Sigma_iG
              } else {
                G0_sample_tilde_I <- as.vector(I_mat %*% G0_sample_tilde)
                Sigma_iG <- if (smw_cached) {
                  applySigmaInvSMW_multiV(G0_sample_tilde_I)
                } else {
                  getSigma_G_multiV(W, tauVecNew, G0_sample_tilde_I, maxiterPCG, tolPCG, LOCO = FALSE)
                }
                var2_a <- t(G0_sample_tilde_I) %*% Sigma_iG
              }
              var2sparseGRM <- var2_a[1, 1]
              varRatio_sparseGRM_vec[n_acc] <- var1 / var2sparseGRM
            } else {
              varRatio_sparseGRM_vec[n_acc] <- 1
            }
            # }


            if (obj.glmm.null$traitType == "binary") {
              var2null <- innerProduct(mu * (1 - mu) * var_weights, G * G)
              var2null_sample <- innerProduct(as.vector(t(mu * (1 - mu) * var_weights) %*% I_mat), G0_sample_tilde * G0_sample_tilde)
              var2null_noXadj <- innerProduct(as.vector(t(mu * (1 - mu) * var_weights) %*% I_mat), G_noXadj * G_noXadj)
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  GE_tilde <- getildeMat[, ne]
                  # GE_sample_tilde = getilde_sample0_Mat[,ne]
                  var22nullGE <- innerProduct(mu * (1 - mu) * var_weights, GE_tilde * GE_tilde)
                  # var22nullGE = innerProduct(as.vector(t(mu*(1-mu)*var_weights) %*% I_mat), as.vector(GE_sample_tilde*GE_sample_tilde))


                  # var22nullGE = innerProduct(mu*(1-mu)*var_weights, GE_tilde*GE_tilde)
                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            } else if (obj.glmm.null$traitType == "quantitative") {
              var2null <- innerProduct(G, G * var_weights)
              var2null_sample <- innerProduct(G0_sample_tilde, G0_sample_tilde * var_weights)
              var2null_noXadj <- innerProduct(G_noXadj, G_noXadj * as.vector(t(var_weights) %*% I_mat))
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  GE_tilde <- getildeMat[, ne]
                  # GE_sample_tilde = getilde_sample0_Mat[,ne]
                  var22nullGE <- innerProduct(GE_tilde, GE_tilde * var_weights)
                  # var22nullGE = innerProduct(as.vector(GE_sample_tilde), as.vector(GE_sample_tilde)*as.vector(t(var_weights) %*% I_mat))
                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            } else if (obj.glmm.null$traitType == "count") {
              var2null <- innerProduct(mu * var_weights, G * G)
              # cat("mean(G0_sample_tilde) ", mean(G0_sample_tilde), "\n")
              # var2null_new = innerProduct(as.vector(t(mu*var_weights) %*% I_mat), G0_sample_tilde*G0_sample_tilde)
              var2null_sample <- innerProduct(as.vector(t(mu * var_weights) %*% I_mat), G0_sample_tilde * G0_sample_tilde)
              muI <- as.vector(t(mu) %*% I_mat) * as.vector(var_weights)
              var2null_noXadj <- innerProduct(as.vector(t(mu * var_weights) %*% I_mat), G_noXadj * G_noXadj)
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  GE_tilde <- getildeMat[, ne]
                  # GE_sample_tilde = getilde_sample0_Mat[,ne]
                  # cat("mean(GE_sample_tilde) ", mean(GE_sample_tilde), "\n")
                  var22nullGE <- innerProduct(mu * var_weights, GE_tilde * GE_tilde)
                  #var22nullGE = innerProduct(as.vector(t(mu*var_weights) %*% I_mat), as.vector(GE_sample_tilde*GE_sample_tilde))
                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            } else if (obj.glmm.null$traitType == "count_nb") {
              var2null <- innerProduct(W, G * G) ## To update
              var2null_sample <- innerProduct(as.vector(t(W) %*% I_mat), G0_sample_tilde * G0_sample_tilde)
              var2null_noXadj <- innerProduct(as.vector(t(W) %*% I_mat), G_noXadj * G_noXadj)
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  # GE_tilde = getildeMat[,ne]
                  GE_sample_tilde <- getilde_sample0_Mat[, ne]
                  # var22nullGE = innerProduct(W, GE_tilde*GE_tilde)
                  var22nullGE <- innerProduct(as.vector(t(W) %*% I_mat), as.vector(GE_sample_tilde * GE_sample_tilde))
                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            }

            cat("AC ", AC, "  var1 ", var1, "  var2null ", var2null, "\n")
            # cat("p_approx ", p_approx, "\n")
            # cat("p_approx_true ", p_approx_true, "\n")
            varRatio_NULL_vec[n_acc]        <- var1 / var2null
            varRatio_NULL_sample_vec[n_acc] <- var1 / var2null_sample
            varRatio_NULL_noXadj_vec[n_acc] <- var1 / var2null_noXadj
            if (!is.null(obj.glmm.null$eMat)) {
              varRatio_NULL_eg_mat[n_acc, ]   <- var1GE_vec / var2nullGE_vec
              varRatio_sparse_eg_mat[n_acc, ] <- var1GE_vec / var2sparseGE_vec
              varModel_egcondg_mat[n_acc, ]   <- varModelGEcond_vec
              varSW_egcondg_mat[n_acc, ]      <- varSWGEcond_vec
            }
            # indexInMarkerList = indexInMarkerList + 1
            numTestedMarker <- numTestedMarker + 1
          } else {
            indexInMarkerList <- indexInMarkerList + 1
          }

          if (indexInMarkerList - 1 == length(listOfMarkersForVarRatio[[k]])) {
            numTestedMarker <- numMarkers0
          }
        } # end of while(numTestedMarker < numMarkers)


        ## Vectors are preallocated; pass only the populated [1:n_acc] view to calCV
        ratioCV <- calCV(varRatio_NULL_noXadj_vec[seq_len(n_acc)])

        if (ratioCV > ratioCVcutoff) {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " > ", ratioCVcutoff, "\n")
          numMarkers0 <- numMarkers0 + 10
          cat("try ", numMarkers0, " markers\n")
        } else {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " < ", ratioCVcutoff, "\n")
        }

        if (indexInMarkerList - 1 == length(listOfMarkersForVarRatio[[k]])) {
          ratioCV <- ratioCVcutoff
          cat("no more markers are available in the MAC category ", k, "\n")
          print(indexInMarkerList - 1)
        }
      } # end of while(ratioCV > ratioCVcutoff)

      if (verbose) {
        tcat_loop <- as.numeric(difftime(Sys.time(), tcat_t0, units = "secs"))
        .nE_eff   <- if (!is.null(obj.glmm.null$eMat)) ncol(obj.glmm.null$eMat) else 0L
        cat(sprintf(
          "[TIMING:VR] cat=%d markers=%d nE=%d  total=%.2fs  geno=%.2fs  sigma_G=%.2fs  ctx_loop=%.2fs  sigma_GE_sum=%.2fs  other=%.2fs\n",
          k, n_acc, .nE_eff,
          tcat_loop, tcat_geno, tcat_sig_G, tcat_ctx, tcat_sig_GE,
          max(0, tcat_loop - tcat_geno - tcat_sig_G - tcat_ctx)
        ))
        if (n_acc > 0) {
          cat(sprintf(
            "[TIMING:VR] cat=%d per-marker:  geno=%.4fs  sigma_G=%.4fs  ctx_loop=%.4fs  sigma_GE_per_ctx=%.4fs\n",
            k,
            tcat_geno   / n_acc,
            tcat_sig_G  / n_acc,
            tcat_ctx    / n_acc,
            if (.nE_eff > 0) tcat_sig_GE / (n_acc * .nE_eff) else 0
          ))
        }
      }

      ## Trim preallocated accumulators to the n_acc actually-used rows before
      ## summarizing so mean()/colMeans() see only real values, not the zero-padded tail.
      varRatio_sparseGRM_vec   <- varRatio_sparseGRM_vec[seq_len(n_acc)]
      varRatio_NULL_vec        <- varRatio_NULL_vec[seq_len(n_acc)]
      varRatio_NULL_sample_vec <- varRatio_NULL_sample_vec[seq_len(n_acc)]
      varRatio_NULL_noXadj_vec <- varRatio_NULL_noXadj_vec[seq_len(n_acc)]
      if (!is.null(obj.glmm.null$eMat)) {
        varRatio_NULL_eg_mat   <- varRatio_NULL_eg_mat[seq_len(n_acc), , drop = FALSE]
        varRatio_sparse_eg_mat <- varRatio_sparse_eg_mat[seq_len(n_acc), , drop = FALSE]
        varModel_egcondg_mat   <- varModel_egcondg_mat[seq_len(n_acc), , drop = FALSE]
        varSW_egcondg_mat      <- varSW_egcondg_mat[seq_len(n_acc), , drop = FALSE]
      }

      if (n_acc > 0) {
        varRatio_sparse <- mean(varRatio_sparseGRM_vec)
        cat("varRatio_sparse", varRatio_sparse, "\n")
        varRatioTable <- rbind(varRatioTable, c(varRatio_sparse, "sparse", k))
      }
      varRatio_null <- mean(varRatio_NULL_vec)
      varRatio_null_sample <- mean(varRatio_NULL_sample_vec)
      cat("varRatio_null", varRatio_null, "\n")

      varRatio_null_noXadj <- mean(varRatio_NULL_noXadj_vec)
      cat("varRatio_null_noXadj", varRatio_null_noXadj, "\n")

      if (!is.null(obj.glmm.null$eMat)) {
        varRatio_NULL_eg_vec <- as.vector(colMeans(varRatio_NULL_eg_mat))
        varRatio_sparse_eg_vec <- as.vector(colMeans(varRatio_sparse_eg_mat))
        for (ne in 1:ncol(obj.glmm.null$eMat)) {
          varRatioTable <- rbind(varRatioTable, c(varRatio_NULL_eg_vec[ne], "null", 0))
          varRatioTable <- rbind(varRatioTable, c(varRatio_sparse_eg_vec[ne], "sparse", 0))
        }
        print(varRatio_NULL_eg_vec)
        print(varRatio_NULL_eg_mat)
      }

      varRatioTable <- rbind(varRatioTable, c(varRatio_null, "null", k))
      varRatioTable <- rbind(varRatioTable, c(varRatio_null_noXadj, "null_noXadj", k))
      varRatioTable <- rbind(varRatioTable, c(varRatio_null_sample, "null_sample", k))


      # varRatioTable = rbind(varRatioTable, c(varRatio_null_noXadj, "null", k))
    } else { # if(cateVarRatioVec[k] == 1)
      varRatioTable <- rbind(varRatioTable, c(1, "null", k))
      varRatioTable <- rbind(varRatioTable, c(1, "null_noXadj", k))
      varRatioTable <- rbind(varRatioTable, c(1, "null_sample", k))
      if (length(varRatio_sparseGRM_vec) > 0) {
        varRatioTable <- rbind(varRatioTable, c(1, "sparse", k))
      }
    }
  } # for(k in 1:length(listOfMarkersForVarRatio)){
  write.table(varRatioTable, varRatioOutFile, quote = F, col.names = F, row.names = F)
  
  if (!is.null(obj.glmm.null$eMat)) {
    #write.table(varModel_egcondg_mat, paste0(varRatioOutFile, ".varModel_egcondg_mat.txt"),quote = F, col.names = F, row.names = F)
    #write.table(varSW_egcondg_mat, paste0(varRatioOutFile, ".varSW_egcondg_mat.txt"),quote = F, col.names = F, row.names = F)
    varSWtoModel_egcondg_mat = varSW_egcondg_mat/varModel_egcondg_mat
    #write.table(varSWtoModel_egcondg_mat, paste0(varRatioOutFile, ".varSWtoModel_egcondg_mat.txt"),quote = F, col.names = F, row.names = F)
    ##threshold
    threshold = 1 + 2 * sqrt(2 / n_donors)
    median_ratios = apply(varSWtoModel_egcondg_mat, 2, median, na.rm = TRUE)
    cat("threshold for using sandwich variance:", threshold, "\n")
    cat("median sandwich/model variance ratios per context:\n")
    print(median_ratios)
    use_sandwich = median_ratios > threshold
    #for donor level context, do not use sandwich variance
    use_sandwich[obj.glmm.null$eMatIsSample] = FALSE
    cat("use_sandwich per context:\n")
    print(use_sandwich)
}

  data <- read.table(varRatioOutFile, header = F)
  print(data)
  return(use_sandwich)
}


# Fits the null glmm
glmmkin.ai_PCG_Rcpp_multiV <- function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau = c(0, 0), fixtau = c(0, 0), maxiter = 20, tol = 0.02, verbose = TRUE, nrun = 30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne, isLowMemLOCO, covarianceIdxMat = NULL, isStoreSigma = FALSE, useSparseGRMtoFitNULL = TRUE, useGRMtoFitNULL = TRUE, isSparseGRMIdentity = FALSE) {
  # Fits the null generalized linear mixed model for a poisson, binomial, and gaussian
  # Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the glmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  # Returns:
  #  model output for the null glmm

  t_begin <- proc.time()
  print(t_begin)
  subSampleInGeno <- subPheno$IndexGeno
  if (is.null(subPheno$IndexGeno)) {
    subSampleInGeno <- subPheno$IndexPheno
  }
  if (verbose) {
    print("Start reading genotype plink file here")
  }

  # print("subSampleInGeno")
  # print(subSampleInGeno)

  print(head(subPheno))
  set_dup_sample_index(as.numeric(factor(subPheno$IID, levels = unique(subPheno$IID))))


  print("length(indicatorGenoSamplesWithPheno)")
  print(length(indicatorGenoSamplesWithPheno))
  # if((!useSparseGRMtoFitNULL & useGRMtoFitNULL) | (skipVarianceRatioEstimation)){
  if (bedFile != "" & !useSparseGRMtoFitNULL & useGRMtoFitNULL) {
    print("HEREHRE")

    re1 <- system.time({
      setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    })
  }
  if (verbose) {
    print("Genotype reading is done")
  }

  if (LOCO) {
    MsubIndVec <- getQCdMarkerIndex()
    chrVec <- data.table:::fread(bimFile, header = F)[, 1]
    chrVec <- chrVec[which(MsubIndVec == TRUE)]
    updatechrList <- updateChrStartEndIndexVec(chrVec)
    LOCO <- updatechrList$LOCO
    chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
    chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
  }

  y <- fit0$y
  n <- length(y)
  X <- model.matrix(fit0)
  offset <- fit0$offset
  if (is.null(offset)) {
    offset <- rep(0, n)
  }

  var_weights <- weights(fit0)


  family <- fit0$family
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu) / mu.eta

  if (is.null(var_weights)) {
    var_weights <- rep(1, length(mu.eta))
  }
  # sqrtW = mu.eta/sqrt(family$variance(mu))
  # }else{
  sqrtW <- mu.eta / sqrt(1 / as.vector(var_weights) * family$variance(mu))
  # }
  W <- sqrtW^2

  alpha0 <- fit0$coef
  eta0 <- eta


  cat("tauInit")
  print(tauInit)

  tau[1:length(tau)] <- 0
  if (family$family %in% c("poisson", "binomial")) {
    # if(family$family %in% c("binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
    tauInit[1] <- 1
    idxtau <- which(fixtau == 0)
    cat("fixtau ", fixtau, "\n")
    cat("tauInit ", tauInit, "\n")
    cat("idxtau ", idxtau, "\n")
    if (sum(tauInit[idxtau]) == 0) {
      tau[idxtau] <- 0.1
    } else {
      tau[idxtau] <- tauInit[idxtau]
    }
  } else { #  if(family$family %in% c("poisson", "binomial")) {
    idxtau <- which(fixtau == 0)
    if (sum(tauInit[idxtau]) == 0) {
      tau[1] <- 1
      # tauInit[1] = 1
      tau[idxtau] <- var(Y) / (length(tau))
      # tau[2] = 0
      # tau[2:length(tau)] = 0
      if (abs(var(Y)) < 0.1) {
        stop("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T\n")
      }
    } else {
      tau[fixtau == 0] <- tauInit[fixtau == 0]
    }
  }

  cat("inital tau is ", tau, "\n")

  if (!is.null(covarianceIdxMat)) {
    idxtau2 <- intersect(covarianceIdxMat[, 1], idxtau)
    print("covarianceIdxMat")
    print(covarianceIdxMat)
    print("idxtau2")
    print(idxtau2)
    if (length(idxtau2) > 0) {
      tau[idxtau2] <- 0
    }
    # i_kmat = get_numofV()
    # if(i_kmat > 0){
    Kmatdiag <- getMeanDiagofKmat(LOCO)
    print(Kmatdiag)
    # }
    tau[2:length(tau)] <- tau[2:length(tau)] / Kmatdiag
  }

  print("tau")
  print(tau)


  #### set up weights for variance
  # if(!is.null(var_weights)){
  # 	set_var_weights(var_weights)
  # }


  # if(isStoreSigma){
  #  gen_sp_Sigma_multiV(W, tau)
  # }

  if (isSparseGRMIdentity) {
    tau[2] <- 0
  }

  # print("Here isStoreSigma")
  re.coef <- Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)

  # if(isStoreSigma){
  #  gen_sp_Sigma_multiV(re.coef$W, tau)
  # }
  # print("Here isStoreSigma 2")

  re <- getAIScore_multiV(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG = tolPCG, traceCVcutoff = traceCVcutoff, LOCO = FALSE)
  tau0 <- tau
  tau0_q2 <- tau[fixtau == 0]
  # tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)
  print("tau0_q2 a")
  print(tau0_q2)
  print("idxtau")
  print(idxtau)


  tau_q2 <- pmax(0, tau0_q2 + tau0_q2^2 * (re$YPAPY - re$Trace) / n)
  tau[idxtau] <- tau_q2

  if (!is.null(covarianceIdxMat)) {
    tau[idxtau[which(idxtau %in% idxtau2)]] <- 0
  }
  print("re$YPAPY")
  print(re$YPAPY)
  print("re$Trace")
  print(re$Trace)

  if (verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  maxiter_in <- maxiter
  if (isSparseGRMIdentity) {
    tau[2] <- 0
    maxiter_in <- 0
    alpha <- re.coef$alpha
    tau0 <- tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 <- eta
  }

  for (i in seq_len(maxiter_in)) {
    # W = sqrtW^2

    if (verbose) cat("\nIteration ", i, tau, ":\n")
    alpha0 <- re.coef$alpha
    tau0 <- tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 <- eta
    # use Get_Coef before getAIScore
    t_begin_Get_Coef <- proc.time()
    # if(isStoreSigma){
    #  gen_sp_Sigma_multiV(W, tau)
    # }
    # cat("eta0 ", eta0, "\n")
    re.coef <- Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)
    t_end_Get_Coef <- proc.time()
    cat("t_end_Get_Coef - t_begin_Get_Coef\n")
    print(t_end_Get_Coef - t_begin_Get_Coef)
    # if(isStoreSigma){
    #  gen_sp_Sigma_multiV(re.coef$W, tau)
    # }


    fit <- fitglmmaiRPCG_multiV(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol, traceCVcutoff = traceCVcutoff, LOCO = FALSE)

    t_end_fitglmmaiRPCG <- proc.time()
    cat("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
    print(t_end_fitglmmaiRPCG - t_end_Get_Coef)

    tau <- as.numeric(fit$tau)
    cov <- re.coef$cov
    alpha <- re.coef$alpha
    eta <- re.coef$eta
    Y <- re.coef$Y
    mu <- re.coef$mu

    mu.eta <- family$mu.eta(eta)

    if (is.null(var_weights)) {
      sqrtW <- mu.eta / sqrt(family$variance(mu))
    } else {
      sqrtW <- mu.eta / sqrt(1 / as.vector(var_weights) * family$variance(mu))
    }
    W <- sqrtW^2


    print(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol))
    cat("tau: ", tau, "\n")
    cat("tau0: ", tau0, "\n")

    # if(family$family == "gaussian"){
    # if(tau[1]<=0){
    #  tau[1] = tau[1] + 0.1
    #  #stop("ERROR! The first variance component parameter estimate is 0\n")
    # }
    # }

    # if(sum(tau[2:length(tau)]) == 0) break
    # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
    if (sum(tau[2:length(tau)]) == 0) {
      break
      # tau[2:length(tau)] = rep(0.1,length(tau)-1)
    } else {
      if (max(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol)) < tol) break

      if (max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
        i <- maxiter
        break
      }
    }
  }

  if (verbose) cat("\nFinal ", tau, ":\n")


  # if(isStoreSigma){
  #   gen_sp_Sigma_multiV(W, tau)
  # }

  # added these steps after tau is estimated 04-14-2018

  re.coef <- Get_Coef_multiV(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)

  cov <- re.coef$cov
  alpha <- re.coef$alpha
  eta <- re.coef$eta
  Y <- re.coef$Y
  mu <- re.coef$mu
  converged <- ifelse(i < maxiter, TRUE, FALSE)

  # var_weights = NULL
  # if(!is.null(var_weights)){
  #  res = (y - mu) * sqrt(var_weights)
  # }else{
  res <- y - mu
  # }

  if (family$family == "binomial") {
    mu2 <- mu * (1 - mu)
    traitType <- "binary"
  } else if (family$family == "poisson") {
    mu2 <- mu
    traitType <- "count"
  } else if (family$family == "gaussian") {
    mu2 <- rep((1 / (tau[1])), length(res))
    # mu2 = rep(1,length(res))
    traitType <- "quantitative"
  }

  # if(isCovariateTransform & hasCovariate){
  # if(!is.null(out.transform) & is.null(fit0$offset)){
  #  coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  # }else{
  #  coef.alpha = alpha
  # }

  # mu2 = mu * (1-mu)

  # if(!is.null(var_weights)){
  mu2_rescaled <- mu2 * var_weights
  y_rescaled <- y * var_weights
  mu_rescaled <- mu * var_weights
  # }else{


  # }

  if (!isCovariateOffset) {
    obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = X, traitType = traitType, isCovariateOffset = isCovariateOffset)
  } else {
    obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = Xorig, traitType = traitType, isCovariateOffset = isCovariateOffset)
  }

  glmmResult$varWeights <- var_weights
  # LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme
  # glmmResult$Sigma_iX = re.coef$Sigma_iX

  glmmResult$LOCO <- LOCO
  t_end_null <- proc.time()
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)
  # if(isStoreSigma){
  #   gen_sp_Sigma_multiV(re.coef$W, tau)
  # }
  if (!isLowMemLOCO & LOCO) {
    # if(isStoreSigma){
    #   gen_sp_Sigma_multiV(re.coef$W, tau)
    # }
    set_Diagof_StdGeno_LOCO()
    glmmResult$LOCOResult <- list()
    for (j in 1:22) {
      startIndex <- chromosomeStartIndexVec[j]
      endIndex <- chromosomeEndIndexVec[j]
      if (!is.na(startIndex) && !is.na(endIndex)) {
        cat("leave chromosome ", j, " out\n")
        setStartEndIndex(startIndex, endIndex, j - 1)
        t_begin_Get_Coef_LOCO <- proc.time()
        re.coef_LOCO <- Get_Coef_multiV(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = TRUE, var_weights = varWeights)
        t_end_Get_Coef_LOCO <- proc.time()
        cat("t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO\n")
        print(t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO)
        cov <- re.coef_LOCO$cov
        alpha <- re.coef_LOCO$alpha
        eta <- re.coef_LOCO$eta
        Y <- re.coef_LOCO$Y
        mu <- re.coef_LOCO$mu
        # mu2 = mu * (1-mu)
        # mu2 = mu

        res <- y - mu


        if (family$family == "binomial") {
          mu2 <- mu * (1 - mu)
        } else if (family$family == "poisson") {
          mu2 <- mu
        } else if (family$family == "gaussian") {
          mu2 <- rep((1 / (tau[1])), length(res))
        }


        if (!is.null(out.transform) & is.null(fit0$offset)) {
          coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
        } else {
          coef.alpha <- alpha
        }

        mu2_rescaled <- mu2 * var_weights
        mu_rescaled <- mu * var_weights
	
        if (!isCovariateOffset) {
          obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
        } else {
          obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
        }
        glmmResult$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK)
      } else {
        glmmResult$LOCOResult[[j]] <- list(isLOCO = FALSE)
      }
    }
  }

  if (isLowMemLOCO & LOCO) {
    glmmResult$chromosomeStartIndexVec <- chromosomeStartIndexVec
    glmmResult$chromosomeEndIndexVec <- chromosomeEndIndexVec
  }
  return(glmmResult)
}




# Fits the null glmm
glmmkin.ai_PCG_Rcpp_multiV_eMat <- function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau = c(0, 0), fixtau = c(0, 0), maxiter = 20, tol = 0.02, verbose = TRUE, nrun = 30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne, isLowMemLOCO, covarianceIdxMat = NULL, isStoreSigma = FALSE, useSparseGRMtoFitNULL = TRUE, useGRMtoFitNULL = TRUE, isSparseGRMIdentity = FALSE) {
  # Fits the null generalized linear mixed model for a poisson, binomial, and gaussian
  # Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the glmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  # Returns:
  #  model output for the null glmm

  t_begin <- proc.time()
  print(t_begin)
  subSampleInGeno <- subPheno$IndexGeno
  if (is.null(subPheno$IndexGeno)) {
    subSampleInGeno <- subPheno$IndexPheno
  }
  if (verbose) {
    print("Start reading genotype plink file here")
  }

  # print("subSampleInGeno")
  # print(subSampleInGeno)

  print(head(subPheno))
  set_dup_sample_index(as.numeric(factor(subPheno$IID, levels = unique(subPheno$IID))))


  print("length(indicatorGenoSamplesWithPheno)")
  print(length(indicatorGenoSamplesWithPheno))
  # if((!useSparseGRMtoFitNULL & useGRMtoFitNULL) | (skipVarianceRatioEstimation)){
  if (bedFile != "" & !useSparseGRMtoFitNULL & useGRMtoFitNULL) {
    print("HEREHRE")

    re1 <- system.time({
      setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    })
  }
  if (verbose) {
    print("Genotype reading is done")
  }

  if (LOCO) {
    MsubIndVec <- getQCdMarkerIndex()
    chrVec <- data.table:::fread(bimFile, header = F)[, 1]
    chrVec <- chrVec[which(MsubIndVec == TRUE)]
    updatechrList <- updateChrStartEndIndexVec(chrVec)
    LOCO <- updatechrList$LOCO
    chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
    chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
  }

  y <- fit0$y
  n <- length(y)
  X <- model.matrix(fit0)
  offset <- fit0$offset
  if (is.null(offset)) {
    offset <- rep(0, n)
  }

  var_weights <- weights(fit0)


  family <- fit0$family
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu) / mu.eta

  if (is.null(var_weights)) {
    var_weights <- rep(1, length(mu.eta))
  }
  # sqrtW = mu.eta/sqrt(family$variance(mu))
  # }else{
  sqrtW <- mu.eta / sqrt(1 / as.vector(var_weights) * family$variance(mu))
  # }
  W <- sqrtW^2

  alpha0 <- fit0$coef
  eta0 <- eta


  cat("tauInit")
  print(tauInit)

  tau[1:length(tau)] <- 0
  if (family$family %in% c("poisson", "binomial")) {
    # if(family$family %in% c("binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
    tauInit[1] <- 1
    idxtau <- which(fixtau == 0)
    cat("fixtau ", fixtau, "\n")
    cat("tauInit ", tauInit, "\n")
    cat("idxtau ", idxtau, "\n")
    if (sum(tauInit[idxtau]) == 0) {
      tau[idxtau] <- 0.1
    } else {
      tau[idxtau] <- tauInit[idxtau]
    }
  } else { #  if(family$family %in% c("poisson", "binomial")) {
    idxtau <- which(fixtau == 0)
    if (sum(tauInit[idxtau]) == 0) {
      tau[1] <- 1
      # tauInit[1] = 1
      tau[idxtau] <- var(Y) / (length(tau))
      # tau[2] = 0
      # tau[2:length(tau)] = 0
      if (abs(var(Y)) < 0.1) {
        stop("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T\n")
      }
    } else {
      tau[fixtau == 0] <- tauInit[fixtau == 0]
    }
  }

  cat("inital tau is ", tau, "\n")

  if (!is.null(covarianceIdxMat)) {
    idxtau2 <- intersect(covarianceIdxMat[, 1], idxtau)
    print("covarianceIdxMat")
    print(covarianceIdxMat)
    print("idxtau2")
    print(idxtau2)
    if (length(idxtau2) > 0) {
      tau[idxtau2] <- 0
    }
    # i_kmat = get_numofV()
    # if(i_kmat > 0){
    Kmatdiag <- getMeanDiagofKmat(LOCO)
    print(Kmatdiag)
    # }
    tau[2:length(tau)] <- tau[2:length(tau)] / Kmatdiag
  }

  print("tau")
  print(tau)


  #### set up weights for variance
  # if(!is.null(var_weights)){
  #     set_var_weights(var_weights)
  # }


   #if(isStoreSigma){
   #  gen_sp_Sigma_multiV(W, tau)
   #}

  if (isSparseGRMIdentity) {
    tau[2] <- 0
  }

  # print("Here isStoreSigma")
  re.coef <- Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)

  # if(isStoreSigma){
  #  gen_sp_Sigma_multiV(re.coef$W, tau)
  # }
  # print("Here isStoreSigma 2")

  re <- getAIScore_multiV_eMat(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG = tolPCG, traceCVcutoff = traceCVcutoff, LOCO = FALSE)
  tau0 <- tau
  tau0_q2 <- tau[fixtau == 0]
  # tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)
  print("tau0_q2 a")
  print(tau0_q2)
  print("idxtau")
  print(idxtau)


  tau_q2 <- pmax(0, tau0_q2 + tau0_q2^2 * (re$YPAPY - re$Trace) / n)
  tau[idxtau] <- tau_q2

  if (!is.null(covarianceIdxMat)) {
    tau[idxtau[which(idxtau %in% idxtau2)]] <- 0
  }
  print("re$YPAPY")
  print(re$YPAPY)
  print("re$Trace")
  print(re$Trace)

  if (verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  maxiter_in <- maxiter
  if (isSparseGRMIdentity) {
    tau[2] <- 0
    maxiter_in <- 0
    alpha <- re.coef$alpha
    tau0 <- tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 <- eta
  }




  for (i in seq_len(maxiter_in)) {
    # W = sqrtW^2

    if (verbose) cat("\nIteration ", i, tau, ":\n")
    alpha0 <- re.coef$alpha
    tau0 <- tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 <- eta
    # use Get_Coef before getAIScore
    t_begin_Get_Coef <- proc.time()
    # if(isStoreSigma){
    #  gen_sp_Sigma_multiV(W, tau)
    # }
    # cat("eta0 ", eta0, "\n")
    re.coef <- Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)
    t_end_Get_Coef <- proc.time()
    cat("t_end_Get_Coef - t_begin_Get_Coef\n")
    print(t_end_Get_Coef - t_begin_Get_Coef)
    # if(isStoreSigma){
    #  gen_sp_Sigma_multiV(re.coef$W, tau)
    # }


    fit <- fitglmmaiRPCG_multiV_eMat(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol, traceCVcutoff = traceCVcutoff, LOCO = FALSE)

    t_end_fitglmmaiRPCG <- proc.time()
    cat("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
    print(t_end_fitglmmaiRPCG - t_end_Get_Coef)

    tau <- as.numeric(fit$tau)
    cov <- re.coef$cov
    alpha <- re.coef$alpha
    eta <- re.coef$eta
    Y <- re.coef$Y
    mu <- re.coef$mu

    mu.eta <- family$mu.eta(eta)

    if (is.null(var_weights)) {
      sqrtW <- mu.eta / sqrt(family$variance(mu))
    } else {
      sqrtW <- mu.eta / sqrt(1 / as.vector(var_weights) * family$variance(mu))
    }
    W <- sqrtW^2


    print(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol))
    cat("tau: ", tau, "\n")
    cat("tau0: ", tau0, "\n")

    # if(family$family == "gaussian"){
    # if(tau[1]<=0){
    #  tau[1] = tau[1] + 0.1
    #  #stop("ERROR! The first variance component parameter estimate is 0\n")
    # }
    # }

    # if(sum(tau[2:length(tau)]) == 0) break
    # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
    if (sum(tau[2:length(tau)]) == 0) {
      break
      # tau[2:length(tau)] = rep(0.1,length(tau)-1)
    } else {
      if (max(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol)) < tol) break

      if (max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
        i <- maxiter
        break
      }
    }
  }

  if (verbose) cat("\nFinal ", tau, ":\n")


  # if(isStoreSigma){
  #   gen_sp_Sigma_multiV(W, tau)
  # }

  # added these steps after tau is estimated 04-14-2018

  re.coef <- Get_Coef_multiV(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)

  cov <- re.coef$cov
  alpha <- re.coef$alpha
  eta <- re.coef$eta
  Y <- re.coef$Y
  mu <- re.coef$mu
  converged <- ifelse(i < maxiter, TRUE, FALSE)

  # var_weights = NULL
  # if(!is.null(var_weights)){
  #  res = (y - mu) * sqrt(var_weights)
  # }else{
  res <- y - mu
  # }

  if (family$family == "binomial") {
    mu2 <- mu * (1 - mu)
    traitType <- "binary"
  } else if (family$family == "poisson") {
    mu2 <- mu
    traitType <- "count"
  } else if (family$family == "gaussian") {
    mu2 <- rep((1 / (tau[1])), length(res))
    # mu2 = rep(1,length(res))
    traitType <- "quantitative"
  }

  # if(isCovariateTransform & hasCovariate){
  # if(!is.null(out.transform) & is.null(fit0$offset)){
  #  coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  # }else{
  #  coef.alpha = alpha
  # }

  # mu2 = mu * (1-mu)

  # if(!is.null(var_weights)){
  mu2_rescaled <- mu2 * var_weights
  y_rescaled <- y * var_weights
  mu_rescaled <- mu * var_weights
  # }else{


  # }

  if (!isCovariateOffset) {
    obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = X, traitType = traitType, isCovariateOffset = isCovariateOffset)
  } else {
    obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = Xorig, traitType = traitType, isCovariateOffset = isCovariateOffset)
  }

  glmmResult$varWeights <- var_weights
  # LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme
  # glmmResult$Sigma_iX = re.coef$Sigma_iX

  glmmResult$LOCO <- LOCO
  t_end_null <- proc.time()
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)
  # if(isStoreSigma){
  #   gen_sp_Sigma_multiV(re.coef$W, tau)
  # }
  if (!isLowMemLOCO & LOCO) {
    # if(isStoreSigma){
    #   gen_sp_Sigma_multiV(re.coef$W, tau)
    # }
    set_Diagof_StdGeno_LOCO()
    glmmResult$LOCOResult <- list()
    for (j in 1:22) {
      startIndex <- chromosomeStartIndexVec[j]
      endIndex <- chromosomeEndIndexVec[j]
      if (!is.na(startIndex) && !is.na(endIndex)) {
        cat("leave chromosome ", j, " out\n")
        setStartEndIndex(startIndex, endIndex, j - 1)
        t_begin_Get_Coef_LOCO <- proc.time()
        re.coef_LOCO <- Get_Coef_multiV(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = TRUE, var_weights = varWeights)
        t_end_Get_Coef_LOCO <- proc.time()
        cat("t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO\n")
        print(t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO)
        cov <- re.coef_LOCO$cov
        alpha <- re.coef_LOCO$alpha
        eta <- re.coef_LOCO$eta
        Y <- re.coef_LOCO$Y
        mu <- re.coef_LOCO$mu
        # mu2 = mu * (1-mu)
        # mu2 = mu

        res <- y - mu


        if (family$family == "binomial") {
          mu2 <- mu * (1 - mu)
        } else if (family$family == "poisson") {
          mu2 <- mu
        } else if (family$family == "gaussian") {
          mu2 <- rep((1 / (tau[1])), length(res))
        }


        if (!is.null(out.transform) & is.null(fit0$offset)) {
          coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
        } else {
          coef.alpha <- alpha
        }

        mu2_rescaled <- mu2 * var_weights
        mu_rescaled <- mu * var_weights


        if (!isCovariateOffset) {
          obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
        } else {
          obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
        }
        glmmResult$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK)
      } else {
        glmmResult$LOCOResult[[j]] <- list(isLOCO = FALSE)
      }
    }
  }

  if (isLowMemLOCO & LOCO) {
    glmmResult$chromosomeStartIndexVec <- chromosomeStartIndexVec
    glmmResult$chromosomeEndIndexVec <- chromosomeEndIndexVec
  }
  return(glmmResult)
}
