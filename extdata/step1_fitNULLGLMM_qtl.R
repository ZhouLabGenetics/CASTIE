#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

cat("===== SAIGE-QTL Dynamic wrapper start =====\n")


## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile", type="character",default="",
    help="Path to plink file for creating the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF and maxMissingRate can be used to specify the maximum missing rates  of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios"),
  make_option("--bedFile", type="character",default="",
    help="Path to bed file. If plinkFile is specified, 'plinkFile'.bed will be used"),
  make_option("--bimFile", type="character",default="",
    help="Path to bim file. If plinkFile is specified, 'plinkFile'.bim will be used"),
  make_option("--famFile", type="character",default="",
    help="Path to fam file. If plinkFile is specified, 'plinkFile'.fam will be used"),
  make_option("--phenoFile", type="character", default="",
    help="Required. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarColList. All categorical covariates must also be included in covarColList."),
  make_option("--phenoCol", type="character", default="",
    help="Required. Column name for phenotype to be tested in the phenotype file, e.g CAD"),
  make_option("--isRemoveZerosinPheno", type="logical", default=FALSE,  
    help="Optional. Whether to remove zeros in the phenotype"),
  make_option("--traitType", type="character", default="binary",help="Required. binary or quantitative [default=binary]"),
  make_option("--invNormalize", type="logical",default=FALSE,
    help="Optional. Only for quantitative. Whether to perform the inverse normalization for the phenotype [default='FALSE']"),
  make_option("--covarColList", type="character", default="",
    help="List of covariates (comma separated)"),
  make_option("--sampleCovarColList", type="character", default="",
    help="List of covariates that are on sample level (comma separated)"),
  make_option("--dynamicCovarColList", type="character", default="",
    help="List of covariates that are for dynamic qtls (comma separated)"), 
  make_option("--longlCol", type="character", default="",
    help=""),	      
  make_option("--qCovarColList", type="character", default="",
    help="List of categorical covariates (comma separated). All categorical covariates must also be in covarColList"),
  make_option("--sampleIDColinphenoFile", type="character", default="IID",
    help="Required. Column name of sample IDs in the phenotype file, e.g. IID"),
  make_option("--tol", type="numeric", default=0.02,
    help="Optional. Tolerance for fitting the null GLMM to converge [default=0.02]."),
  make_option("--maxiter", type="integer", default=500,
    help="Optional. Maximum number of iterations used to fit the null GLMM [default=20]."),
  make_option("--tolPCG", type="numeric", default=1e-5,
    help="Optional. Tolerance for PCG to converge [default=1e-5]."),
  make_option("--maxiterPCG", type="integer", default=500,
    help="Optional. Maximum number of iterations for PCG [default=500]."),		
  make_option("--nThreads", type="integer", default=1,
    help="Optional. Number of threads (CPUs) to use [default=1]."),
  make_option("--SPAcutoff", type="numeric", default=2,
    help="Optional. Cutoff for the deviation of score test statistics from mean in the unit of sd to perform SPA [default=2]."),
  make_option("--numRandomMarkerforVarianceRatio", type="integer", default=30,
    help="Optional. An integer greater than 0. Number of markers to be randomly selected for estimating the variance ratio. The number will be automatically added by 10 until the coefficient of variantion (CV) for the variance ratio estimate is below ratioCVcutoff [default=30]."),
  make_option("--skipModelFitting", type="logical", default=FALSE,
    help="Optional. Whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"),
  make_option("--skipVarianceRatioEstimation", type="logical", default=FALSE,
    help="Optional. Whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"),
  make_option("--memoryChunk", type="numeric", default=2,
   help="Optional. Size (Gb) for each memory chunk [default=2]"),
  make_option("--tauInit", type="character", default="0,0",
   help="Optional. Initial values for tau. [default=0,0]"),
  make_option("--LOCO", type="logical", default=TRUE,
    help="Whether to apply the leave-one-chromosome-out (LOCO) approach when fitting the null model using the full GRM [default=TRUE]."),
  make_option("--isLowMemLOCO", type="logical", default=FALSE,
    help="Whehter to output the model file by chromosome when LOCO=TRUE. If TRUE, the memory usage in Step 1 and Step 2 will be lower [default=FALSE]"),
  make_option("--traceCVcutoff", type="numeric", default=0.0025,
    help="Optional. Threshold for coefficient of variation (CV) for the trace estimator. Number of runs for trace estimation will be increased until the CV is below the threshold [default=0.0025]."),
  make_option("--nrun", type="numeric", default=30,
    help="Number of rums in trace estimation. [default=30]"),	       
  make_option("--ratioCVcutoff", type="numeric", default=0.001,
    help="Optional. Threshold for coefficient of variation (CV) for estimating the variance ratio. The number of randomly selected markers will be increased until the CV is below the threshold [default=0.001]"),
  make_option("--outputPrefix", type="character", default="~/",
    help="Required. Path and prefix of the output files [default='~/']"),
  make_option("--outputPrefix_varRatio", type="character", default="",
    help="Optional. Path and prefix of the output the variance ratio file. if not specified, it will be the same as the outputPrefix"),
  make_option("--IsOverwriteVarianceRatioFile", type="logical", default=FALSE,
    help="Optional. Whether to overwrite the variance ratio file if the file exist.[default='FALSE']"),
  make_option("--sparseGRMFile", type="character", default="",
   help="Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]"),
  make_option("--sparseGRMSampleIDFile", type="character", default="",
   help="Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM [default=NULL]"),
  make_option("--isCateVarianceRatio", type="logical", default=FALSE,
    help="Required. Whether to estimate variance ratio based on different MAC categories. If yes, variance ratio will be estiamted for multiple MAC categories corresponding to cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude. Currently, if isCateVarianceRatio=TRUE, then LOCO=FALSE [default=FALSE]"),
  make_option("--relatednessCutoff", type="numeric", default=0,
    help="Optional. Threshold (minimum relatedness coefficient) to treat two samples as unrelated when the sparse GRM is used [default=0]"),
  make_option("--cateVarRatioMinMACVecExclude",type="character", default="10,20.5",
    help="Optional. vector of float. Lower bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='10,20.5']"),
  make_option("--cateVarRatioMaxMACVecInclude",type="character", default="20.5",
    help="Optional. vector of float. Higher bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='20.5']"),    
  make_option("--isCovariateTransform", type="logical", default=TRUE,
    help="Optional. Whether use qr transformation on covariates [default='TRUE']."),
  make_option("--isDiagofKinSetAsOne", type="logical", default=FALSE,
    help="Optional. Whether to set the diagnal elements in GRM to be 1 [default='FALSE']."),
  make_option("--useSparseGRMtoFitNULL", type="logical", default=FALSE,
    help="Optional. Whether to use sparse GRM to fit the null model [default='FALSE']."),
  make_option("--useSparseGRMforVarRatio", type="logical", default=FALSE, 
    help="Optional. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE"),
  make_option("--minMAFforGRM", type="numeric", default=0.01,
    help="Optional. Minimum MAF of markers used for GRM"),
  make_option("--maxMissingRateforGRM", type="numeric", default=0.15,
    help="Optional. Maximum missing rate of markers used for GRM"),
  make_option("--minCovariateCount", type="numeric", default=-1,
    help="Optional. Binary covariates with a count less than minCovariateCount will be excluded from the model to avoid convergence issues [default=-1] (no covariates will be excluded)."),
  make_option("--includeNonautoMarkersforVarRatio", type="logical", default=FALSE,
    help="Optional. Whether to allow for non-autosomal markers for variance ratio. [default, 'FALSE']"),
  make_option("--FemaleOnly", type="logical", default=FALSE,
    help="Optional. Whether to run Step 1 for females only [default=FALSE]. if TRUE, --sexCol and --FemaleCode need to be specified"), 
  make_option("--MaleOnly", type="logical", default=FALSE,
    help="Optional. Whether to run Step 1 for males only [default=FALSE]. if TRUE, --sexCol and --MaleCode need to be specified"),   
  make_option("--sexCol", type="character", default="",
   help="Optional. Column name for sex in the phenotype file, e.g Sex"),
  make_option("--FemaleCode", type="character", default="1",
   help="Optional. Values in the column for sex in the phenotype file are used for females [default, '1']"),
  make_option("--MaleCode", type="character", default="0",
   help="Optional. Values in the column for sex in the phenotype file are used for males [default, '0']"),
  make_option("--isCovariateOffset", type="logical", default=FALSE,
   help="Optional. Whether to estimate fixed effect coeffciets. [default, 'TRUE']"),
  make_option("--SampleIDIncludeFile", type="character",default="",
    help="Path to the file that contains one column for IDs of samples who will be include for null model fitting."),
 make_option("--VmatFilelist", type="character", default="",
    help="List of additional V (comma separated)"),		    
 make_option("--VmatSampleFilelist", type="character", default="",
    help="List of additional V (comma separated)"),
   make_option("--useGRMtoFitNULL", type="logical", default=TRUE,help=""),
  make_option("--offsetCol", type="character", default=NULL,
   help="offset column"),
       make_option("--varWeightsCol", type="character", default=NULL,
   help="variance weight column"),
  make_option("--isStoreSigma", type="logical", default=TRUE,
   help="Optional. Whether to store the inv Sigma matrix. [default, 'TRUE']"),
  make_option("--isShrinkModelOutput", type="logical", default=TRUE,
   help="Optional. Whether to remove unnecessary objects for step2 from the model output. [default, 'TRUE']"),
  make_option("--use_PCG", type="logical", default=FALSE,
   help="Optional. Whether to force PCG (Case 3) solver for null model fitting. Automatically set to TRUE when variance components are out of bounds. [default, 'FALSE']"),
  make_option("--isWriteReport", type="logical", default=FALSE,
   help="Optional. Whether to save a fitting report (solver used, convergence, offset flag) to ./report/. [default, 'FALSE']"),
  make_option("--initialSubSampleProp",
  type="numeric", default=1,
      help="Optional. The proportion of subsamples used for estimate the inital values for variance covariance parameters"),
  make_option("--library", type="character", default="",
      help="Optional. Path to the library directory where SAIGEQTL is installed")
)


## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

## Load SAIGEQTL with optional library path
if(opt$library != ""){
  suppressPackageStartupMessages({
    library(SAIGEQTL, lib.loc=opt$library)
  })
  cat("Loaded SAIGEQTL from library path:", opt$library, "\n")
} else {
  suppressPackageStartupMessages({
    library(SAIGEQTL)
  })
  cat("Loaded SAIGEQTL from default library path\n")
}

print(sessionInfo())
print(opt)

covars <- strsplit(opt$covarColList,",")[[1]]
qcovars <- strsplit(opt$qCovarColList,",")[[1]]
scovars <- strsplit(opt$sampleCovarColList,",")[[1]]
ecovars <- strsplit(opt$dynamicCovarColList,",")[[1]]


convertoNumeric = function(x,stringOutput){
        y= tryCatch(expr = as.numeric(x),warning = function(w) {return(NULL)})
        if(is.null(y)){
                stop(stringOutput, " is not numeric\n")
        }else{
                cat(stringOutput, " is ", y, "\n")
        }
        return(y)
}

tauInit <- convertoNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")
cateVarRatioMinMACVecExclude <- convertoNumeric(x=strsplit(opt$cateVarRatioMinMACVecExclude,",")[[1]], "cateVarRatioMinMACVecExclude")
cateVarRatioMaxMACVecInclude <- convertoNumeric(x=strsplit(opt$cateVarRatioMaxMACVecInclude,",")[[1]], "cateVarRatioMaxMACVecInclude")

BLASctl_installed <- require(RhpcBLASctl)
##by Alex Petty @pettyalex
if (BLASctl_installed){
  # Set number of threads for BLAS to 1, this step does not benefit from multithreading or multiprocessing
    original_num_threads <- blas_get_num_procs()
    blas_set_num_threads(1)
}

#initialSubSampleProp=0.01
#initialSubSampleTimes=10
initialSubSampleProp=opt$initialSubSampleProp

initialSubSampleTimes=1
    tauTable = NULL	
if(initialSubSampleProp < 1){
       
       phenoFile=opt$phenoFile
       sampleIDColinphenoFile=opt$sampleIDColinphenoFile 
       phenoCol=opt$phenoCol
       outputPrefix=opt$outputPrefix
       eCovarCol = ecovars 
       covarColList = covars
       longlCol = opt$longlCol 
       offsetCol = opt$offsetCol
       qCovarCol = qcovars
       varWeightsCol = opt$varWeightsCol	

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


    ## check whether the phenotype file is large
    cmd <- paste0("du ", phenoFile, "| awk '{print $1}' > ", outputPrefix, "_", phenoCol, "_size_temp")
    system(cmd)
    datasize <- data.table:::fread(paste0(outputPrefix, "_", phenoCol, "_size_temp"), header = F, data.table = F)
    isphenoFileLarge <- FALSE
    if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
      if (datasize[1, 1] > 200000) {
        isphenoFileLarge <- TRUE
      }
    } else {
      if (datasize[1, 1] > 500000) {
        isphenoFileLarge <- TRUE
      }
    }

    if (isphenoFileLarge) {
      if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
        cmd <- paste0("gunzip -c ", phenoFile, "| head -n 1 | sed 's/\\t/\\n/g' | sed 's/\ /\\n/g' | awk '{print $1\"\\t\"NR}' > ", outputPrefix, "_", phenoCol, "_lineNum_temp")
        system(cmd)
      } else {
        cmd <- paste0("cat ", phenoFile, "| head -n 1 | sed 's/\\t/\\n/g' | sed 's/\ /\\n/g' | awk '{print $1\"\\t\"NR}' > ", outputPrefix, "_", phenoCol, "_lineNum_temp")
        system(cmd)
      }

      # print(cmd)

      checkColListDataFrame <- data.frame(colna = checkColList)
      phenoFilephenoCol_lineNum <- data.table::fread(paste0(outputPrefix, "_", phenoCol, "_lineNum_temp"), header = F, data.table = F)

      phenoFilephenoCol_lineNum_checkColList <- merge(checkColListDataFrame, phenoFilephenoCol_lineNum, by.x = 1, by.y = 1)

      write.table(phenoFilephenoCol_lineNum_checkColList[, 2], paste0(outputPrefix, "_", phenoCol, "_colnames_subset_temp"), quote = F, col.names = F, row.names = F)


      if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
        cmdb <- paste0("cut -f $(tr '\\n' ',' < ", outputPrefix, "_", phenoCol, "_colnames_subset_temp | sed 's/,$//') <(gunzip -c", phenoFile, ")> ", outputPrefix, "_", phenoCol, "_subcols_temp")
      } else {
        cmdb <- paste0("cut -f $(tr '\\n' ',' < ", outputPrefix, "_", phenoCol, "_colnames_subset_temp | sed 's/,$//') ", phenoFile, "> ", outputPrefix, "_", phenoCol, "_subcols_temp")
      }

      # print(cmdb)
      system(cmdb)

      phenoFiletemp <- paste0(outputPrefix, "_", phenoCol, "_subcols_temp")


      data <- data.table:::fread(phenoFiletemp,
        header = T,
        stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table = F
      )
      # data = data.frame(ydat)

      file.remove(paste0(outputPrefix, "_", phenoCol, "_colnames_subset_temp"))
      file.remove(paste0(outputPrefix, "_", phenoCol, "_lineNum_temp"))
      #file.remove(paste0(outputPrefix, "_", phenoCol, "_subcols_temp"))
    } else { # !isphenoFileLarge

      if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
        data <- data.table:::fread(
          cmd = paste0(
            "gunzip -c ",
            phenoFile
          ), header = T, stringsAsFactors = FALSE,
          colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
        )
      } else {
        data <- data.table:::fread(phenoFile,
          header = T,
          stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
        )
      }
    }

    file.remove(paste0(outputPrefix, "_", phenoCol, "_size_temp"))
    ID = data[,which(colnames(data) == sampleIDColinphenoFile)]



     for(init_i in 1:initialSubSampleTimes){
       set.seed(init_i)
        	
	# Check for duplicates in the IID column
	if (anyDuplicated(ID) > 0) {
  	# If IDs are duplicated, subset within each unique ID group
  	  sampled_indices <- unlist(
          lapply(unique(ID), function(id) {
          # Subset rows for the current ID
          id_rows <- which(ID == id)
          # Sample 1% of rows for the current ID
          sample(id_rows, size = max(1, floor(0.01 * length(id_rows)), na.rm = TRUE))
          })
          )
       } else {
       # If IDs are not duplicated, directly bootstrap 1% of rows
        sampled_indices <- sample(seq_len(nrow(data)), size = max(1, floor(0.01 * nrow(data)), na.rm = TRUE))
       }
       # Save the indices to a file
	datasubset = data[sampled_indices, ]
	datasubsetFile = paste0(outputPrefix, "_", phenoCol, "_phenosubset_", init_i,".txt")
	data.table::fwrite(datasubset, file=datasubsetFile, quote=F, sep="\t", row.names=F, col.names=T)

#fit_success <- TRUE
#tryCatch({
#set seed
outputprefixsubset = paste0(outputPrefix, "_", phenoCol, "_phenosubset_", init_i)
fitNULLGLMM_multiV(plinkFile=opt$plinkFile,
            bedFile=opt$bedFile,
            bimFile=opt$bimFile,
            famFile=opt$famFile,
            useSparseGRMtoFitNULL=opt$useSparseGRMtoFitNULL,
            sparseGRMFile=opt$sparseGRMFile,
            sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
            phenoFile = datasubsetFile,
            phenoCol = opt$phenoCol,
            isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
            sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
            traitType = opt$traitType,
            outputPrefix = outputprefixsubset,
            isCovariateOffset=opt$isCovariateOffset,
            nThreads = opt$nThreads,
            useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
            invNormalize = opt$invNormalize,
            covarColList = covars,
            qCovarCol = qcovars,
            tol=opt$tol,
            maxiter=opt$maxiter,
            tolPCG=opt$tolPCG,
            maxiterPCG=opt$maxiterPCG,
            SPAcutoff = opt$SPAcutoff,
            numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
            skipModelFitting = opt$skipModelFitting,
            skipVarianceRatioEstimation = TRUE,
            memoryChunk = opt$memoryChunk,
            tauInit = tauInit,
            LOCO = opt$LOCO,
            isLowMemLOCO = opt$isLowMemLOCO,
            traceCVcutoff = opt$traceCVcutoff,
            nrun = opt$nrun,
            ratioCVcutoff = opt$ratioCVcutoff,
            outputPrefix_varRatio = opt$outputPrefix_varRatio,
            IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
            relatednessCutoff = opt$relatednessCutoff,
            isCateVarianceRatio = opt$isCateVarianceRatio,
            cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
            cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
            isCovariateTransform = opt$isCovariateTransform,
            isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
            minMAFforGRM = opt$minMAFforGRM,
            maxMissingRateforGRM = opt$maxMissingRateforGRM,
            minCovariateCount=opt$minCovariateCount,
            includeNonautoMarkersforVarRatio=opt$includeNonautoMarkersforVarRatio,
            sexCol=opt$sexCol,
            FemaleCode=opt$FemaleCode,
            FemaleOnly=opt$FemaleOnly,
            MaleCode=opt$MaleCode,
            MaleOnly=opt$MaleOnly,
            SampleIDIncludeFile=opt$SampleIDIncludeFile,
            VmatFilelist=opt$VmatFilelist,
            VmatSampleFilelist=opt$VmatSampleFilelist,
            longlCol=opt$longlCol,
            useGRMtoFitNULL=opt$useGRMtoFitNULL,
            offsetCol=opt$offsetCol,
            varWeightsCol=opt$varWeightsCol,
            sampleCovarCol=scovars,
            isStoreSigma=opt$isStoreSigma,
      isShrinkModelOutput=opt$isShrinkModelOutput,
      eCovarCol=ecovars
        )


#}, error = function(e) {
#  message("Initial model failed with error: ", e$message)
#    fit_success <<- FALSE  # Track failure
#}) 


  my_env = new.env()
  load(paste0(outputprefixsubset, ".rda"), envir = my_env)
  modglmm = my_env$modglmm
  print(modglmm$theta)

if(!opt$isCovariateOffset){
  if(sum(modglmm$theta[2:length(modglmm$theta)]) <= 0 || sum(modglmm$theta[2:length(modglmm$theta)]) > 1){
        cat("Theta out of bounds, trying PCG solver first...\n")
        pcg_success <- tryCatch({
          set_use_PCG(TRUE)
          set.seed(1)
          fitNULLGLMM_multiV(plinkFile=opt$plinkFile,
              bedFile=opt$bedFile,
              bimFile=opt$bimFile,
              famFile=opt$famFile,
              useSparseGRMtoFitNULL=opt$useSparseGRMtoFitNULL,
              sparseGRMFile=opt$sparseGRMFile,
              sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
              phenoFile = datasubsetFile,
              phenoCol = opt$phenoCol,
              isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
              sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
              traitType = opt$traitType,
              outputPrefix = paste0(outputprefixsubset, ".pcg"),
              isCovariateOffset = FALSE,
              nThreads = opt$nThreads,
              useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
              invNormalize = opt$invNormalize,
              covarColList = covars,
              qCovarCol = qcovars,
              tol=opt$tol,
              maxiter=opt$maxiter,
              tolPCG=opt$tolPCG,
              maxiterPCG=opt$maxiterPCG,
              SPAcutoff = opt$SPAcutoff,
              numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
              skipModelFitting = opt$skipModelFitting,
              skipVarianceRatioEstimation = TRUE,
              memoryChunk = opt$memoryChunk,
              tauInit = tauInit,
              LOCO = opt$LOCO,
              isLowMemLOCO = opt$isLowMemLOCO,
              traceCVcutoff = opt$traceCVcutoff,
              nrun = opt$nrun,
              ratioCVcutoff = opt$ratioCVcutoff,
              outputPrefix_varRatio = opt$outputPrefix_varRatio,
              IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
              relatednessCutoff = opt$relatednessCutoff,
              isCateVarianceRatio = opt$isCateVarianceRatio,
              cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
              cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
              isCovariateTransform = opt$isCovariateTransform,
              isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
              minMAFforGRM = opt$minMAFforGRM,
              maxMissingRateforGRM = opt$maxMissingRateforGRM,
              minCovariateCount=opt$minCovariateCount,
              includeNonautoMarkersforVarRatio=opt$includeNonautoMarkersforVarRatio,
              sexCol=opt$sexCol,
              FemaleCode=opt$FemaleCode,
              FemaleOnly=opt$FemaleOnly,
              MaleCode=opt$MaleCode,
              MaleOnly=opt$MaleOnly,
              SampleIDIncludeFile=opt$SampleIDIncludeFile,
              VmatFilelist=opt$VmatFilelist,
              VmatSampleFilelist=opt$VmatSampleFilelist,
              longlCol=opt$longlCol,
              useGRMtoFitNULL=opt$useGRMtoFitNULL,
              offsetCol=opt$offsetCol,
              varWeightsCol=opt$varWeightsCol,
              sampleCovarCol=scovars,
              isStoreSigma=opt$isStoreSigma,
              isShrinkModelOutput=opt$isShrinkModelOutput,
              eCovarCol=ecovars
          )
          TRUE
        }, error = function(e) {
          cat("PCG refit failed:", e$message, "\n")
          FALSE
        })
        set_use_PCG(FALSE)

        pcg_resolved <- FALSE
        if(pcg_success){
          my_env_pcg = new.env()
          load(paste0(outputprefixsubset, ".pcg.rda"), envir = my_env_pcg)
          modglmm_pcg = my_env_pcg$modglmm
          print(modglmm_pcg$theta)
          if(!is.null(modglmm_pcg$theta) &&
             !(sum(modglmm_pcg$theta[2:length(modglmm_pcg$theta)]) <= 0 ||
               sum(modglmm_pcg$theta[2:length(modglmm_pcg$theta)]) > 1)){
            file.rename(paste0(outputprefixsubset, ".pcg.rda"), paste0(outputprefixsubset, ".rda"))
            tauVec = modglmm_pcg$theta
            tauTable = rbind(tauTable, tauVec)
            cat("PCG solver succeeded, theta in bounds\n")
            pcg_resolved <- TRUE
          } else {
            cat("PCG solver: theta still out of bounds\n")
            if(file.exists(paste0(outputprefixsubset, ".pcg.rda"))) file.remove(paste0(outputprefixsubset, ".pcg.rda"))
          }
        }

        if(!pcg_resolved){
          cat("All variance component parameter estiamtes are out of bounds, now try including all covariates as offset\n")
          opt$isCovariateOffset = TRUE
          set.seed(1)
          fitNULLGLMM_multiV(plinkFile=opt$plinkFile,
              bedFile=opt$bedFile,
              bimFile=opt$bimFile,
              famFile=opt$famFile,
              useSparseGRMtoFitNULL=opt$useSparseGRMtoFitNULL,
              sparseGRMFile=opt$sparseGRMFile,
              sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
              phenoFile = opt$phenoFile,
              phenoCol = opt$phenoCol,
              isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
              sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
              traitType = opt$traitType,
              outputPrefix = paste0(outputprefixsubset, ".offset"),
              isCovariateOffset=opt$isCovariateOffset,
              nThreads = opt$nThreads,
              useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
              invNormalize = opt$invNormalize,
              covarColList = covars,
              qCovarCol = qcovars,
              tol=opt$tol,
              maxiter=opt$maxiter,
              tolPCG=opt$tolPCG,
              maxiterPCG=opt$maxiterPCG,
              SPAcutoff = opt$SPAcutoff,
              numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
              skipModelFitting = opt$skipModelFitting,
              skipVarianceRatioEstimation = TRUE,
              memoryChunk = opt$memoryChunk,
              tauInit = tauInit,
              LOCO = opt$LOCO,
              isLowMemLOCO = opt$isLowMemLOCO,
              traceCVcutoff = opt$traceCVcutoff,
              nrun = opt$nrun,
              ratioCVcutoff = opt$ratioCVcutoff,
              outputPrefix_varRatio = opt$outputPrefix_varRatio,
              IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
              relatednessCutoff = opt$relatednessCutoff,
              isCateVarianceRatio = opt$isCateVarianceRatio,
              cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
              cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
              isCovariateTransform = opt$isCovariateTransform,
              isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
              minMAFforGRM = opt$minMAFforGRM,
              maxMissingRateforGRM = opt$maxMissingRateforGRM,
              minCovariateCount=opt$minCovariateCount,
              includeNonautoMarkersforVarRatio=opt$includeNonautoMarkersforVarRatio,
              sexCol=opt$sexCol,
              FemaleCode=opt$FemaleCode,
              FemaleOnly=opt$FemaleOnly,
              MaleCode=opt$MaleCode,
              MaleOnly=opt$MaleOnly,
              SampleIDIncludeFile=opt$SampleIDIncludeFile,
              VmatFilelist=opt$VmatFilelist,
              VmatSampleFilelist=opt$VmatSampleFilelist,
              longlCol=opt$longlCol,
              useGRMtoFitNULL=opt$useGRMtoFitNULL,
              offsetCol=opt$offsetCol,
              varWeightsCol=opt$varWeightsCol,
              sampleCovarCol=scovars,
              isStoreSigma=opt$isStoreSigma,
              isShrinkModelOutput=opt$isShrinkModelOutput,
              eCovarCol=ecovars
          )
          my_env = new.env()
          load(paste0(outputprefixsubset, ".offset.rda"), envir = my_env)
          modglmm = my_env$modglmm
          print(modglmm$theta)
          if(sum(modglmm$theta[2:length(modglmm$theta)]) <= 0 || sum(modglmm$theta[2:length(modglmm$theta)]) > 1){
            cat("All variance component parameter estiamtes are out of bounds.\n")
            file.remove(paste0(outputprefixsubset, ".offset.rda"))
            if (file.exists(paste0(outputprefixsubset, ".offset.varianceRatio.txt"))) {
              file.remove(paste0(outputprefixsubset, ".offset.varianceRatio.txt"))
            }else{
              if (file.exists(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))) {
                file.remove(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))
              }
            }
          }else{
            tauVec = modglmm$theta
            tauTable = rbind(tauTable, tauVec)
            file.rename(paste0(outputprefixsubset, ".offset.rda"), paste0(outputprefixsubset, ".rda"))
            if (file.exists(paste0(outputprefixsubset, ".offset.varianceRatio.txt"))) {
              file.rename(paste0(outputprefixsubset, ".offset.varianceRatio.txt"), paste0(outputprefixsubset, ".varianceRatio.txt"))
            }else{
              if (file.exists(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))) {
                file.rename(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"), paste0(opt$outputPrefix_varRatio, ".varianceRatio.txt"))
              }
            }
          }
        }


  }else{
  	tauVec = modglmm$theta
	tauTable = rbind(tauTable, tauVec)
  }
}else{

  if(!(sum(modglmm$theta[2:length(modglmm$theta)]) <= 0 || sum(modglmm$theta[2:length(modglmm$theta)]) > 1)){
	tauVec = modglmm$theta
	tauTable = rbind(tauTable, tauVec)
  }
}


file.remove(paste0(outputprefixsubset, ".rda"))
file.remove(datasubsetFile)


     }
}



if(!is.null(tauTable)){
tauInit = colMeans(tauTable)
#opt$maxiter = 1

print("tauInit")
print(tauInit)
}else{
tauInit <- convertoNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")
}


# Tracking variables for the final report
step1_used_pcg <- FALSE
step1_used_offset <- opt$isCovariateOffset

fit_success <- TRUE
tryCatch({

#set seed
set.seed(1)
fitNULLGLMM_multiV(plinkFile=opt$plinkFile,
	    bedFile=opt$bedFile,
	    bimFile=opt$bimFile,
	    famFile=opt$famFile,
	    useSparseGRMtoFitNULL=opt$useSparseGRMtoFitNULL, 
            sparseGRMFile=opt$sparseGRMFile,
            sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
            phenoFile = opt$phenoFile,
            phenoCol = opt$phenoCol,
	    isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
            sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
            traitType = opt$traitType,
            outputPrefix = opt$outputPrefix,
	    isCovariateOffset=opt$isCovariateOffset,
            nThreads = opt$nThreads,
	    useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
            invNormalize = opt$invNormalize,
            covarColList = covars,
            qCovarCol = qcovars,
	    tol=opt$tol,
	    maxiter=opt$maxiter,
            tolPCG=opt$tolPCG,
            maxiterPCG=opt$maxiterPCG,
            SPAcutoff = opt$SPAcutoff,
            numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
            skipModelFitting = opt$skipModelFitting,
	    skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
            memoryChunk = opt$memoryChunk,
            tauInit = tauInit,
            LOCO = opt$LOCO,
	    isLowMemLOCO = opt$isLowMemLOCO,
            traceCVcutoff = opt$traceCVcutoff,
	    nrun = opt$nrun,
            ratioCVcutoff = opt$ratioCVcutoff,
	    outputPrefix_varRatio = opt$outputPrefix_varRatio,
	    IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
            relatednessCutoff = opt$relatednessCutoff,
            isCateVarianceRatio = opt$isCateVarianceRatio,
            cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
            cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
            isCovariateTransform = opt$isCovariateTransform,
            isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
	    minMAFforGRM = opt$minMAFforGRM,
	    maxMissingRateforGRM = opt$maxMissingRateforGRM,
	    minCovariateCount=opt$minCovariateCount,
	    includeNonautoMarkersforVarRatio=opt$includeNonautoMarkersforVarRatio,
	    sexCol=opt$sexCol,
    	    FemaleCode=opt$FemaleCode,
	    FemaleOnly=opt$FemaleOnly,
	    MaleCode=opt$MaleCode,
	    MaleOnly=opt$MaleOnly,
	    SampleIDIncludeFile=opt$SampleIDIncludeFile,
	    VmatFilelist=opt$VmatFilelist,
	    VmatSampleFilelist=opt$VmatSampleFilelist,
	    longlCol=opt$longlCol,
	    useGRMtoFitNULL=opt$useGRMtoFitNULL,
	    offsetCol=opt$offsetCol,
	    varWeightsCol=opt$varWeightsCol,
	    sampleCovarCol=scovars,
	    isStoreSigma=opt$isStoreSigma,
      isShrinkModelOutput=opt$isShrinkModelOutput,
      eCovarCol=ecovars
	)

}, error = function(e) {
  message("Initial model failed with error: ", e$message)
    fit_success <<- FALSE  # Track failure
})


# ── Helper functions for the 4-stage fallback ──────────────────────────────
.fitStage <- function(suffix, isCovOff, usePCG) {
  if (usePCG) set_use_PCG(TRUE)
  set.seed(1)
  result <- tryCatch(
    fitNULLGLMM_multiV(
      plinkFile=opt$plinkFile, bedFile=opt$bedFile,
      bimFile=opt$bimFile, famFile=opt$famFile,
      useSparseGRMtoFitNULL=opt$useSparseGRMtoFitNULL,
      sparseGRMFile=opt$sparseGRMFile,
      sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
      phenoFile=opt$phenoFile, phenoCol=opt$phenoCol,
      isRemoveZerosinPheno=opt$isRemoveZerosinPheno,
      sampleIDColinphenoFile=opt$sampleIDColinphenoFile,
      traitType=opt$traitType,
      outputPrefix=paste0(opt$outputPrefix, suffix),
      isCovariateOffset=isCovOff,
      nThreads=opt$nThreads,
      useSparseGRMforVarRatio=opt$useSparseGRMforVarRatio,
      invNormalize=opt$invNormalize,
      covarColList=covars, qCovarCol=qcovars,
      tol=opt$tol, maxiter=opt$maxiter,
      tolPCG=opt$tolPCG, maxiterPCG=opt$maxiterPCG,
      SPAcutoff=opt$SPAcutoff,
      numMarkersForVarRatio=opt$numRandomMarkerforVarianceRatio,
      skipModelFitting=opt$skipModelFitting,
      skipVarianceRatioEstimation=opt$skipVarianceRatioEstimation,
      memoryChunk=opt$memoryChunk, tauInit=tauInit,
      LOCO=opt$LOCO, isLowMemLOCO=opt$isLowMemLOCO,
      traceCVcutoff=opt$traceCVcutoff, nrun=opt$nrun,
      ratioCVcutoff=opt$ratioCVcutoff,
      outputPrefix_varRatio=opt$outputPrefix_varRatio,
      IsOverwriteVarianceRatioFile=opt$IsOverwriteVarianceRatioFile,
      relatednessCutoff=opt$relatednessCutoff,
      isCateVarianceRatio=opt$isCateVarianceRatio,
      cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude,
      cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude,
      isCovariateTransform=opt$isCovariateTransform,
      isDiagofKinSetAsOne=opt$isDiagofKinSetAsOne,
      minMAFforGRM=opt$minMAFforGRM,
      maxMissingRateforGRM=opt$maxMissingRateforGRM,
      minCovariateCount=opt$minCovariateCount,
      includeNonautoMarkersforVarRatio=opt$includeNonautoMarkersforVarRatio,
      sexCol=opt$sexCol, FemaleCode=opt$FemaleCode,
      FemaleOnly=opt$FemaleOnly, MaleCode=opt$MaleCode, MaleOnly=opt$MaleOnly,
      SampleIDIncludeFile=opt$SampleIDIncludeFile,
      VmatFilelist=opt$VmatFilelist, VmatSampleFilelist=opt$VmatSampleFilelist,
      longlCol=opt$longlCol, useGRMtoFitNULL=opt$useGRMtoFitNULL,
      offsetCol=opt$offsetCol, varWeightsCol=opt$varWeightsCol,
      sampleCovarCol=scovars, isStoreSigma=opt$isStoreSigma,
      isShrinkModelOutput=opt$isShrinkModelOutput, eCovarCol=ecovars
    ),
    error = function(e) { cat("Fit failed:", e$message, "\n"); NULL }
  )
  if (usePCG) set_use_PCG(FALSE)
  !is.null(result)
}

.thetaOK <- function(suffix) {
  rda <- paste0(opt$outputPrefix, suffix, ".rda")
  if (!file.exists(rda) || file.size(rda) == 0) return(FALSE)
  e <- new.env()
  ok <- tryCatch({ load(rda, envir = e); TRUE }, error = function(e) FALSE)
  if (!ok) return(FALSE)
  m <- e$modglmm; print(m$theta)
  !is.null(m$theta) &&
    !(sum(m$theta[2:length(m$theta)]) <= 0 || sum(m$theta[2:length(m$theta)]) > 1)
}

.promote <- function(suffix) {
  from_rda <- paste0(opt$outputPrefix, suffix, ".rda")
  to_rda   <- paste0(opt$outputPrefix, ".rda")
  if (file.exists(from_rda)) file.rename(from_rda, to_rda)
  from_vr <- paste0(opt$outputPrefix, suffix, ".varianceRatio.txt")
  to_vr   <- paste0(opt$outputPrefix, ".varianceRatio.txt")
  if (file.exists(from_vr)) {
    file.rename(from_vr, to_vr)
  } else {
    from_vr2 <- paste0(opt$outputPrefix_varRatio, suffix, ".varianceRatio.txt")
    to_vr2   <- paste0(opt$outputPrefix_varRatio, ".varianceRatio.txt")
    if (file.exists(from_vr2)) file.rename(from_vr2, to_vr2)
  }
}

.cleanup <- function(suffix) {
  for (f in c(paste0(opt$outputPrefix, suffix, ".rda"),
              paste0(opt$outputPrefix, suffix, ".varianceRatio.txt"),
              paste0(opt$outputPrefix_varRatio, suffix, ".varianceRatio.txt")))
    if (file.exists(f)) file.remove(f)
}

# ── 4-stage fallback ─────────────────────────────────────────────────────────
# Stage 1: SMW + isCovariateOffset=FALSE  (initial fit, already done)
# Stage 2: SMW + isCovariateOffset=TRUE
# Stage 3: PCG + isCovariateOffset=FALSE
# Stage 4: PCG + isCovariateOffset=TRUE   (final, kept regardless)

if(fit_success){

  if(!opt$isCovariateOffset){
    # Stage 1 result is in .rda — check theta
    if(!.thetaOK("")){
      # Stage 2: SMW + offset=TRUE
      cat("Stage 1 theta OOB. Trying Stage 2: SMW + isCovariateOffset=TRUE\n")
      .fitStage(".stage2", isCovOff=TRUE, usePCG=FALSE)
      if(.thetaOK(".stage2")){
        .promote(".stage2")
        step1_used_offset <- TRUE
        cat("Stage 2 succeeded\n")
      } else {
        .cleanup(".stage2")
        # Stage 3: PCG + offset=FALSE
        cat("Stage 2 theta OOB. Trying Stage 3: PCG + isCovariateOffset=FALSE\n")
        .fitStage(".stage3", isCovOff=FALSE, usePCG=TRUE)
        if(.thetaOK(".stage3")){
          .promote(".stage3")
          step1_used_pcg <- TRUE
          cat("Stage 3 succeeded\n")
        } else {
          .cleanup(".stage3")
          # Stage 4: PCG + offset=TRUE (final — keep regardless)
          cat("Stage 3 theta OOB. Trying Stage 4: PCG + isCovariateOffset=TRUE\n")
          .fitStage(".stage4", isCovOff=TRUE, usePCG=TRUE)
          .promote(".stage4")
          step1_used_pcg    <- TRUE
          step1_used_offset <- TRUE
          cat("Stage 4 complete (final fallback)\n")
        }
      }
    }
    # else: Stage 1 theta OK, .rda kept as-is
  }

}else{
  # fit_success=FALSE: initial fit errored — start at Stage 2
  # Stage 2: SMW + offset=TRUE
  cat("Initial fit errored. Trying Stage 2: SMW + isCovariateOffset=TRUE\n")
  .fitStage(".stage2", isCovOff=TRUE, usePCG=FALSE)
  if(.thetaOK(".stage2")){
    .promote(".stage2")
    step1_used_offset <- TRUE
    cat("Stage 2 succeeded\n")
  } else {
    .cleanup(".stage2")
    # Stage 3: PCG + offset=FALSE
    cat("Stage 2 theta OOB. Trying Stage 3: PCG + isCovariateOffset=FALSE\n")
    .fitStage(".stage3", isCovOff=FALSE, usePCG=TRUE)
    if(.thetaOK(".stage3")){
      .promote(".stage3")
      step1_used_pcg <- TRUE
      cat("Stage 3 succeeded\n")
    } else {
      .cleanup(".stage3")
      # Stage 4: PCG + offset=TRUE (final — keep regardless)
      cat("Stage 3 theta OOB. Trying Stage 4: PCG + isCovariateOffset=TRUE\n")
      .fitStage(".stage4", isCovOff=TRUE, usePCG=TRUE)
      .promote(".stage4")
      step1_used_pcg    <- TRUE
      step1_used_offset <- TRUE
      cat("Stage 4 complete (final fallback)\n")
    }
  }
}


# ── Report ──────────────────────────────────────────────────────────────────

write_step1_report <- function(outputPrefix, outputPrefix_varRatio,
                               used_pcg, used_offset, fit_success, opt) {
  dir.create("./report", showWarnings = FALSE, recursive = TRUE)
  report_file <- file.path("./report",
                           paste0(basename(outputPrefix), ".step1_report.txt"))

  model_file <- paste0(outputPrefix, ".rda")

  # Find variance ratio file
  var_ratio_file <- ""
  for (f in c(paste0(outputPrefix_varRatio, ".varianceRatio.txt"),
              paste0(outputPrefix, ".varianceRatio.txt"))) {
    if (file.exists(f)) { var_ratio_file <- f; break }
  }

  # Load final model to read theta and convergence
  theta_vals          <- NA
  converged           <- FALSE
  is_cov_offset_final <- used_offset
  if (file.exists(model_file)) {
    tmp_env <- new.env()
    load(model_file, envir = tmp_env)
    m <- tmp_env$modglmm
    theta_vals <- m$theta
    converged  <- isTRUE(m$converged)
    if (!is.null(m$isCovariateOffset)) is_cov_offset_final <- m$isCovariateOffset
  }

  theta_oob <- is.null(theta_vals) || is.na(theta_vals[1]) ||
               sum(theta_vals[2:length(theta_vals)]) <= 0 ||
               sum(theta_vals[2:length(theta_vals)]) > 1

  overall_ok <- fit_success && file.exists(model_file) && converged && !theta_oob
  status_str <- if (overall_ok) "SUCCESS" else "FAILURE"
  solver_str <- if (used_pcg) "PCG" else "SMW"
  theta_str  <- if (all(!is.na(theta_vals))) {
    paste(round(theta_vals, 6), collapse = ", ")
  } else { "NA" }

  writeLines(c(
    "=== SAIGE-QTL Step 1 Fitting Report ===",
    paste0("Date              : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Output prefix     : ", outputPrefix),
    paste0("Phenotype         : ", opt$phenoCol),
    paste0("Trait type        : ", opt$traitType),
    "",
    "--- Convergence ---",
    paste0("Status            : ", status_str),
    paste0("fit_success       : ", fit_success),
    paste0("converged flag    : ", converged),
    paste0("Theta in bounds   : ", !theta_oob),
    paste0("Final theta       : ", theta_str),
    "",
    "--- Solver ---",
    paste0("Solver used       : ", solver_str),
    paste0("isCovariateOffset : ", is_cov_offset_final),
    "",
    "--- Output files ---",
    paste0("Model file        : ", model_file,
           " (exists: ", file.exists(model_file), ")"),
    paste0("Var ratio file    : ",
           if (var_ratio_file != "") var_ratio_file else "not found")
  ), report_file)

  cat("\n=== FINAL SUMMARY ===\n")
  cat(sprintf("Analysis Status  : %s\n", status_str))
  cat(sprintf("Solver used      : %s\n", solver_str))
  cat(sprintf("isCovariateOffset: %s\n", is_cov_offset_final))
  if (overall_ok) {
    cat("Model file       :", model_file, "\n")
    if (var_ratio_file != "") cat("Var ratio file   :", var_ratio_file, "\n")
  } else {
    cat("Step 1 did not converge cleanly. Check report for details.\n")
  }
  cat(sprintf("Report saved     : %s\n", report_file))

  invisible(list(
    convergence_status  = status_str,
    final_model_file    = model_file,
    variance_ratio_file = var_ratio_file,
    solver              = solver_str,
    is_cov_offset       = is_cov_offset_final
  ))
}

if (opt$isWriteReport) {
  write_step1_report(
    outputPrefix          = opt$outputPrefix,
    outputPrefix_varRatio = opt$outputPrefix_varRatio,
    used_pcg              = step1_used_pcg,
    used_offset           = step1_used_offset,
    fit_success           = fit_success,
    opt                   = opt
  )
}
