#' @title Perform voxelwise genome-wide association study.
#'
#' @description
#' Loads and prepares the data, runs preliminary (if necessary) and full analyses.
#' Then visualises the results.
#'
#' @export
#' @import fslr
#' @import MatrixEQTL
#' @param genePath Path to a directory containing plink files with genomic data.
#' @param brainPath Path to a directory containing imaging data.
#' @param mdtPath Path to a referance image.
#' @param maskPath Path to an image mask.
#' @param infoPath Path to .csv containing subject demographics.
#' @param outPath Path to output directory.
#' @param subFactor Downsampling factor.
#' @param out.subFactor Output images downsampling factor. Default: \code{max(1, subFactor-1)}.
#' @param force.snps \code{character} vector of SNPs forced to be analysed even if not passing the quality control.
#' @param useModel Regression model to use.
#' @param top.thresh Number of top SNPs to be analysed and visualised.
#' @param log.cutoff Negative log p-value cutoff value for results visualisation.
#' @param eff.no.tests Effective number of tests (SNPs) for p-value correction.
#' @param sampleSize Subject sample size.
#' @param mockPath Path to input directory with saved objects for vGWAS mocking. Mock parameters should match  original parameters. \code{NULL} for real analyses.
performVGWAS <- function(genePath, brainPath, mdtPath, maskPath, infoPath, subFactor, out.subFactor = max(1, subFactor-1), outPath = NULL, force.snps = NULL, useModel = modelLINEAR, top.thresh = 5, log.cutoff = 4, eff.no.tests = 275575, sampleSize = NULL, mockPath = NULL){
  cat("Preparing data...\n")

  # Create output directory if needed
  if(!is.null(outPath) & !dir.exists(outPath)){
    dir.create(outPath)
  }

  # List imaging files
  niiFiles <- list.files(brainPath, full.names = TRUE)[]
  niiIDs <- substr(list.files(brainPath), 6, 15)

  # Load subject personal data
  subject.info <- read.csv(infoPath)

  # Select subjects from ADNI1, with their respective images available
  subject.info <- subject.info[which(subject.info$COLPROT == "ADNI1" & subject.info$PTID %in% niiIDs), c("PTID", "PTGENDER", "AGE", "PTETHCAT", "PTRACCAT")]
  subject.info <- unique(subject.info)
  rownames(subject.info) <- subject.info$PTID

  # Select non-hispanic/latino whites
  subject.info <- subject.info[which(subject.info$PTETHCAT=="Not Hisp/Latino" & subject.info$PTRACCAT=="White"),]
  covar <- subject.info[,c("PTGENDER", "AGE")]
  levels(covar$PTGENDER)[levels(covar$PTGENDER)=="Male"] <- 1
  levels(covar$PTGENDER)[levels(covar$PTGENDER)=="Female"] <- 2
  covar <- t(covar)
  covar <- covar[, order(colnames(covar))]



  ############################################################ study-spec ^^^^

  # Load genomic data
  plinkFiles <- list.files(genePath, full.names = TRUE)
  cat("Loading genomic data... ")
  time <- system.time(snps <- readSNPs(plinkFiles = plinkFiles, force.snps = force.snps))
  cat(sprintf("%.2fs\n", time[3]))
  no.snps <- nrow(snps@.Data)

  # List subjects
  subjects <- colnames(snps[, colnames(snps) %in% subject.info$PTID])
  niiFiles <- niiFiles[niiIDs %in% subjects]
  # Set sample
  sample <- 1:min(sampleSize,length(niiFiles))

  # Mask: whole brain
  cat("Loading image mask...\n")
  mask <- fslr::readNIfTI2(maskPath)
  if (subFactor > 0){
    mask <- downsample(mask, subFactor, bin = TRUE)
  }

  # Load imaging data
  cat("Loading imaging data... ")
  if(is.null(mockPath)){
    # Real loading
    time <- system.time(pheno <- readFlatROIs(paths =  niiFiles[sample], ids = subjects[sample], subFactor = subFactor, mask = mask))
    cat(sprintf("%.2fs\n", time[3]))
    if(!is.null(outPath))
      save(pheno, file = paste(outPath, "/flatROIs-", 2^subFactor, ".R" , sep = ""))
  }
  else{
    # Load saved data
    cat("\n")
    load(file = paste(mockPath, "/flatROIs-", 2^subFactor, ".R" , sep = ""))
  }

  # Create SlicedData objects
  voxelData <- MatrixEQTL::SlicedData$new(pheno[, subjects[sample]])
  voxelData$ResliceCombined(500)
  rm(pheno)

  # Create SlicedData objects
  snpData <- MatrixEQTL::SlicedData$new(snps@.Data[, subjects[sample]])
  snpData$ResliceCombined(500)
  rm(snps)

  # Create SlicedData objects
  cvrt <- MatrixEQTL::SlicedData$new(covar[, subjects[sample]])
  cvrt$ResliceCombined(500)


  cat("Running preliminary analysis...\n")
  # Perform preliminary analysis
  if(top.thresh > 0){
    if(is.null(mockPath)){
      # Real analysis
      meh <- vGWAS(snpData = snpData, voxelData = voxelData, cvrt = cvrt, useModel = useModel, prescan = TRUE)
      if(!is.null(outPath))
        save(meh, file = paste(outPath, "/pre-", 2^subFactor, ".R" , sep = ""))
    }
    else{
      # Load saved data
      load(file = paste(mockPath, "/pre-", 2^subFactor, ".R" , sep = ""))
    }
    # Top p-values
    top.snps <- names(sort(meh$all$min.pv.snps)[1:top.thresh])
  }
  else{
    top.snps <- c()
  }

  # Select snps of interest
  sel.snps <- unique(c(top.snps, force.snps))
  # Save selected p-values
  snp.min.pv <- meh$all$min.pv.snps[sel.snps]
  # Save p-values by voxel
  voxel.min.pv <- meh$all$min.pv.gene

  # Create SlicedData objects
  snpData.selected <- selectSNPs(sel.snps = sel.snps, snpData = snpData)
  snpData.selected <- SlicedData$new(snpData.selected)
  snpData.selected$ResliceCombined(500)

  # Run full analysis
  cat("Running full analysis...\n")
  meh <- vGWAS(snpData = snpData.selected, voxelData = voxelData, cvrt = cvrt, useModel = useModel, prescan = FALSE)

  # Get results for selected SNPs
  results.by.snp <- getSNPresults(meh, sel.snps)

  # Calculate FDR for the best SNP
  if(eff.no.tests == 0){
    eff.no.tests <- no.snps
  }
  p.value <- sort(voxel.min.pv)
  voxel.results <- getSNPfdr(p.value, eff.no.tests)

  # Visualise snp.results
  cat("Results visualisation...\n")
  mdt <- fslr::readNIfTI2(mdtPath)
  mdt <- fslr::fslthresh(mdt, verbose = FALSE)
  mdt <- downsample(mdt, out.subFactor)

  if(!is.null(outPath))
    pdf(paste(outPath, "/report-", 2^subFactor, ".pdf", sep = ""))
  hist(voxel.results[,"p.value"], main ="Histogram of raw p-values", xlab = "p-value")
  hist(voxel.results[,"pc.value"], main ="Histogram of pc-values", xlab = "pc-value")
  result.images.by.snp.cut <- lapply(X = sel.snps, results = results.by.snp, mask = mask, mdt = mdt, log.cutoff = log.cutoff, afterTitle = paste("\ncut-off at ", log.cutoff), outPath = outPath, FUN = visualiseSNP)
  result.images.by.snp <- lapply(X = sel.snps, results = results.by.snp, mask = mask, mdt = mdt, afterTitle = "\nno cut-off", outPath = outPath, FUN = visualiseSNP)
  names(result.images.by.snp.cut) <- sel.snps
  names(result.images.by.snp) <- sel.snps
  result.image.by.vox <- visualiseVox(voxel.min.pv, mask = mask, mdt = mdt, log.cutoff = log.cutoff, outPath = outPath)
  if(!is.null(outPath)){
    dev.off()
    writeNIfTI2(result.image.by.vox, paste(outPath, "/result-image-by-vox-", 2^subFactor, ".nii", sep = ""))
    lapply(X = sel.snps, FUN = function(x){writeNIfTI2(result.images.by.snp.cut[[x]], paste(outPath, "/result-image-", x, "-", 2^subFactor, ".nii", sep = ""))})
  }
}
