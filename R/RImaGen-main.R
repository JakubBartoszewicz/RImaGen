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
#' @param niiFiles Paths to files containing imaging data.
#' @param niiIDs Subject IDs for each image. Shuffle for permutation tests.
#' @param ref.imgPath Path to a referance image.
#' @param maskPath Path to an image mask.
#' @param infoPath Path to .csv containing subject demographics.
#' @param subFactor Downsampling factor.
#' @param covar Covariates matrix with subject IDs as column names.
#' @param errorCovariance Covariance matrix for the error term. Set to \code{numeric()} for multiple of identity (default, most cases).
#' @param outPath Path to output directory.
#' @param out.subFactor Output images downsampling factor. Default: equal to \code{subFactor}.
#' @param matPath Path to convolution matrix for coregistration of results to the reference image. Set to \code{NULL} if in the same space.
#' @param force.snps \code{character} vector of SNPs forced to be analysed even if not passing the quality control.
#' @param useModel Regression model to use.
#' @param top.thresh Number of top SNPs to be analysed and visualised.
#' @param log.cutoff Negative log p-value cutoff value for results visualisation.
#' @param eff.no.tests Effective number of tests (SNPs) for p-value correction.
#' @param sampleSize Subject sample size.
#' @param randomSample \code{logical}. Set to \code{TRUE} for subject sample randomisation.
#' @param mockPath.flatROIs Path to input file with saved flatROIs object for vGWAS mocking. Mock parameters should match  original parameters.
#' Useful for registering results to different standard spaces and resolutions, selecting alternative SNPs for the full analysys, of using different genetic models.
#' @param mockPath.pre Path to input file with saved preliminary analysis results object for vGWAS mocking. Mock parameters should match  original parameters.
#' Useful for registering results to different standard spaces and resolutions, selecting alternative SNPs for the full analysys, of using different genetic models.
#' \code{NULL} for real analyses.
#' @return A list containing all the resulting statistical parametric maps.
performVGWAS <- function(genePath, niiFiles, niiIDs, ref.imgPath, maskPath, infoPath, subFactor,
                         covar = character(), errorCovariance = numeric(), out.subFactor = subFactor, matPath = NULL,
                         outPath = NULL, force.snps = NULL, useModel = MatrixEQTL::modelLINEAR, top.thresh = 5, log.cutoff = 4,
                         eff.no.tests = 275575, sampleSize = NULL, randomSample = FALSE, mockPath.flatROIs = NULL, mockPath.pre = NULL){

  cat("Preparing data...\n")

  # Create output directories if needed
  if(!is.null(outPath)){
    # Create output directory if needed
    if(!dir.exists(outPath)){
      dir.create(outPath)
    }
    # Create png output directory if needed
    if(!dir.exists(paste(outPath, "/png", sep = ""))){
      dir.create(paste(outPath, "/png", sep = ""))
    }
    # Create pdf output directory if needed
    if(!dir.exists(paste(outPath, "/pdf", sep = ""))){
      dir.create(paste(outPath, "/pdf", sep = ""))
    }
  }

# Load data ---------------------------------------------------------------


  # Load genomic data
  plinkFiles <- list.files(genePath, full.names = TRUE)
  cat("Loading genomic data... ")
  time <- system.time(snps <- readSNPs(plinkFiles = plinkFiles, force.snps = force.snps))
  cat(sprintf("%.2fs\n", time[3]))
  # Count SNPs
  no.snps <- nrow(snps@.Data)

  # List subjects
  subjects <- colnames(snps[, colnames(snps) %in% colnames(covar)])
  niiFiles <- niiFiles[niiIDs %in% subjects]

  # Set sample
  if(randomSample){
    # Randomize sample
    sub.sample <- sort(sample(x = 1:length(niiFiles), size = min(sampleSize,length(niiFiles))))
  }
  else{
    # Select first few
    sub.sample <- 1:min(sampleSize,length(niiFiles))
  }
  if(!is.null(outPath)){
    # Write sampled subjects info
    write.csv(t(covar[, subjects[sub.sample]]), file = paste(outPath, "/subjects.csv", sep = ""))
  }

  # Create SlicedData objects
  snpData <- MatrixEQTL::SlicedData$new(snps@.Data[, subjects[sub.sample]])
  snpData$ResliceCombined(500)
  rm(snps)
  gc()
  # Recode plink values
  for(i in 1:snpData$nSlices()){
    snpData[[i]][which(snpData[[i]] == 0)] <- NA
    snpData[[i]] <- snpData[[i]] - 1
  }

  # Mask: whole brain
  cat("Loading image mask...\n")
  mask <- fslr::readNIfTI2(maskPath)
  if (subFactor > 0){
    mask <- downsample(mask, subFactor, bin = TRUE)
  }

  # Load imaging data
  cat("Loading imaging data... ")
  if(is.null(mockPath.flatROIs)){
    # Real loading
    time <- system.time(pheno <- readFlatROIs(paths =  niiFiles[sub.sample], ids = subjects[sub.sample], subFactor = subFactor, mask = mask))
    cat(sprintf("%.2fs\n", time[3]))
    if(!is.null(outPath))
      save(pheno, file = paste(outPath, "/flatROIs-", 2^subFactor, ".R" , sep = ""))
  }
  else{
    # Load saved data
    cat("\n")
    load(file = mockPath.flatROIs)
  }

  # Create SlicedData objects
  voxelData <- MatrixEQTL::SlicedData$new(pheno[, subjects[sub.sample]])
  voxelData$ResliceCombined(500)
  rm(pheno)
  gc()


  # Create SlicedData objects
  cvrt <- MatrixEQTL::SlicedData$new(covar[, subjects[sub.sample]])
  cvrt$ResliceCombined(500)

# Run analysis ---------------------------------------------------------------

  cat("Running preliminary analysis...\n")
  # Perform preliminary analysis
  if(top.thresh > 0){
    if(is.null(mockPath.pre)){
      # Real analysis
      meh <- vGWAS(snpData = snpData, voxelData = voxelData, cvrt = cvrt, errorCovariance = errorCovariance, useModel = useModel, prescan = TRUE)
      if(!is.null(outPath))
        save(meh, file = paste(outPath, "/pre-", 2^subFactor, ".R" , sep = ""))
    }
    else{
      # Load saved data
      load(file = mockPath.pre)
    }
    # Top p-values
    min.pv.snps <- sort(meh$all$min.pv.snps)
    ranks.snps <- names(min.pv.snps)
    top.snps <- ranks.snps[1:top.thresh]
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
  rm(snpData)
  gc()

  # Run full analysis
  cat("Running full analysis...\n")
  meh <- vGWAS(snpData = snpData.selected, voxelData = voxelData, cvrt = cvrt, errorCovariance = errorCovariance, useModel = useModel, prescan = FALSE)
  if(top.thresh <= 0){
    min.pv.snps <- sort(meh$all$min.pv.snps)
    ranks.snps <- names(min.pv.snps)
  }

  # Get results for selected SNPs
  results.by.snp <- getSNPresults(meh, sel.snps)

  # Calculate FDR for the best SNP
  if(eff.no.tests == 0){
    eff.no.tests <- no.snps
  }
  p.value <- sort(voxel.min.pv)
  voxel.results <- getSNPfdr(p.value, eff.no.tests)

  # Extract ranks of selected SNPs
  ranks.sel.snps <- match(sel.snps, ranks.snps)
  # Extract min. p-values
  min.p.value <- min.pv.snps[sel.snps]
  # Extract test statistic
  statistic <- as.numeric(lapply(X = sel.snps, FUN = function(x){
    results.by.snp[[x]]$statistic[1]
  }))
  # Extract beta (slope of the regression line)
  if(useModel == MatrixEQTL::modelLINEAR){
  beta <- as.numeric(lapply(X = sel.snps, FUN = function(x){
    results.by.snp[[x]]$beta[1]
  }))}
  else{
    beta <- c()
  }
  # Prepare report
  report.snps <- cbind(ranks.sel.snps, min.p.value, statistic, beta)
  report.snps <- report.snps[order(report.snps[,"ranks.sel.snps"]),]

# Visualise results ---------------------------------------------------------------

  # Visualise SNPS results
  cat("Results visualisation...\n")
  ref.img <- fslr::readNIfTI2(ref.imgPath)
  ref.img <- fslr::fslthresh(ref.img, verbose = FALSE)
  ref.img <- downsample(ref.img, out.subFactor)

  if(!is.null(outPath)){
    write.csv(report.snps, file = paste(outPath, "/sel-snps-", 2^subFactor, ".csv", sep = ""))
    # Save histogram pdfs
    pdf(paste(outPath, "/pdf/hist-p-val-", 2^subFactor, ".pdf", sep = ""))
    hist(voxel.results[,"p.value"], main ="Histogram of winning SNP p-values", xlab = "p-value")
    dev.off()
    pdf(paste(outPath, "/pdf/hist-pc-val-", 2^subFactor, ".pdf", sep = ""))
    hist(voxel.results[,"pc.value"], main ="Histogram of winning SNP pc-values", xlab = "pc-value")
    dev.off()

    # Save histogram pngs
    png(filename = paste(outPath, "/png/hist-p-val-", 2^subFactor, ".pdf", sep = ""))
    hist(voxel.results[,"p.value"], main ="Histogram of winning SNP p-values", xlab = "p-value")
    dev.off()
    png(filename = paste(outPath, "/png/hist-pc-val-", 2^subFactor, ".pdf", sep = ""))
    hist(voxel.results[,"pc.value"], main ="Histogram of winning SNP pc-values", xlab = "pc-value")
    dev.off()

    # Open report file
    pdf(paste(outPath, "/report-", 2^subFactor, "-", ref.img@pixdim[2], "mm.pdf", sep = ""))
  }

  hist(voxel.results[,"p.value"], main ="Histogram of winning SNP p-values", xlab = "p-value")
  hist(voxel.results[,"pc.value"], main ="Histogram of winning SNP pc-values", xlab = "pc-value")

  result.images.by.snp.cut <- lapply(X = sel.snps, results = results.by.snp, mask = mask, ref.img = ref.img, log.cutoff = log.cutoff, afterTitle = paste("\ncut-off at", log.cutoff), outPath = outPath, matPath = matPath, FUN = visualiseSNP)
  result.images.by.snp <- lapply(X = sel.snps, results = results.by.snp, mask = mask, ref.img = ref.img, afterTitle = "\nno cut-off", outPath = outPath, matPath = matPath, FUN = visualiseSNP)
  names(result.images.by.snp.cut) <- sel.snps
  names(result.images.by.snp) <- sel.snps
  result.image.by.vox <- visualiseVox(voxel.min.pv, mask = mask, ref.img = ref.img, log.cutoff = log.cutoff, outPath = outPath, matPath = matPath)
  if(!is.null(outPath)){
    # Close report file
    dev.off()
    # Save .nii images
    writeNIfTI2(result.image.by.vox, paste(outPath, "/result-image-by-vox-", 2^subFactor, "-", ref.img@pixdim[2], "mm.nii", sep = ""))
    lapply(X = sel.snps, FUN = function(x){writeNIfTI2(result.images.by.snp.cut[[x]], paste(outPath, "/result-image-", x, "-", 2^subFactor, "-", ref.img@pixdim[2], "mm.nii", sep = ""))})
  }

  return (list(result.images.by.snp, result.images.by.snp.cut, result.image.by.vox))
}
