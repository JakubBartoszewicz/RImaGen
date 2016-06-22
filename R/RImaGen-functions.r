#' Get False Discovery Rates.
#'
#' @export
#' @import qvalue
#' @param p.value A vector of raw p-values.
#' @param eff.no.tests Efficient number of performed tests.
#' @return A summary of FDR q-values.
getSNPfdr <- function(p.value, eff.no.tests)
{
  # Calculate corrected p-values using Beta distribution
  pc.value <- pbeta(p.value, 1, eff.no.tests)
  # Calculate FDR
  FDR <- p.adjust(pc.value)
  # Calculate FDR q-values
  q.value.FDR.obj <- qvalue(pc.value, pfdr = FALSE)
  q.value.FDR <- q.value.FDR.obj$qvalues
  # Calculate pFDR q-values
  q.value.pFDR.obj <- qvalue(pc.value, pfdr = TRUE)
  q.value.pFDR <- q.value.pFDR.obj$qvalues
  # Bind into one data.frame
  voxel.results <- cbind(p.value, pc.value, FDR, q.value.FDR, q.value.pFDR)
  return (voxel.results)
}

#' Downsample an image.
#'
#' @export
#' @import fslr
#' @param img An image of class \code{nifti}.
#' @param subFactor Downsampling factor.
#' @param method Downsampling method.
#' @param bin \code{logical}. Set to \code{TRUE} for downsampling of binary images or masks.
#' @param thresh Treshhold for re-binarisation of downsampled binary images and masks. A threshold near 1 (say 0.9) is conservative according to FLIRT FAQ.
#' @return A downsampled image.
downsample <- function(img, subFactor = 1, method = c("FLIRT", "SUBSAMP2"), bin = FALSE, thresh = 0.9){
  # Ensure correct mathod
  type = match.arg(method)
  if (subFactor > 0){
    # If subFactor is correct
    if(type == "SUBSAMP2")
    {
      # Either use fslmaths -subsamp2 if you really must
      for (i in 1:subFactor){
        img <- fslr::fslsub2(img, verbose = TRUE)
      }
    }
    if(type == "FLIRT")
    {
      # Or use flirt (recommended)
      img <- fslr::flirt(img, img, opts = paste("-applyisoxfm ", 2^subFactor))
    }
    if(bin)
    {
      # If binary mask, rebinarise
      img <- binariseMask(img, thresh)
    }
  }
  return (img)
}

#' Get vGWAS results for selected SNPs.
#'
#' @export
#' @param meh vGWAS results.
#' @param sel.snps \code{character} vector of SNP names.
#' @return A named list of results for each SNP.
getSNPresults <- function(meh, sel.snps)
{
  # Split MatrixEQTL results by SNP name
  results.by.snp <- split(meh$all$eqtls, meh$all$eqtls$snps)
  # Select SNPs of interest
  results.by.snp <- results.by.snp[sel.snps]
  return (results.by.snp)
}

#' Visualise by-voxel statistics for a chosen SNP, e.g that SNP's p-value per voxel.
#'
#' @export
#' @param snp.name Chosen SNP name.
#' @param results A named list of results for each SNP.
#' @param mask Image mask (a \code{nifti} object or an array).
#' @param mdt Reference image.
#' @param log.cutoff Negative log p-value cutoff value.
#' @param plot \code{logical}. Set to \code{TRUE} for automated plotting of the results.
#' @param title Plot title.
#' @param titleSize Title font size.
#' @param beforeTitle Text to print before the title.
#' @param afterTitle Text to print after the title.
#' @return Results visualisation image.
visualiseSNP <- function(snp.name, results, mask, mdt, log.cutoff = 0, plot = TRUE, title = snp.name, titleSize = 1, beforeTitle = "-log10 p-values by\n", afterTitle =""){
  # Select SNP results by SNP name
  result <- results[[snp.name]][2:6]
  # Set row names to voxel numbers
  row.names(result) <- result$gene
  # Cast data types
  result <- result[as.character(1:nrow(result)),]
  # Visualise results
  img <- visualiseVox(result$pvalue, mask, mdt, log.cutoff, plot, title = title, titleSize = titleSize, beforeTitle = beforeTitle, afterTitle = afterTitle)
  return (img)
}

#' Visualise by-voxel statistics for each voxel, e.g winning SNP p-value per voxel.
#'
#' @export
#' @param result Vector of target statistics.
#' @param mask Image mask.
#' @param mdt Reference image.
#' @param log.cutoff Negative log p-value cutoff value.
#' @param plot \code{logical}. Set to \code{TRUE} for automated plotting of the results.
#' @param title Plot title.
#' @param titleSize Title font size.
#' @param beforeTitle Text to print before the title.
#' @param afterTitle Text to print after the title.
#' @return Results visualisation image.
visualiseVox <- function(result, mask, mdt, log.cutoff = 0, plot = TRUE, title = "winning SNP", titleSize = 1, beforeTitle = "-log10 p-values by\n", afterTitle =""){
  # Calculate log-negative results
  img <- -log10(result)
  # Find coordinates of voxels to be filled
  fit <- which(mask == 1)
  # Fill mask with results
  mask[fit] <- img
  # Upsample image
  mask <- fslr::flirt(mask, mask, opts = paste("-applyisoxfm ", mdt@pixdim[1]))
  # Create nifti image from the filled mask
  img <- fslr::niftiarr(img = mdt, arr = mask)
  # Cutoff low results
  img[which(img < log.cutoff)] <- 0
  if(plot){
    # Plot results with crosshairs centered at maximum voxel intensity
    fslr::ortho2(mdt, img, xyz = which(img == max(img), arr.ind = TRUE), text = paste(beforeTitle, title, afterTitle), text.cex = titleSize)
  }
  return (img)
}

#' Select SNPs by name from a \code{\linkS4class{SlicedData}} object.
#'
#' @export
#' @import MatrixEQTL
#' @param sel.snps \code{character} vector of SNP names.
#' @param snpData \code{\linkS4class{SlicedData}} object containing SNP data.
#' @return Array containing selected SNP data.
selectSNPs <- function(sel.snps, snpData){

  # Select rows by SNP name
  selected <- lapply(sel.snps, function(x){as.vector(snpData$FindRow(x)$row)})
  if(is.null(selected))
    warning("No SNPs selected!")
  # Transform selected data into an array
  selected <- t(simplify2array(selected))
  # Use original column names
  colnames(selected) <- snpData$columnNames
  # Name rows by SNP names
  rownames(selected) <- sel.snps

  return(selected)
}

#' @title Perform voxelwise genome-wide association analysis.
#'
#' @description
#' Checks for associations prepared \code{\linkS4class{SlicedData}} objects containing
#' flattened image ROIs, genomic data and covariates. Wrapper for \code{\link[MatrixEQTL]{Matrix_eQTL_engine}}.
#'
#' @export
#' @import MatrixEQTL
#' @param snpData \code{\linkS4class{SlicedData}} object containing selected SNPs.
#' @param voxelData \code{\linkS4class{SlicedData}} object containing flat voxel ROIs.
#' @param cvrt \code{\linkS4class{SlicedData}} objact containing covariates.
#' @param useModel Regression model to use.
#' @param prescan \code{logical}. \code{TRUE} for preliminary analysis allowing "winning SNP" ranking.
#' @param errorCovariance \code{numeric}. The error covariance matrix. Use \code{numeric()} for homoskedastic independent errors.
#' @param pvalue.hist \code{logical}, \code{numeric} or \code{"qq-plot"}.
#' Defines whether and how the distribution of p-values is recorded in the returned object.
#' Set to \code{FALSE} for a faster analysis. To record information for a histogram set to the number of bins.
#' To record information for a qq-plot, set to \code{qq-plot}.
#' Use \code{plot} on the returned object.
#' @return vGWAS results.
vGWAS <- function(snpData, voxelData, cvrt, useModel, prescan, errorCovariance = numeric(), pvalue.hist = FALSE){

  if(prescan){
    # Do not record full results for prescan
    p.val <- 1e-100
  }
  else{
    # Record all results for full analysis (even the non-significant ones)
    p.val <- 1
  }

  # Set temp outfile
  tmp = tempfile()

  # Run MatrixEQTL analysis
  meh <- MatrixEQTL::Matrix_eQTL_engine(
    snps = snpData,
    gene = voxelData,
    cvrt = cvrt,
    output_file_name = tmp,
    pvOutputThreshold = p.val,
    useModel = useModel,
    verbose = TRUE,
    pvalue.hist = pvalue.hist,
    min.pv.by.genesnp = prescan,
    noFDRsaveMemory = prescan)

  #unlink temp outfile
  unlink(tmp)
  return(meh)
}

#' Read SNPs performing genome quality control.
#'
#' @export
#' @import snpStats
#' @param plinkFiles Paths to the plink files
#' @param call.rate.cutoff Call rate cutoff value
#' @param maf.cutoff Minor allele frequency cutoff value
#' @param hwe.pval Hardy-Weinberg equilibrium p-value cutoff value
#' @param force.snps \code{character} vector of SNPs forced to be analysed even if not passing the quality control.
#' @param forcedOnly \code{logical}. If \code{TRUE}, only the SNPs of interest will be analysed, regardless of their data quality.
#' @return \code{\linkS4class{SnpMatrix}} of SNPs passing the quality control
readSNPs <- function(plinkFiles, call.rate.cutoff = 0.95, maf.cutoff = 0.1, hwe.pval = 5.7e-07, force.snps = NULL, forcedOnly = FALSE){
  # Read plink files
  genomes <- snpStats::read.plink(plinkFiles[1], plinkFiles[2], plinkFiles[3])

  if(forcedOnly){
    # If forcedOnly select only the chosen SNPs
    snps <- genomes$genotypes[,force.snps]
  }
  else{
    # Else, select all SNPs
    snps <- genomes$genotypes

    if(call.rate.cutoff > 0 | maf.cutoff > 0 | hwe.pval > 0){
      # Perform quality control using provided cutoff values
      csum <- snpStats::col.summary(snps)
      snps <- snps [,which(csum$Call.rate>=call.rate.cutoff & csum$MAF>=maf.cutoff & 2*pnorm(-abs(csum$z.HWE))>hwe.pval)]

      # Check if forced SNPs passed the QC
      forced <- force.snps[!(force.snps %in% colnames(snps))]

      if(length(forced > 0)){
        # If some of the forced SNPs did not pass the QC, add them to the data set
        warning(paste("Forcing SNPs not passing the quality control:\n", paste(capture.output(csum[forced,c("Call.rate", "MAF", "z.HWE")]), collapse = "\n")))
        snps <- cbind(snps, genomes$genotypes[,forced])
      }
    }
  }
  # Transpose snps matrix
  snps <- t(snps)
  # Set column names to subject IDs
  colnames(snps) <- genomes$fam[,2]
  # Sort columns by subject ID
  snps <- snps[, order(colnames(snps))]
  rm(genomes)
  gc()
  return(snps)
}

#' Read flattened ROIs from a list of NIfTI files using a cluster of choice.
#'
#' @export
#' @import parallel
#' @param paths A list of paths to the NIfTI files.
#' @param ids A list of subject IDs.
#' @param subFactor Downsampling factor.
#' @param mask Image mask (a \code{nifti} object or an array).
#' @param method Downsampling method.
#' @param nCores number of cores to use. Default: all physical cores available
#' @param clType cluster type. Default: \code{PSOCK}
#' @param makeCluster \code{logical}. \code{TRUE} if a cluster should be made (default).
#' @param userCluster Cluster to use specified by the user. If \code{NULL} and makeCluster is \code{FALSE}, a default cluster will be used.
#' @return Array of chosen voxel intensities for all the subjects.
readFlatROIs <- function(paths, ids, subFactor = 0, mask = NULL, method = c("FLIRT", "SUBSAMP2"), nCores = getOption("mc.cores", detectCores(logical = FALSE)), clType = "PSOCK", makeCluster = TRUE, userCluster = NULL)
{
  # Get the coordinates of all "on" voxels in the mask resulting in a list of considered voxels
  voxels <- which(mask != 0, arr.ind = T)
  if (makeCluster){
    # Make cluster explicitly (assuming all the data stored locally anyway)
    cl <- parallel::makeCluster(nCores, clType)
    parallel::clusterExport(cl, c("readNIfTI2", "fslsub2", "flirt", "downsample"))
  }
  else{
    # Use cluster provided by user
    cl <- userCluster
  }
  # Load flat ROIs
  rois <- parallel::parLapply(cl, paths, readFlatROI, subFactor, voxels, method)
  # Stop cluster
  parallel::stopCluster(cl)
  # Transform ROIs to an array
  rois <- simplify2array(rois)
  # Set column names to subject IDs
  colnames(rois) <- ids
  # Set row names to voxel numbers
  rownames(rois) <- 1:dim(voxels)[1]
  return(rois)
}
