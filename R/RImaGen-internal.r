#' Read flattened ROI from a NIfTI file.
#'
#' @keywords internal
#' @param path A path to the NIfTI file.
#' @param subFactor Downsampling factor.
#' @param voxelList List of chosen voxel coordinates. If NULL, return the whole file.
#' @param method Downsampling method.
#' @return Vector of chosen voxel intensities.
readFlatROI <- function(path, subFactor = 0 , voxelList = NULL, method = c("FLIRT", "SUBSAMP2"))
{
  # Read NIfTI image
  img <- fslr::readNIfTI2(path)
  # Downsample
  img <- downsample(img, subFactor, method)
  if(length(voxelList)>0){
    # Apply voxelList mask
    roi <- img[voxelList]
    rm(img)
    return(roi)
  }
  else{
    # Or flatten the original image
    return(as.vector(img))
  }
}

#' Binarise mask.
#'
#' @keywords internal
#' @import fslr
#' @param mask An image of class \code{nifti}.
#' @param thresh Treshhold for binarisation of images and masks. A threshold near 1 (say 0.9) is conservative according to FLIRT FAQ.
#' @return A binarised image.
binariseMask <- function(mask, thresh = 0.9)
{
  # Threshold image
  mask <- fslr::fslthresh(mask, thresh = thresh)
  # Binarase image
  mask <- fslr::fslbin(mask)
  return(mask)
}
