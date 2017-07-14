#' @title Identify Target DEGs Associated with Each Significant Driver
#'
#' @description If a DEG is a target of a driver, the causal links between the driver and the DEG has to present in at least 20% (by default) of the tumors where the driver is present.
#'
#' @param driver identify the target DEGs for the given list of drivers
#' @param tumorset
#' @param TDIresults TDI Tumor-SGA-DEG triplets
#' @param postprobnoiselv posterior probability cutoff per SGA
#' @param numoftargenes.cutoff minimal number of target DEGs for an SGA to be called a driver in the tumor
#' @param DEGcallrate.thd minimal call rate for a DEG to be associated with a driver
#'
#' @export
#'
#' @return Driver.DEGlist
#'
#' @examples \dontrun{
#'
#' }
#'

identify.tarDEGlist.givendriverNtumorset <- function(driver, tumorset, TDIresults, postprobnoiselv, numoftargenes.cutoff = 5, DEGcallrate.thd = 0.2) {

  # Update the TDIresult
  idx2keep = intersect(intersect(which(is.element(TDIresults[, "patient_name"], tumorset)>0), which(TDIresults[, "cause_gene_name"]==driver)), which(TDIresults[, "posterior"]>=postprobnoiselv))
  TDIresults = TDIresults[idx2keep, ]
  numoftumors = length(tumorset)

  # Get DEG names
  DEGs = unique(TDIresults[, "result_gene_name"])
  numofDEGs = length(DEGs)

  # Create DriverDEGpair.countmatrix
  DriverDEGpair.countmatrix = matrix(0, 1, numofDEGs)
  rownames(DriverDEGpair.countmatrix) = driver
  colnames(DriverDEGpair.countmatrix) = DEGs

  # Get all the tumor TDI result files and update the counts in DriverDEGpair.countmatrix
  drivertumors = 0
  for (i in 1:numoftumors) {
    # TDI results for tumor i
    tumor.i = tumorset[i]
    # cat("Processing tumor", i, tumor.i, "...\n")
    idx.i = which(TDIresults[, "patient_name"]==tumor.i)
    TDIresults.i = TDIresults[idx.i, ]

    if (nrow(TDIresults.i)<numoftargenes.cutoff){
      next
    } else {
      DEGs.i = TDIresults.i[, "result_gene_name"]
      DriverDEGpair.countmatrix[driver, DEGs.i] = DriverDEGpair.countmatrix[driver, DEGs.i] + 1
      drivertumors = drivertumors + 1
    }
  }

  DriverDEGcallrate = DriverDEGpair.countmatrix/drivertumors
  idx.tardegs = which(DriverDEGcallrate >= DEGcallrate.thd)
  if (length(idx.tardegs)>0) {
    driver.tarDEGs = DriverDEGcallrate[driver, idx.tardegs]
    driver.tarDEGs = sort(driver.tarDEGs, decreasing=T)
  } else {
    driver.tarDEGs = c()
  }

  # Return the Driver.DEGlist
  return(driver.tarDEGs)
}



