#' @title Check the Target DEGs Overlaps between drivers KEAP1 and NFE2L2 for the New Version
#'
#' @description Count the number tumors that KEAP1 and NFE2L2 are called drivers and calculate their driver call rate. 
#'
#' @param drivercallpertumor.new a driver call per tumor matrix where each row represents a tumor and each column represents a driver
#' @param sigdrivers.new a driver call matrix which contains the driver call frequency and driver call rate in the new version
#' @param sigdrivers.tarDEGs.new a list that contains TDI predicted target DEGs for each driver in the new version
#' @param fname.output a output file to keep the comparison record
#'
#' @export NULL
#'
#' @return NULL
#'
#' @examples \dontrun{
#'
#' }
#'

Check.tarDEGsoverlapsbetweenKEAP1andNFE2L2 <- function(drivercallpertumor.new, sigdrivers.new, sigdrivers.tarDEGs.new, fname.output) {
  # KEAP1: driver call frequency and call rate
  KEAP1.drivercalls = sigdrivers.new["KEAP1", ]

  # NFE2L2: driver call frequency and call rate
  NFE2L2.drivercalls = sigdrivers.new["NFE2L2", ]
  
  # Number of tumors that are driven by both KEAP1 and NFE2L2
  KEAP1.drivingtumors = names(which(drivercallpertumor.new[, "KEAP1"]==1))
  NFE2L2.drivingtumors = names(which(drivercallpertumor.new[, "NFE2L2"]==1))
  KEAP1NFE2L2.drivingtumors = intersect(KEAP1.drivingtumors, NFE2L2.drivingtumors)

  # TDI predicted target DEGs for KEAP1 and NFE2L2, and their overlaps
  KEAP1.tarDEGs = names(sigdrivers.tarDEGs.new[["KEAP1"]])
  NFE2L2.tarDEGs = names(sigdrivers.tarDEGs.new[["NFE2L2"]])
  KEAP1NFE2L2.tarDEGs = intersect(KEAP1.tarDEGs, NFE2L2.tarDEGs)

  # Write the comparison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("3. Driver KEAP1 and NFE2L2: \n")
  cat("Driver call frequency and driver call rate for KEAP1:", KEAP1.drivercalls["#DriverEvents"], KEAP1.drivercalls["Drivercallrate"], "\n")
  cat("Driver call frequency and driver call rate for NFE2L2:", NFE2L2.drivercalls["#DriverEvents"], NFE2L2.drivercalls["Drivercallrate"], "\n")
  cat("Number of tumors driven by both KEAP1 and NFE2L2:", length(KEAP1NFE2L2.drivingtumors), "\n")
  cat("Number of target DEGs regulated by KEAP1:", length(KEAP1.tarDEGs), "\n")
  cat("Number of target DEGs regulated by NFE2L2:", length(NFE2L2.tarDEGs), "\n")
  cat("Number of target DEGs regulated by KEAP1 both NFE2L2:", length(KEAP1NFE2L2.tarDEGs), "\n\n\n")
  sink() # Close the file
}
