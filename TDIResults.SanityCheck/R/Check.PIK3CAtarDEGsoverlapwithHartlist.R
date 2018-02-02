#' @title Check whether TDI Predicted PIK3CA Target DEGs Overlaps with the Hart list
#'
#' @description Extract and compare the TDI predicted PIK3CA regulated DEGs with the Hart list. Find out how significant two list overlap with each other.
#'
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

Check.PIK3CAtarDEGsoverlapwithHartlist <- function(sigdrivers.new, sigdrivers.tarDEGs.new, fname.output) {
  # PIK3CA: driver call frequency and call rate
  PIK3CA.drivercalls = sigdrivers.new["PIK3CA", ]
  
  # TDI predicted PIK3CA target DEGs
  PIK3CA.tarDEGs = names(sigdrivers.tarDEGs.new[["PIK3CA"]])

  # PIK3CA affected DEGs in Hart list
  numofDEGs.Hartlist = length(Hartlist)

  # DEG overlaps between TDI predictions and Hart list
  DEGoverlaps = intersect(PIK3CA.tarDEGs, Hartlist)

  # Write the comparison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("4. TDI predicted PIK3CA target DEGs and Hart predicted PIK3CA DEGs: \n")
  cat("Driver call frequency and driver call rate for PIK3CA:", PIK3CA.drivercalls["#DriverEvents"], PIK3CA.drivercalls["Drivercallrate"], "\n")
  cat("Number of TDI predicted target DEGs regulated by PIK3CA:", length(PIK3CA.tarDEGs), "\n")
  cat("Number of PIK3CA affected DEGs in Hart list:", numofDEGs.Hartlist, "\n")
  cat("Number of DEG overlaps between TDI predictions and Hart list:", length(DEGoverlaps), "\n\n\n")
  sink() # Close the file
}

