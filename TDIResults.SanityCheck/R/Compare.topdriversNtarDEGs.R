#' @title Compare TDI Identified Most Significant Drivers and Their Target DEGs between the New Version and the Standard Version
#'
#' @description Compare TDI triplet results with the previously most reliable version. Compare the overlaps of the most significant drivers between the new version and standard version. For the overlaped top drivers, compare their target DEGs overlap.
#'
#' @param sigdrivers.new a driver call matrix which contains the driver call frequency and driver call rate in the new version
#' @param sigdrivers.tarDEGs.new a list that contains TDI predicted target DEGs for each driver in the new version
#' @param sigdrivers.standard a driver call matrix which contains the driver call frequency and driver call rate in the standard version
#' @param sigdrivers.tarDEGs.standard a list that contains TDI predicted target DEGs for each driver in the standard version
#' @param topdriverfreqcutoff driver frequency cutoff for selecting top drivers
#' @param topdrivercallratecutoff driver call rate cutoff for selecting top drivers
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

Compare.topdriversNtarDEGs <- function(sigdrivers.new, sigdrivers.tarDEGs.new, sigdrivers.standard, sigdrivers.tarDEGs.standard, topdriverfreqcutoff, topdrivercallratecutoff, fname.output) {
  # Find the top drivers for the standard version
  idx.topdrivers.new = intersect(which(sigdrivers.new[, "#DriverEvents"]>=topdriverfreqcutoff), which(sigdrivers.new[, "Drivercallrate"]>=topdrivercallratecutoff))
  topdrivers.new = rownames(sigdrivers.new)[idx.topdrivers.new]

  # Find the top drivers for the standard version
  idx.topdrivers.standard = intersect(which(sigdrivers.standard[, "#DriverEvents"]>=topdriverfreqcutoff), which(sigdrivers.standard[, "Drivercallrate"]>=topdrivercallratecutoff))
  topdrivers.standard = rownames(sigdrivers.standard)[idx.topdrivers.standard]

  # Intersections and mismatch
  topdrivers.intersect = intersect(topdrivers.new, topdrivers.standard)
  numoftopdrivers.intersect = length(topdrivers.intersect)
  topdrivers.onlyinstandard = setdiff(topdrivers.standard, topdrivers.new)
  topdrivers.onlyinnew = setdiff(topdrivers.new, topdrivers.standard)

  # Write the comparison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat(paste("2.1. Top drivers(driver call frequencty>=", topdriverfreqcutoff, "; driver call rate>=", topdrivercallratecutoff, "):", sep=""), "\n")
  cat("Number of top drivers in the new version:", length(topdrivers.new), "\n")
  cat("Number of top drivers in the standard version:", length(topdrivers.standard), "\n")
  cat("Number of top dirvers overlap between the new version and the standard version:", length(topdrivers.intersect), "\n")
  if (length(topdrivers.onlyinstandard)>0) {
  	cat(paste(topdrivers.onlyinstandard, collapse=";"), "are not identified as the top drivers in the new version!", "\n")
  }
  if (length(topdrivers.onlyinnew)>0) {
  	cat(paste(topdrivers.onlyinnew, collapse=";"), "are identified as the top drivers in the new version only!", "\n")
  }
  cat("\n")
  sink() # Close the file

  # Write the overlaped top driver target DEGs comparison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("2.2. Target DEGs of the Overlaped Top drivers: \n")
  cat("Driver name\t# target DEGs in new version\t# target DEGs in standard version\t# overlaps\n")
  for (i in 1:numoftopdrivers.intersect) {
  	driver.i = topdrivers.intersect[i]
  	drivertarDEGs.new.i = names(sigdrivers.tarDEGs.new[[driver.i]])
  	drivertarDEGs.standard.i = names(sigdrivers.tarDEGs.standard[[driver.i]])
  	numoftarDEGs.new.i = length(drivertarDEGs.new.i)
  	numoftarDEGs.standard.i = length(drivertarDEGs.standard.i)
  	numoftarDEGs.overlap.i = length(intersect(drivertarDEGs.new.i, drivertarDEGs.standard.i))

  	cat(driver.i, "\t", numoftarDEGs.new.i, "\t", numoftarDEGs.standard.i, "\t", numoftarDEGs.overlap.i, "\n")
  }
  cat("\n\n")
  sink() # Close the file
}

