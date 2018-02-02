#' @title Compare TDI Results for Significant Drivers and Their Target DEGs between the New Version and the Standard Version
#'
#' @description Compare TDI triplet results with the previously most reliable version. Check the top SGA changes for DEGs.
#'
#' @param TDIresults.sigtriplet.new Tumor-SGA-DEG triplets for significant drivers and their target DEGs extracted from the new TDI results version 
#' @param TDIresults.sigtriplet.standard Tumor-SGA-DEG triplets for significant drivers and their target DEGs extracted from the standard TDI results version 
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

Compare.TDIresults.sigTSDtriplet <- function(TDIresults.sigtriplet.new, TDIresults.sigtriplet.standard, fname.output) {
  # Tumor_DEG pairs
  TumorDEGpairs.new = paste(TDIresults.sigtriplet.new[, "patient_name"], TDIresults.sigtriplet.new[, "result_gene_name"])
  TumorDEGpairs.standard = paste(TDIresults.sigtriplet.standard[, "patient_name"], TDIresults.sigtriplet.standard[, "result_gene_name"])
  TumorDEGpairs.overlaps = intersect(TumorDEGpairs.new, TumorDEGpairs.standard)

  # Triplets
  tripletmerged.new = paste(TDIresults.sigtriplet.new[, "patient_name"], paste(TDIresults.sigtriplet.new[, "cause_gene_name"], TDIresults.sigtriplet.new[, "result_gene_name"]))
  tripletmerged.standard = paste(TDIresults.sigtriplet.standard[, "patient_name"], paste(TDIresults.sigtriplet.standard[, "cause_gene_name"], TDIresults.sigtriplet.standard[, "result_gene_name"]))
  tripletmerged.overlaps = intersect(tripletmerged.new, tripletmerged.standard)
  
  # Write the comarison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("1.2. Top driver changes for their target DEGs:\n")
  cat("Number of TDI triplet records for significant drivers and their target DEGs in the new version:", nrow(TDIresults.sigtriplet.new), "\n")
  cat("Number of TDI triplet records for significant drivers and their target DEGs in the standard version:", nrow(TDIresults.sigtriplet.standard), "\n")
  cat("Number of overlaped Tumor-DEG pairs:", length(TumorDEGpairs.overlaps), "\n")
  cat("Number of TDI triplet records overlaps between both versions:", length(tripletmerged.overlaps), "\n")
  cat("Percentage of changes for significant drivers and their target DEGs from standard version to new version:", (length(TumorDEGpairs.overlaps) - length(tripletmerged.overlaps))/length(TumorDEGpairs.overlaps), "\n\n\n")
  sink() # Close the file
}
