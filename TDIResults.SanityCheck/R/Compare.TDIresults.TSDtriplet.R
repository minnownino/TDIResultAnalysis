#' @title Compare TDI Results between the New Version and the Standard Version
#'
#' @description Compare TDI triplet results with the previously most reliable version. Check the top SGA changes for DEGs.
#'
#' @param TDIresults.triplet.new Tumor-SGA-DEG triplets extracted from the new TDI results version 
#' @param TDIresults.triplet.standard Tumor-SGA-DEG triplets extracted from the standard TDI results version 
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

Compare.TDIresults.TSDtriplet <- function(TDIresults.triplet.new, TDIresults.triplet.standard, fname.output) {
  # Tumor_DEG pairs
  TumorDEGpairs.new = paste(TDIresults.triplet.new[, "patient_name"], TDIresults.triplet.new[, "result_gene_name"])
  TumorDEGpairs.standard = paste(TDIresults.triplet.standard[, "patient_name"], TDIresults.triplet.standard[, "result_gene_name"])
  TumorDEGpairs.overlaps = intersect(TumorDEGpairs.new, TumorDEGpairs.standard)

  # Triplets
  tripletmerged.new = paste(TDIresults.triplet.new[, "patient_name"], paste(TDIresults.triplet.new[, "cause_gene_name"], TDIresults.triplet.new[, "result_gene_name"]))
  tripletmerged.standard = paste(TDIresults.triplet.standard[, "patient_name"], paste(TDIresults.triplet.standard[, "cause_gene_name"], TDIresults.triplet.standard[, "result_gene_name"]))
  tripletmerged.overlaps = intersect(tripletmerged.new, tripletmerged.standard)
  
  # Write the comarison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("1.1. Top SGA changes:\n")
  cat("Number of TDI triplet records in the new version:", nrow(TDIresults.triplet.new), "\n")
  cat("Number of TDI triplet records in the standard version:", nrow(TDIresults.triplet.standard), "\n")
  cat("Number of overlaped Tumor-DEG pairs:", length(TumorDEGpairs.overlaps), "\n")
  cat("Number of TDI triplet records overlaps between both versions:", length(tripletmerged.overlaps), "\n")
  cat("Percentage of top SGA changes for DEGs from standard version to new version:", (length(TumorDEGpairs.overlaps) - length(tripletmerged.overlaps))/length(TumorDEGpairs.overlaps), "\n\n\n")
  sink() # Close the file
}
