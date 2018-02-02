#' @title Extract the TDI Predicted RB1 Target DEGs
#'
#' @description Extract the TDI predicted RB1 target DEGs in order to check the PASTAA program (Thomas-Chollier et al., 2011): trap.molgen.mpg.de/PASTAA.htm. Find out whether there is E2F/E2F1 motif binding sites enrichment in the promoters of the RB1 target DEGs that TDI predicted. 
#'
#' @param sigdrivers.new a driver call matrix which contains the driver call frequency and driver call rate in the new version
#' @param sigdrivers.tarDEGs.new a list that contains TDI predicted target DEGs for each driver in the new version
#' @param fname.output a output file to keep the comparison record
#' @param fname.RB1tarDEGs a output file containing the TDI predicted RB1 regulated DEG list 
#'
#' @export NULL
#'
#' @return NULL
#'
#' @examples \dontrun{
#'
#' }
#'

Extract.RB1tarDEGs <- function(sigdrivers.new, sigdrivers.tarDEGs.new, fname.output, fname.RB1tarDEGs) {
  # RB1: driver call frequency and call rate
  RB1.drivercalls = sigdrivers.new["RB1", ]
  
  # TDI predicted target DEGs for RB1 and output to file.
  RB1.tarDEGs = names(sigdrivers.tarDEGs.new[["RB1"]])
  write.csv(RB1.tarDEGs, file=fname.RB1tarDEGs, quote=F, header=F)

  # Write the comparison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("7. RB1 target DEGs: \n")
  cat("Driver call frequency and driver call rate for RB1:", RB1.drivercalls["#DriverEvents"], RB1.drivercalls["Drivercallrate"], "\n")
  cat("Number of TDI predicted target DEGs regulated by RB1:", length(RB1.tarDEGs), "\n")
  cat("The most significant motif enriched in the DEGs: \t (Check the PASTAA program: trap.molgen.mpg.de/PASTAA.htm) \n") # Check the PASTAA program
  cat("Number of DEGs that contains the motif: \n\n\n")
  sink() # Close the file  
}
