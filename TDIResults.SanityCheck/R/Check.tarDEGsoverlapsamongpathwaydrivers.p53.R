#' @title Check whether Drivers Involved in P53 Pathway Share Common Target DEGs
#'
#' @description Check whether drivers involved in the p53 pathway, i.e. TP53, MDM4, PTEN, CCNE1, CDKN2A, PPM1D, CASP8, CCND1 and CDK4, share common target DEGs.
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

Check.tarDEGsoverlapsamongpathwaydrivers.p53 <- function(sigdrivers.new, sigdrivers.tarDEGs.new, fname.output) {
  sink(fname.output, append=T) # Open the file for writing
  cat("6.1. p53 pathway involved drivers and their target DEGs:\n")

  # p53 pathway genes
  drivers.p53pathway = c("TP53", "MDM4", "PTEN", "CCNE1", "CDKN2A", "PPM1D", "CASP8", "CCND1", "CDK4")
  # Genes that are not identified as significant drivers in p53 pathway
  drivers.nosig = setdiff(drivers.p53pathway, rownames(sigdrivers.new))
  if (length(drivers.nosig)>0) {
    cat("p53 pathway genes that are not identified as significant drivers:", paste(drivers.nosig, collapse=";"), "\n")
  }
  drivers.p53pathway = intersect(rownames(sigdrivers.new), drivers.p53pathway)
  numofdrivers.p53pathway = length(drivers.p53pathway)

  # Driver call frequency and call rate
  drivercalls = sigdrivers.new[drivers.p53pathway, ]

  # TDI predicted target DEGs for drivers involved in p53 pathway
  tarDEGs.p53pathway = sigdrivers.tarDEGs.new[drivers.p53pathway]

  # Table: Number of DEG overlaps between different drivers
  numofDEGoverlaps = matrix(NA, numofdrivers.p53pathway, numofdrivers.p53pathway)
  rownames(numofDEGoverlaps) = drivers.p53pathway
  colnames(numofDEGoverlaps) = drivers.p53pathway
  for (i in 1:(numofdrivers.p53pathway-1)) {
    driver.i = drivers.p53pathway[i]
    tarDEGs.i = names(tarDEGs.p53pathway[[driver.i]])
    for (j in (i+1):numofdrivers.p53pathway) {
      driver.j = drivers.p53pathway[j]
      tarDEGs.j = names(tarDEGs.p53pathway[[driver.j]])
      numofDEGoverlaps[driver.i, driver.j] = length(intersect(tarDEGs.i, tarDEGs.j))
    }
  }

  # Write the comparison results to the file
  # Driver call frequency and call rate
  cat("\n")
  for (i in 1:numofdrivers.p53pathway) {
    driver.i = drivers.p53pathway[i]
    cat(paste("Driver call frequency and driver call rate for ", driver.i, ": ", drivercalls[driver.i, "#DriverEvents"], " ", drivercalls[driver.i, "Drivercallrate"], sep=""), "\n")
  }

  # Number of target DEGs for each driver in p53 pathway
  cat("\n")
  for (i in 1:numofdrivers.p53pathway) {
    driver.i = drivers.p53pathway[i]
    cat("Number of TDI predicted target DEGs regulated by", driver.i, ":", length(tarDEGs.p53pathway[[driver.i]]), "\n")
  }

  # Number of target DEGs overlaps between drivers in p53 pathway
  cat("\nNumber of target DEGs overlaps between drivers in p53 pathway:\n")
  cat("\t", paste(drivers.p53pathway, collapse="\t"), "\n")
  for (i in 1:numofdrivers.p53pathway) {
  	driver.i = drivers.p53pathway[i]
  	cat(driver.i, "\t", paste(numofDEGoverlaps[driver.i, ], collapse="\t"), "\n")
  }
  cat("\n\n")
  sink() # Close the file
}



