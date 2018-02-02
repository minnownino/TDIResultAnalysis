#' @title Check whether Drivers Involved in Notch Pathway Share Common Target DEGs
#'
#' @description Check whether drivers involved in the Notch pathway, i.e. Notch1, EP300, CREBBP and NCOR2, share common target DEGs.
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

Check.tarDEGsoverlapsamongpathwaydrivers.notch <- function(sigdrivers.new, sigdrivers.tarDEGs.new, fname.output) {
  # Write the comparison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("6.2. Notch pathway involved drivers and their target DEGs:\n")

  # notch pathway genes
  drivers.notchpathway = c("NOTCH1", "EP300", "CREBBP", "NCOR2")
  # Genes that are not identified as significant drivers in notch pathway
  drivers.nosig = setdiff(drivers.notchpathway, rownames(sigdrivers.new))
  if (length(drivers.nosig)>0) {
    cat("Notch pathway genes that are not identified as significant drivers:", paste(drivers.nosig, collapse=";"), "\n")
  }
  drivers.notchpathway = intersect(rownames(sigdrivers.new), drivers.notchpathway)
  numofdrivers.notchpathway = length(drivers.notchpathway)

  # Driver call frequency and call rate
  drivercalls = sigdrivers.new[drivers.notchpathway, ]

  # TDI predicted target DEGs for drivers involved in Notch pathway
  tarDEGs.notchpathway = sigdrivers.tarDEGs.new[drivers.notchpathway]

  # Table: Number of DEG overlaps between different drivers
  numofDEGoverlaps = matrix(NA, numofdrivers.notchpathway, numofdrivers.notchpathway)
  rownames(numofDEGoverlaps) = drivers.notchpathway
  colnames(numofDEGoverlaps) = drivers.notchpathway
  for (i in 1:(numofdrivers.notchpathway-1)) {
    driver.i = drivers.notchpathway[i]
    tarDEGs.i = names(tarDEGs.notchpathway[[driver.i]])
    for (j in (i+1):numofdrivers.notchpathway) {
      driver.j = drivers.notchpathway[j]
      tarDEGs.j = names(tarDEGs.notchpathway[[driver.j]])
      numofDEGoverlaps[driver.i, driver.j] = length(intersect(tarDEGs.i, tarDEGs.j))
    }
  }

  # Driver call frequency and call rate
  cat("\n")
  for (i in 1:numofdrivers.notchpathway) {
    driver.i = drivers.notchpathway[i]
    cat(paste("Driver call frequency and driver call rate for ", driver.i, ": ", drivercalls[driver.i, "#DriverEvents"], " ", drivercalls[driver.i, "Drivercallrate"], sep=""), "\n")
  }

  # Number of target DEGs for each driver in Notch pathway
  cat("\n")
  for (i in 1:numofdrivers.notchpathway) {
    driver.i = drivers.notchpathway[i]
    cat("Number of TDI predicted target DEGs regulated by", driver.i, ":", length(tarDEGs.notchpathway[[driver.i]]), "\n")
  }

  # Number of target DEGs overlaps between drivers in Notch pathway
  cat("\nNumber of target DEGs overlaps between drivers in Notch pathway:\n")
  cat("\t", paste(drivers.notchpathway, collapse="\t"), "\n")
  for (i in 1:numofdrivers.notchpathway) {
  	driver.i = drivers.notchpathway[i]
  	cat(driver.i, "\t", paste(numofDEGoverlaps[driver.i, ], collapse="\t"), "\n")
  }
  cat("\n\n")
  sink() # Close the file
}



