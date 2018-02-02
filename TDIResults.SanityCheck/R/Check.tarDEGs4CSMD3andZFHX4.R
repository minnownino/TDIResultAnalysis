#' @title Check whether Some Experimental Tested Target DEGs are Found in the New TDI Result List
#'
#' @description Check Check whether CDKN3, GJB2, MAD2L1, MMP11 and WNT11 are DEG targets of CSMD3; check whether DES, NPR1, CDCA3, GSTM5 and FHL1 are DEG targets of ZFHX4.
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

Check.tarDEGs4CSMD3andZFHX4 <- function(sigdrivers.new, sigdrivers.tarDEGs.new, fname.output) {
  # CSMD3 and ZFHX4: driver call frequency and call rate
  CSMD3.drivercalls = sigdrivers.new["CSMD3", ]
  ZFHX4.drivercalls = sigdrivers.new["ZFHX4", ]

  # TDI predicted target DEGs for CSMD3 and ZFHX4
  CSMD3.tarDEGs = names(sigdrivers.tarDEGs.new[["CSMD3"]])
  ZFHX4.tarDEGs = names(sigdrivers.tarDEGs.new[["ZFHX4"]])

  # Experimental tested target DEGs for CSMD3 and ZFHX4
  CSMD3.exptest.tarDEGs = c("CDKN3", "GJB2", "MAD2L1", "MMP11", "WNT11")
  ZFHX4.exptest.tarDEGs = c("DES", "NPR1", "CDCA3", "GSTM5", "FHL1")

  # Write the comparison results to the file
  sink(fname.output, append=T) # Open the file for writing
  cat("5. CSMD3 and ZFHX4 target DEGs: \n")
  cat("Driver call frequency and driver call rate for CSMD3:", CSMD3.drivercalls["#DriverEvents"], CSMD3.drivercalls["Drivercallrate"], "\n")
  cat("Driver call frequency and driver call rate for ZFHX4:", ZFHX4.drivercalls["#DriverEvents"], ZFHX4.drivercalls["Drivercallrate"], "\n")
  cat("Experimental tested target DEGs found in CSMD3:", paste(intersect(CSMD3.exptest.tarDEGs, CSMD3.tarDEGs), collapse="; "), "\n")
  cat("Experimental tested target DEGs not found in CSMD3:", paste(setdiff(CSMD3.exptest.tarDEGs, CSMD3.tarDEGs), collapse="; "), "\n")
  cat("Experimental tested target DEGs found in ZFHX4:", paste(intersect(ZFHX4.exptest.tarDEGs, ZFHX4.tarDEGs), collapse="; "), "\n")
  cat("Experimental tested target DEGs not found in ZFHX4:", paste(setdiff(ZFHX4.exptest.tarDEGs, ZFHX4.tarDEGs), collapse="; "), "\n\n\n")
  sink() # Close the file
}
