#' @title Extract the significant TDI Tumor-SGA-DEG triplet
#'
#' @description Use only the significant TDI drivers and their target DEGs to extract the Tumor-SGA-DEG triplet.
#'
#'
#' @param TDIresults.triplet original TDI Tumor-SGA-DEG triplet
#' @param sigdrivers significant drivers
#' @param sigdrivers.tarDEGs target DEGs of the significant drivers
#' @param numofcores Number of CPU cores to be requested for parallel computing
#'
#' @export
#'
#' @return sigTDIresults.triplet
#'
#' @examples \dontrun{
#'
#' }
#'

extract.sigDriversTarDEGs.TSDtriplet.parallel <- function(TDIresults.triplet, sigdrivers, sigdrivers.tarDEGs, numofcores = 6) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # Names of the significant drivers
  sigdrivers.names = rownames(sigdrivers)
  numofsigdrivers = length(sigdrivers.names)

  # Extract the significant TSD triplet for each significant driver
  sigTDIresults.triplet <- foreach(driver.i = sigdrivers.names, .combine="rbind") %dopar% {
    tarDEGs.i = names(sigdrivers.tarDEGs[[driver.i]])
    idx.i = intersect(which(TDIresults.triplet[, "cause_gene_name"]==driver.i), which(is.element(TDIresults.triplet[, "result_gene_name"], tarDEGs.i)>0))

    # Update the sigTDIresults.triplet
    if (length(idx.i)>1) {
      sigTDIresults.triplet.i = TDIresults.triplet[idx.i, ]
    } else {
      sigTDIresults.triplet.i = c()
    }

    sigTDIresults.triplet.i  # Returned value for each loop
  }

  #stop cluster
  stopCluster(cl)

  return(sigTDIresults.triplet)
}



