#' @title Extract the Posterior Probabilities for the Given SGA-DEG Pair
#'
#' @description Extract the posterior probabilities for the given SGA-DEG pair in a tumor set
#'
#' @param sgadegpair the given SGA-DEG pair
#' @param tumorset the tumor set that may contain the SGA-DEG pair
#' @param fpath.TDI the file path to the TDI results
#' @param numofcores Number of CPU cores to be requested for parallel computing (6 by default)
#'
#' @export
#'
#' @return postprobs
#'
#' @examples \dontrun{
#'
#' }
#'

extract.postprobs.givenSGADEGpair.parallel <- function(sgadegpair, tumorset, fpath.TDI, numofcores = 6) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # Go through all tumors in the tumor set
  postprobs <- foreach(tumor.i = tumorset, .combine=c) %dopar% {
    fname.i = paste(fpath.TDI, tumor.i, ".csv", sep="")
    tdiresults.i = read.data(fname.i)
    postprob.i = tdiresults.i[sgadegpair[1], sgadegpair[2]]
    postprob.i # Returned value for each loop
  }
  names(postprobs) = tumorset

  #stop cluster
  stopCluster(cl)

  # Return value for the function
  return(postprobs)
}
