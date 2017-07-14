#' @title Extract the Posterior Probabilities for the Given SGA to the Tumor-DEG pairs.
#'
#' @description
#'
#' @param SGA the name of the given SGA
#' @param tumordegpairs the name of tumor and deg pair to search for the posterior probability
#' @param fpath.TDI the name of the folder that contains the TDI outputs
#' @param numofcores Number of CPU cores to be requested for parallel computing
#'
#' @export postprobs
#'
#' @return
#'
#' @examples \dontrun{
#'
#' }
#'

extract.tarDEGspostprobfromTDIresults.givenSGA.parallel <- function(SGA, tumordegpairs, fpath.TDI, numofcores = 6) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # Get the tumor names
  tumors = unique(tumordegpairs[, "patient_name"]) # Function "unique" returns a vector which will be in the same order as the order of the first unique elements
  numoftumors = length(tumors)
  tumorsinthefolder = gsub(".csv", "", dir(fpath.TDI))
  tumorsnotinthefolder = setdiff(tumors, tumorsinthefolder)
  if (length(tumorsnotinthefolder)>0) {
    stop("There are tumors not in the TDI results folder!")
  }

  postprobs <- foreach(tumor.i = tumors, .combine = "c") %dopar% {
    # Read the TDI result for tumor.i from the file
    fname.i = paste(fpath.TDI, tumor.i, ".csv", sep="")
    TDIresultmatrix.i = read.data(fname.i)

    # Get the name of the degs
    idx.tumor.i = which(tumordegpairs[, "patient_name"]==tumor.i)
    degs.i = tumordegpairs[idx.tumor.i, "result_gene_name"]

    # Get the posterior probabilies
    postprob.i = TDIresultmatrix.i[SGA, degs.i]
    postprob.i # Returned value for each loop
  }


  #stop cluster
  stopCluster(cl)

  return(postprobs)
}

