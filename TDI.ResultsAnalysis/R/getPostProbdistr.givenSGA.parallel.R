#' @title Get the posterior probability for the given SGA (parallel computing)
#'
#' @description Collect all the posterior probability that assigned to any DEG and get the distribution for the given SGA.
#'
#' @details
#'
#' @param fpath path to the folder that contains the TDI results.
#' @param bins a vector giving the breaking points between histogram cells.
#' @param givenSGA name of the given SGA
#' @param numofcores Number of CPU cores to be requested for parallel computing
#' @param deletefiles a Boolean value to declare whether or not to delete files that do not have the given SGA. Default is FALSE.
#'
#' @export
#'
#'
#' @return postprob.count.givenSGA
#'
#'
#' @examples \dontrun{
#'
#' }
#'

getPostProbdistr.givenSGA.parallel <- function(fpath, bins=seq(0, 1, by=0.001), givenSGA, numofcores=6, deletefiles=FALSE) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  ### Processing TDI results matrix
  # Get the list of tumor names
  files = list.files(fpath)
  files = files[which(regexpr("TCGA", files)>0)]

  # Go through all patient TDI data in parallel
  postprob.count.givenSGA <- foreach(file.i=files, .combine="+") %dopar% {
    # Read the posterior probability for tumor i
    fname.PostProb.i = paste(fpath, "/", file.i, sep="")
    PostProb.i = read.data(fname.PostProb.i)
    sganames.i = rownames(PostProb.i)

    # Check whether the given SGA is in the current patient TDI data
    if (is.element(givenSGA, sganames.i)) { # Found in the current patient TDI data
      # Stop if there is posterior probability larger than 1.01
      if (length(which(PostProb.i>1.01))>0) {
        stop("There are posterior probability larger than 1.01 for ", file.i, "!\n")
      }

      # Histogram counts
      hist.i = hist(PostProb.i[givenSGA, ], bins, plot=FALSE)
      postprob.count.givenSGA.i = hist.i$counts

    } else { # Not found in the current patient TDI data
      # Skip the current patient by assigning all zeros to the counts
      postprob.count.givenSGA.i = matrix(0, 1, length(bins)-1)

      # Delete the current patient file if the deletefile is true
      if (deletefiles) {
        file.remove(fname.PostProb.i)
      }
    }

    postprob.count.givenSGA.i # Return value for each loop
  }

  #stop cluster
  stopCluster(cl)

  ### Return the count matrix
  return(postprob.count.givenSGA)
}


