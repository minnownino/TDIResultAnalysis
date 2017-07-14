#' @title Get the posterior probability for each SGA
#'
#' @description Collect all the posterior probability that assigned to any DEG and get the distribution for each SGA.
#'
#' @details
#'
#'
#' @param fpath path to the folder that contains the TDI results.
#' @param bins a vector giving the breaking points between histogram cells.
#' @param givenSGA name of the given SGA
#'
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

getPostProbdistr.givenSGA <- function(fpath, bins=seq(0, 1, by=0.0001), givenSGA) {
  ### Processing TDI results matrix
  # Get the list of tumor names
  files = list.files(fpath)
  files = files[which(regexpr("TCGA", files)>0)]
  numoftumors = length(files)

  # All counts of instances in different nonzeros posterior probability intervals for the given GT
  postprob.count.givenSGA = matrix(0, 1, length(bins)-1)

  # Go through all tumor samples
  for (j in 1:numoftumors) {
    # Tumor sample i
    tumor.j = files[j]
    cat("Preprocessing ", fpath, " sample.", j, ": ", tumor.j, "...\n", sep="")

    # Read the posterior probability for tumor j
    fname.PostProb.j = paste(fpath, "/", tumor.j, sep="")
    PostProb.j = read.data(fname.PostProb.j)

    # Stop if there is posterior probability larger than 1.01
    if (length(which(PostProb.j>1.01))>0) {
      stop("There are posterior probability larger than 1.01 for ", tumor.j, "!\n")
    }


    # Posterior probability distribution of the given GT in the current tumor
    if (!is.element(givenSGA, rownames(PostProb.j))) {
      next # Go to next one if the given GT is not found in the current tumor
    } else {
      hist.givenSGA.j = hist(PostProb.j[givenSGA, ], bins, plot=FALSE)
      postprob.count.givenSGA = postprob.count.givenSGA + hist.givenSGA.j$counts
    }
  }

  ### Return the count matrix
  return(postprob.count.givenSGA)
}
