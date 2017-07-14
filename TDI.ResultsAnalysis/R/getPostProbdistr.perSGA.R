#' @title Get the posterior probability for each SGA
#'
#' @description Collect all the posterior probability that assigned to any DEG and get the distribution for each SGA.
#'
#' @details 
#'
#'
#' @param fpath path to the folder that contains the TDI results.
#' @param bins a vector giving the breaking points between histogram cells.
#' @param fname.output output .RData file name that contains three variables, i.e., bins, postprob.count.SGAs (posterior probability distribution for all SGAs), postprob.count.A0 (posterior probability distribution for leaking node A0).
#'
#' @export
#'
#' @return none
#'
#' @examples \dontrun{
#'
#' }
#'

getPostProbdistr.perSGA <- function(fpath, bins=seq(0, 1, by=0.0001), fname.output) {
  ### Processing TDI results matrix
  # Get the list of tumor names
  files = list.files(fpath)
  files = files[which(regexpr("TCGA", files)>0)]
  numoftumors = length(files)

  # All counts of instances in different nonzeros posterior probability intervals
  postprob.count.SGAs = c() # All posterior probability frequency counts except leaking node A0
  postprob.count.A0 = matrix(0, 1, length(bins)-1) # Posterior probability frequency of the leaking node A0 and cancer type labels for PANCAN

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

    # Index for the leaking node
    indx.A0 = union(which(is.element(rownames(PostProb.j), "A0")), which(regexpr("CancerType", rownames(PostProb.j))>0))
    if (length(indx.A0)==0) {
      stop("Leaking node A0 was not found. Please double check the data!\n")
    }

    # Posterior probabilities excluding leaking nodes
    PostProb.SGAs.j = PostProb.j[-indx.A0,]
    SGAs.j = rownames(PostProb.SGAs.j)
    numofSGAs.j = length(SGAs.j)
    # Posterior probabilities for leaking nodes
    PostProb.A0.j = PostProb.j[indx.A0, ]

    # Add new SGAs to the postprob.count.SGAs matrix
    if (length(postprob.count.SGAs)==0) {
      # Initialization
      postprob.count.SGAs = matrix(0, numofSGAs.j, length(bins)-1)
      rownames(postprob.count.SGAs) = SGAs.j
    } else {
      # Add new SGAs from tumor j
      SGAs2add.j = setdiff(SGAs.j, rownames(postprob.count.SGAs))
      numofSGAs2add.j = length(SGAs2add.j)
      if (numofSGAs2add.j>0) {
        postprob.count.SGAs = rbind(matrix(0, numofSGAs2add.j, length(bins)-1), postprob.count.SGAs)
        rownames(postprob.count.SGAs)[1:numofSGAs2add.j] = SGAs2add.j
      }
    }

    # Update the postprob.count.SGAs matrix
    for (n in 1:numofSGAs.j) {
      SGA.n = SGAs.j[n]
      hist.SGA.n = hist(PostProb.SGAs.j[SGA.n, ], bins, plot=FALSE)
      postprob.count.SGAs[SGA.n, ] = postprob.count.SGAs[SGA.n, ] + hist.SGA.n$counts
    }

    # Update the postprob.count.A0 matrix
    hist.A0.j = hist(PostProb.A0.j, bins, plot=FALSE)
    postprob.count.A0 = postprob.count.A0 + hist.A0.j$counts
  }

  ### Save the frequency counting to .RData format
  save(bins, postprob.count.SGAs, postprob.count.A0, file=fname.output)
}
