#' @title Get the posterior probability for each SGA (Parallel computing)
#'
#' @description Collect all the posterior probability that assigned to any DEG and get the distribution for each SGA.
#'
#' @details
#'
#'
#' @param fpath path to the folder that contains the TDI results.
#' @param bins a vector giving the breaking points between histogram cells.
#' @param sganames the gene name of SGAs
#' @param numofcores Number of CPU cores to be requested for parallel computing
#' @param fname.output .RData file name for the output
#'
#' @export
#'
#' @return none
#'
#' @examples \dontrun{
#'
#' }
#'

getPostProbdistr.perSGA.parallel <- function(fpath, bins=seq(0, 1, by=0.001), sganames, numofcores = 4, fname.output) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  ### Processing TDI results matrix
  # Get the list of tumor names
  files = list.files(fpath)
  files = files[which(regexpr("TCGA", files)>0)]
  numoftumors = length(files)

  # Add A0 to the sga list
  sganames = c(sganames, "A0")
  numofSGAs = length(sganames)

  # Go through all patient TDI data in parallel
  postprob.count.list <- foreach(file.i=files, .inorder=FALSE) %dopar% {
    # Read the posterior probability for tumor i
    fname.PostProb.i = paste(fpath, "/", file.i, sep="")
    PostProb.i = read.data(fname.PostProb.i)
    sganames.i = rownames(PostProb.i)
    numofSGAs.i = length(sganames.i)

    # Stop if there is posterior probability larger than 1.01
    if (length(which(PostProb.i>1.01))>0) {
      stop("There are posterior probability larger than 1.01 for ", file.i, "!\n")
    }

    # Index for the leaking node
    indx.A0 = union(which(is.element(rownames(PostProb.i), "A0")), which(regexpr("CancerType", rownames(PostProb.i))>0))
    if (length(indx.A0)==0) {
      stop("Leaking node A0 was not found. Please double check the data!\n")
    }

    # Posterior probability distribution for tumor i
    postprob.count.i = matrix(0, numofSGAs.i, length(bins)-1)
    rownames(postprob.count.i) = sganames.i

    # Update the postprob.count.SGAs matrix
    for (n in 1:numofSGAs.i) {
      SGA.n = sganames.i[n]
      hist.SGA.n = hist(PostProb.i[SGA.n, ], bins, plot=FALSE)
      postprob.count.i[SGA.n, ] = hist.SGA.n$counts
    }

    postprob.count.i # Return value for each loop
  }

  #stop cluster
  stopCluster(cl)

  # Combine counts from all tumors into a single matrix
  postprob.count = matrix(0, numofSGAs, length(bins)-1)
  rownames(postprob.count) = sganames
  for (i in 1:length(postprob.count.list)) {
    postprob.count.i = postprob.count.list[[i]]
    sganames.i = rownames(postprob.count.i)
    postprob.count[sganames.i, ] = postprob.count[sganames.i, ] + postprob.count.i
  }
  rm(postprob.count.list)
  gc() # Garbage collection

  # Exclude the SGAs with no counts
  idx2del = which(rowSums(postprob.count)==0)
  if (length(idx2del)>0) {
    postprob.count = postprob.count[-idx2del, ]
  }
  idx.A0 = which(rownames(postprob.count)=="A0") # Index for leaking node A0
  if (length(idx.A0)==0) {
    stop("There is no leaking node 'A0' in posterior probability count matrix!")
  }
  postprob.count.SGAs = postprob.count[-idx.A0, ]
  postprob.count.A0 = postprob.count[idx.A0, ]

  # Output results
  save(bins, postprob.count.SGAs, postprob.count.A0, file=fname.output)
}
