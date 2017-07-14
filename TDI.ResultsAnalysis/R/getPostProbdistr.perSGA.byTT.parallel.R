#' @title Get the posterior probability for each SGA by tissue type (Parallel computing)
#'
#' @description Collect all the posterior probability that assigned to any DEG and get the distribution for each SGA by tissue type.
#'
#' @details
#'
#'
#' @param fpath path to the folder that contains the TDI results.
#' @param bins a vector giving the breaking points between histogram cells.
#' @param sganames the gene name of SGAs
#' @param tumorsgroupedbyTTs tumors are grouped into different tissue types
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

getPostProbdistr.perSGA.byTT.parallel <- function(fpath, bins=seq(0, 1, by=0.001), sganames, tumorsgroupedbyTTs, numofcores = 4, fname.output) {
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

  # Go through all tissue types
  tissuetypes = names(tumorsgroupedbyTTs)
  numoftts = length(tissuetypes)
  postprob.count.SGAs = list()
  postprob.count.A0 = list()
  for (n in 1:numoftts) {
    # Tissue type n
    tt.n = tissuetypes[n]
    tumors.n = tumorsgroupedbyTTs[[tt.n]]
    files.n = intersect(paste(tumors.n, ".csv", sep=""), files) # Tumors that present the current tissue type
    if (length(files.n)==0) {
      warning("There is no files for tissue type:", tt.n)
      next
    }

    # Go through all patient TDI data in parallel for the current tissue type
    postprob.count.list.n <- foreach(file.i=files.n, .inorder=FALSE) %dopar% {
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

    # Combine counts from all tumors into a single matrix
    postprob.count.n = matrix(0, numofSGAs, length(bins)-1)
    rownames(postprob.count.n) = sganames
    for (i in 1:length(postprob.count.list.n)) {
      postprob.count.i = postprob.count.list.n[[i]]
      sganames.i = rownames(postprob.count.i)
      postprob.count.n[sganames.i, ] = postprob.count.n[sganames.i, ] + postprob.count.i
    }
    rm(postprob.count.list.n)
    gc() # Garbage collection

    # Exclude the SGAs with no counts
    idx2del = which(rowSums(postprob.count.n)==0)
    if (length(idx2del)>0) {
      postprob.count.n = postprob.count.n[-idx2del, ]
    }
    idx.A0 = which(rownames(postprob.count.n)=="A0") # Index for leaking node A0
    if (length(idx.A0)==0) {
      stop("There is no leaking node 'A0' in posterior probability count matrix!")
    }
    postprob.count.SGAs.n = postprob.count.n[-idx.A0, ]
    postprob.count.A0.n = postprob.count.n[idx.A0, ]

    postprob.count.SGAs[[tt.n]] = postprob.count.SGAs.n
    postprob.count.A0[[tt.n]] = postprob.count.A0.n
  }

  #stop cluster
  stopCluster(cl)

  # Output results
  save(bins, postprob.count.SGAs, postprob.count.A0, file=fname.output)
}
