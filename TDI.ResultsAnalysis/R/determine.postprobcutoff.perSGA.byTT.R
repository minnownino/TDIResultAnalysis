#' @title Determine the Posterior Probability Cutoff for Each SGA by Tissue Types
#'
#' @description Use the random permutation posterior probability distribution generated from function "getPostProbdistr.perSGA.byTT.parallel" to determine the posterior probability cutoff for each SGA by setting an appropriate significance level. The posterior probability is determined for each tissue type.
#'
#' @param fnames.postprobdistr.randperm.byTT the file names of the posterior probability distribution per SGA for different tissue types from random permutation experiments
#' @param siglv significant level of posterior probability cutoff from random permutation experiments
#' @param numofcores Number of CPU cores to be requested for parallel computing (6 by default)
#' @param postprobthd.bound the lower and upper bound for the posterior probability threshold (NULL by default)
#' @param fname.output .csv file name of the posterior probability cutoff for each tissue type (NULL by default)
#'
#' @export
#'
#' @return postprobcutoff.perSGA.output
#'
#' @examples \dontrun{
#'
#' }
#'

determine.postprobcutoff.perSGA.byTT <- function(fnames.postprobdistr.randperm.byTT, siglv = 0.001, numofcores = 6, postprobthd.bound = NULL, fname.output = NULL) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # Go through each simulation result
  numofsimus = length(fnames.postprobdistr.randperm.byTT)
  postprob.count.SGAs.randperms = list()
  for (i in 1:numofsimus) {
    cat("Loading", fnames.postprobdistr.randperm.byTT[i], "...\n")
    load(fnames.randperms[i])
    postprob.count.SGAs.randperms[[i]] = postprob.count.SGAs
  }

  # Tissue types
  tts = names(postprob.count.SGAs.randperms[[1]])
  numoftts = length(tts)
  postprobthd.byTT.list <- foreach(tt.i = tts) %dopar% {
    # Go through all random permutated data
    for (j in 1:numofsimus) {
      # Random permutated data j for tissue type i
      postprob.count.SGAs.ij = postprob.count.SGAs.randperms[[j]][[tt.i]]

      # Calculate the CDF
      postprob.pdf.SGAs.ij = postprob.count.SGAs.ij/rowSums(postprob.count.SGAs.ij) # PDF
      postprob.cdf.SGAs.ij = t(apply(postprob.pdf.SGAs.ij, 1, cumsum)) # CDF

      # Find the posterior probability cutoff for each SGA
      indx.postprobcutoff.ij = apply(abs(1-siglv-postprob.cdf.SGAs.ij), 1, which.min)
      postprobcutoff.perSGA.ij = as.matrix(bins[indx.postprobcutoff.ij])

      # Update the postprob.cutoff matrix
      if (j==1) {
        postprobcutoff.perSGA.i = postprobcutoff.perSGA.ij
      } else {
        postprobcutoff.perSGA.i = cbind(postprobcutoff.perSGA.i, postprobcutoff.perSGA.ij)
      }
    }

    # Posterior probability threshold for tissue type i
    if (numofsimus==1) {
      postprobthd.byTT.i = postprobcutoff.perSGA.i
    } else {
      postprobthd.byTT.i = rowMeans(postprobcutoff.perSGA.i)
    }
    postprobthd.byTT.i = as.matrix(postprobthd.byTT.i)
    rownames(postprobthd.byTT.i) = rownames(postprob.count.SGAs.ij)

    # Returned value for each loop
    postprobthd.byTT.i
  }

  #stop cluster
  stopCluster(cl)

  # Extract the all the SGAs that appeared in at least one tissue type
  SGAnames.byTT = lapply(postprobthd.byTT.list, rownames)
  SGAnames = unique(unlist(SGAnames.byTT))
  numofSGAs = length(SGAnames)

  # Create and update the posterior probability threshold by tissue type matrix
  postprobthd.byTT = matrix(NA, numofSGAs, numoftts)
  rownames(postprobthd.byTT) = SGAnames
  colnames(postprobthd.byTT) = tts
  for (i in 1:numoftts) {
    tt.i = tts[i]
    postprobthd.byTT.i = postprobthd.byTT.list[[i]]
    sganames.i = rownames(postprobthd.byTT.i)
    postprobthd.byTT[sganames.i, tt.i] = postprobthd.byTT.i
  }

  # Check if there is required upper and lower bounds for the posterior probability threshold
  if (!is.null(postprobthd.bound)) {
    # Apply the posterior probability bound
    postprobthd.byTT[which(postprobthd.byTT<postprobthd.bound[1])] = postprobthd.bound[1]
    postprobthd.byTT[which(postprobthd.byTT>postprobthd.bound[2])] = postprobthd.bound[2]
  }

  # Output the postprobthd.byTT to .csv file
  if (!is.null(fname.output)) {
    write.csv(postprobthd.byTT, file =fname.output, quote=F)
  }

  # Return the posterior probability threshold by tissue types
  return(postprobthd.byTT)
}

