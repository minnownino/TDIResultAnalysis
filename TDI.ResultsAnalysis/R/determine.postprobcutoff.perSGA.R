#' @title Determine the Posterior Probability Cutoff for Each SGA
#'
#' @description Use the random permutation posterior probability distribution generated from function "getPostProbdistr.perSGA" to determine the posterior probability cutoff for each SGA by setting an appropriate significance level.
#'
#' @param fnames.postprobdistr.randperm the file names of the posterior probability distribution per SGA from random permutation experiments
#' @param siglv significant level of posterior probability cutoff from random permutation experiments
#' @param postprobthd.bound the lower and upper bound for the posterior probability threshold (NULL by default)
#' @param fname.output .csv file name of the posterior probability cutoff (NULL by default)
#'
#' @export
#'
#' @return postprobcutoff.perSGA.output
#'
#' @examples \dontrun{
#'
#' }
#'

determine.postprobcutoff.perSGA <- function(fnames.postprobdistr.randperm, siglv = 0.001, postprobthd.bound = NULL, fname.output = NULL) {
  ### Calculate the postperior probability cutoff per SGA
  # Go through each simulation result
  numofrandperm.simus = length(fnames.postprobdistr.randperm)
  for (j in 1:numofrandperm.simus) {
    # Load simulation result j
    cat("Loading", fnames.postprobdistr.randperm[j], "...\n")
    load(fnames.postprobdistr.randperm[j])
    if (j==1) {
      SGAnames = rownames(postprob.count.SGAs)
    }
    postprob.pdf.SGAs = postprob.count.SGAs/rowSums(postprob.count.SGAs) # PDF
    postprob.cdf.SGAs = t(apply(postprob.pdf.SGAs, 1, cumsum)) # CDF

    # Find the posterior probability cutoff for each SGA
    indx.postprobcutoff = apply(abs(1-siglv-postprob.cdf.SGAs), 1, which.min)
    postprobcutoff.perSGA.j = as.matrix(bins[indx.postprobcutoff])
    rownames(postprobcutoff.perSGA.j) = SGAnames

    # Update the postprob.cutoff matrix
    if (j==1) {
      postprobcutoff.perSGA = postprobcutoff.perSGA.j
    } else {
      postprobcutoff.perSGA = cbind(postprobcutoff.perSGA, postprobcutoff.perSGA.j)
    }
    rm(postprob.pdf.SGAs, postprob.cdf.SGAs)
    gc()
  }
  postprobcutoff.perSGA = cbind(postprobcutoff.perSGA, matrix(0,nrow(postprobcutoff.perSGA),2))
  simu.no = seq(1, numofrandperm.simus, by = 1)
  colnames(postprobcutoff.perSGA) = c(paste("PostProbcutoff.Simu.", simu.no, sep=""), "PostProbcutoff.mean", "PostProbcutoff.sd")
  postprobcutoff.perSGA[, "PostProbcutoff.mean"] = rowMeans(postprobcutoff.perSGA[, paste("PostProbcutoff.Simu.", simu.no, sep="")])
  postprobcutoff.perSGA[, "PostProbcutoff.sd"] = apply(postprobcutoff.perSGA[, paste("PostProbcutoff.Simu.", simu.no, sep="")], 1, sd)

  # Posterior Probability cutoff per SGA
  postprobcutoff.perSGA.output = as.matrix(postprobcutoff.perSGA[, c("PostProbcutoff.mean")])
  colnames(postprobcutoff.perSGA.output) = "PostProbcutoff"
  if (!is.null(postprobthd.bound)) {
    # Apply the posterior probability bound
    postprobcutoff.perSGA.output[which(postprobcutoff.perSGA.output<postprobthd.bound[1]), "PostProbcutoff"] = postprobthd.bound[1]
    postprobcutoff.perSGA.output[which(postprobcutoff.perSGA.output>postprobthd.bound[2]), "PostProbcutoff"] = postprobthd.bound[2]
  }
  if (!is.null(fname.output)) {
    write.csv(postprobcutoff.perSGA.output, file =fname.output, quote=F)
  }

  # Return the posterior probability threshold per SGA
  return(postprobcutoff.perSGA)
}

