#' @title Extract the TSD Triplets that Satisfy the Posterior Probability Threshold
#'
#' @description Use the posterior probability threshold per SGA to filter the TSD triplets.
#'
#' @param TSDtriplet the original TSDtriplet
#' @param postprobthd.perSGA posterior probability threshold per SGA
#' @param numofcores Number of CPU cores to be requested for parallel computing
#'
#' @return TSDtriplet.filtered
#'
#' @examples \dontrun{
#'
#' }
#'

extract.TSDtriplet.bypostprobthd <- function(TSDtriplet, postprobthd.perSGA, numofcores = 4) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # SGAs that are found as the top regulators
  topSGAs = unique(TSDtriplet[, "cause_gene_name"])
  topSGAs = setdiff(topSGAs, "A0") # Exclude null node
  numoftopSGAs = length(topSGAs)

  # Keep the triplet that satisfy the posterior probability threshold
  idx2keep <- foreach(sga.i=topSGAs, .combine=c) %dopar% {
    # Find the triplet records to keep for sga.i
    postprobthd.i = postprobthd.perSGA[sga.i, "PostProbcutoff"]
    idx.i = which(TSDtriplet[, "cause_gene_name"]==sga.i)
    if (length(idx.i)==1) {
      if (TSDtriplet[idx.i, "posterior"]>=postprobthd.i) {
        idx2keep.i = idx.i
      } else {
        idx2keep.i = c()
      }
    } else if (length(idx.i)>1) {
      TSDtriplet.i = TSDtriplet[idx.i, ]
      idx2keep.i = idx.i[which(TSDtriplet.i[, "posterior"]>=postprobthd.i)]
    } else {
      stop(sga.i, "not found in the TSD triplet!")
    }

    idx2keep.i # Return value for each loop
  }

  #stop cluster
  stopCluster(cl)

  # Filter the TSD triplet
  TSDtriplet.filtered = TSDtriplet[idx2keep, ]

  return(TSDtriplet.filtered)
}

