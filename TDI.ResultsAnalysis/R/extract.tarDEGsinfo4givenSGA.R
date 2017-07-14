#' @title Extract information of target DEGs for the given SGA/driver.
#'
#' @description Given an SGA, find all the DEGs that are predicted as its target by TDI. Count the number of tumors TDI predict there is a causal link between the given SGA and DEG; collect the posterior probability distribution for each DEG; calculate the mean and variance of the posterior probability.
#'
#' @param givenSGA the name of the given SGA
#' @param drivingtumorset tumors that the given SGA is called a driver
#' @param TDIresults TDI Tumor-SGA-DEG triplets
#' @param postprobnoiselv.givenSGA posterior probability cutoff for the given SGA
#'
#' @export GivenSGA.tarDEGsinfo
#'
#' @return GivenSGA.tarDEGsinfo
#'
#' @examples \dontrun{
#'
#' }
#'

extract.tarDEGsinfo4givenSGA <- function(givenSGA, drivingtumorset, TDIresults, postprobnoiselv.givenSGA) {
  # Filter the TDIresults table to contain only the given SGA and its driving tumors
  idx2keep = intersect(which(TDIresults[, "cause_gene_name"]==givenSGA), which(is.element(TDIresults[, "patient_name"], drivingtumorset)>0))
  if (length(idx2keep)==0) {
    stop("There is no target DEGs found in TDI results for", givenSGA, "and its driving tumors!")
  } else if (length(idx2keep)==1) {
    print(TDIresults[idx2keep, ])
    stop("There is one target DEG information found in TDI results for", givenSGA, "and its driving tumors!")
  } else {
    TDIresults = TDIresults[idx2keep, ]
  }

  # All potential target DEGs
  tarDEGs = unique(TDIresults[, "result_gene_name"])
  numoftarDEGs = length(tarDEGs)

  # Go through and extract the information of all potential target DEGs for the given SGA
  tarDEGs.info = matrix(0, numoftarDEGs, 3)
  rownames(tarDEGs.info) = tarDEGs
  colnames(tarDEGs.info) = c("#tumors", "postprob.mean", "postprob.sd")
  tarDEGs.postprobdistr = list()
  for (i in 1:numoftarDEGs) {
    # TDI results for deg.i
    deg.i = tarDEGs[i]
    idx.i = which(TDIresults[, "result_gene_name"]==deg.i)

    # Count the number of tumors; extract the posterior probability distribution; calculate mean and variance.
    if (length(idx.i)==1) {
      # Only one tumor
      tarDEGs.info[deg.i, "#tumors"] = 1
      tarDEGs.postprobdistr[[deg.i]] = as.numeric(TDIresults[idx.i, "posterior"])
      tarDEGs.info[deg.i, "postprob.mean"] = tarDEGs.postprobdistr[[deg.i]]
      tarDEGs.info[deg.i, "postprob.sd"] = 0
    } else {
      # two or more tumors
      TDIresults.i = TDIresults[idx.i, ]
      tarDEGs.info[deg.i, "#tumors"] = nrow(TDIresults.i)
      tarDEGs.postprobdistr[[deg.i]] = as.numeric(TDIresults.i[, "posterior"])
      tarDEGs.info[deg.i, "postprob.mean"] = mean(tarDEGs.postprobdistr[[deg.i]])
      tarDEGs.info[deg.i, "postprob.sd"] = sd(tarDEGs.postprobdistr[[deg.i]])
    }
  }

  # Return value for the function
  GivenSGA.tarDEGsinfo = list()
  GivenSGA.tarDEGsinfo[["generalinfo"]] = tarDEGs.info
  GivenSGA.tarDEGsinfo[["postprobdistr"]] = tarDEGs.postprobdistr
  return(GivenSGA.tarDEGsinfo)
}
