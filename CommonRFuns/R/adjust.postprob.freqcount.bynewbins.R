#' @title Adjust the frequency count by using new bins for posterior probability frequency distribution
#'
#' @description
#'
#' @details
#'
#'
#' @param bins.original
#' @param postprob.freqcount.old
#' @param bins.new
#'
#'
#' @export
#'
#'
#' @return A vector that contains the the new counts with respect to the new bins.
#'
#'
#' @examples \dontrun{
#'
#' }
#'

adjust.postprob.freqcount.bynewbins <- function(bins.old, postprob.freqcount.old, bins.new) {
  # Create a new posterior probability frequency count matrix
  GTs = rownames(postprob.freqcount.old)
  numofGTs = nrow(postprob.freqcount.old)
  numofbins.new = length(bins.new)
  postprob.freqcount.new = matrix(NA, numofGTs, numofbins.new)
  rownames(postprob.freqcount.new) = GTs

  for (i in 1:numofbins.new) {
    # Find the indices of the original bins to merge
    if (i==1) {
      indx2merge = which(bins.old<=bins.new[i])
    } else {
      bin.lowbound = bins.new[i-1]
      bin.upperbound = bins.new[i]
      indx2merge = intersect(which(bins.old>bin.lowbound), which(bins.old<=bin.upperbound))
    }

    if (length(indx2merge)==1) {
      postprob.freqcount.new[, i] = postprob.freqcount.old[, indx2merge]
    } else if (length(indx2merge)>1) {
      postprob.freqcount.new[, i] = rowSums(postprob.freqcount.old[, indx2merge], na.rm=T)
    } else {
      stop("New bins exceed the range of the old bins!")
    }
  }

  # Return new frequency counting matrix
  return(postprob.freqcount.new)
}
