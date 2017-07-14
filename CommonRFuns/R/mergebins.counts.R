#' @title
#'
#' @description
#'
#' @details
#'
#'
#' @param bins.original
#' @param counts.original
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

mergebins.counts <- function(bins.original, counts.original, bins.new) {
  # Check whether the new bins has the same range as the old bins
  if (bins.new[1]!=bins.original||bins.new[length(bins.new)]!=bins.original[length(bins.original)]) {
  	stop("The new bins and the original bins do not have the same range!")
  }

  # 



  return(counts.new)
}
