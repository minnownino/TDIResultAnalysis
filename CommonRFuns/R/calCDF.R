#' @title Calculate the CDF from PDF
#'
#' @description
#'
#' @details
#'
#'
#' @param data.PDF
#'
#'
#' @export
#'
#'
#' @return
#'
#'
#' @examples \dontrun{
#'
#' }
#'

calCDF <- function(data.PDF) {
  numofpts = length(data.PDF)
  data.CDF = 0*data.PDF
  data.CDF[1] = data.PDF[1]
  for (i in 2:numofpts) {
    data.CDF[i] = data.CDF[i-1] + data.PDF[i]
  }

  return(data.CDF)
}
