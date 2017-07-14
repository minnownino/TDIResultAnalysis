#' @title
#'
#' @description
#'
#' @details
#'
#'
#' @param M
#' @param g1
#' @param g2
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

merge.matrix <- function(M, g1, g2) {
  gnames = colnames(M)
  idx.g1 = which(gnames==g1)
  idx.g2 = which(gnames==g2)

  idx.tumors = which((M[, idx.g1] + M[, idx.g2])>0)
  M[, idx.g1] = 0
  M[idx.tumors, idx.g1] = 1
  M = M[, -idx.g2]

  return(M)
}
