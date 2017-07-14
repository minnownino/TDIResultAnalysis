#' @title Adjust the Driver per Tumor Call by Only Keep Significant Target DEGs
#'
#' @description Use only the significant target DEGs as the target DEGs for each driver. Make the driver call using only the set of DEGs.
#'
#' @param TSDtriplet.filtered TSD triplet filtered by posterior probability threshold
#' @param drivercallpertumor Original driver call per tumor matrix
#' @param sigdriverstarDEGs the list of significant drivers and their significant target DEGs
#' @param numofsigDEGs threshold of number of target DEGs used to identify the SGAs as drivers
#'
#' @return sigdrivercallpertumor.adjusted
#'
#' @examples \dontrun{
#'
#' }
#'

adjusted.sigdrivercallpertumor <- function(TSDtriplet.filtered, drivercallpertumor, sigdriverstarDEGs, numofsigDEGs = 5) {
  # Significant drivers and their significant target DEGs
  sigdrivers = names(sigdriverstarDEGs)
  numofsigdrivers = length(sigdrivers)
  sigdrivercallpertumor = drivercallpertumor[, sigdrivers]
  sigdrivercallpertumor.adjusted = 0*sigdrivercallpertumor

  # Go through all drivers
  for (i in 1:numofsigdrivers) {
    # Driver i and its regulated DEGs
    driver.i = sigdrivers[i]
    tarDEGs.i = names(sigdriverstarDEGs[[driver.i]])

    # Triplets filtered by regulated DEGs
    idx.i = which(TSDtriplet.filtered[, "cause_gene_name"]==driver.i)
    TSDtriplet.i = TSDtriplet.filtered[idx.i, ]
    idx2keep.i = which(is.element(TSDtriplet.i[, "result_gene_name"], tarDEGs.i))
    if (length(idx2keep.i)>1) {
      TSDtriplet.i = TSDtriplet.i[idx2keep.i, ]
    } else {
      next # No triplets with significant target DEGs are found
    }

    # Go through all tumors that originally that the SGA is called as a driver
    tumors.i = names(which(sigdrivercallpertumor[, driver.i]==1))
    numoftumors.i = length(tumors.i)
    for (j in 1:numoftumors.i) {
      # Tumor j
      tumor.j = tumors.i[j]
      idx.ij = which(TSDtriplet.i[, "patient_name"]==tumor.j)
      if (length(idx.ij)==0) {
        # No significant DEGs found in the current tumor
        next
      } else {
        # One or more significant DEGs found in the current tumor
        if (length(idx.ij)>=numofsigDEGs) {
          sigdrivercallpertumor.adjusted[tumor.j, driver.i] = 1
        }
      }
    }
  }

  # Return the adjusted driver call per tumor matrix
  return(sigdrivercallpertumor.adjusted)
}
