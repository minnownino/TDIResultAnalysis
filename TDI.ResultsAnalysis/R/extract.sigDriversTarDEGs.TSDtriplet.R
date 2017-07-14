#' @title Extract the significant TDI Tumor-SGA-DEG triplet
#'
#' @description Use only the significant TDI drivers and their target DEGs to extract the Tumor-SGA-DEG triplet.
#'
#'
#' @param TDIresults.triplet original TDI Tumor-SGA-DEG triplet
#' @param sigdrivers significant drivers
#' @param sigdrivers.tarDEGs target DEGs of the significant drivers
#' @param postprobnoiselv posterior probability cutoff per SGA
#'
#' @export
#'
#' @return sigTDIresults.triplet
#'
#' @examples \dontrun{
#'
#' }
#'

extract.sigDriversTarDEGs.TSDtriplet <- function(TDIresults.triplet, sigdrivers, sigdrivers.tarDEGs, postprobnoiselv) {
  # Names of the significant drivers
  sigdrivers.names = rownames(sigdrivers)
  numofsigdrivers = length(sigdrivers.names)

  # Extract the significant TSD triplet for each significant driver
  sigTDIresults.triplet = c()
  for (i in 1:numofsigdrivers) {
    driver.i = sigdrivers.names[i]
    postprobnoiselv.i = postprobnoiselv[driver.i, "PostProbcutoff"]
    tarDEGs.i = names(sigdrivers.tarDEGs[[driver.i]])
    idx.i = intersect(intersect(which(TDIresults.triplet[, "cause_gene_name"]==driver.i), which(is.element(TDIresults.triplet[, "result_gene_name"], tarDEGs.i)>0)), which(TDIresults.triplet[, "posterior"]>=postprobnoiselv.i))

    # Update the sigTDIresults.triplet
    if (length(idx.i)>1) {
      sigTDIresults.triplet = rbind(sigTDIresults.triplet, TDIresults.triplet[idx.i, ])
    }
  }

  return(sigTDIresults.triplet)
}



