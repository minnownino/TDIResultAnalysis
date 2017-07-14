#' @title Extract the Tumor-SGA-DEG triplets from the TDI results.
#'
#' @description Generic function for extracting the Tumor-SGA-DEG records along with the posterior probability from the TDI results. This function goes through all tumors files, i.e. TDI results, to find out the SGA with the highest posterior probability assigned for each DEG.
#'
#' @details The output file contains four columns: 1. Tumor name, 2. SGA name, 3. DEG name, and 4. Posterior probability. Each row/record represents the SGA with the highest posterior probability regulating the DEG in that tumor.
#'
#' @param fpath.TDIresults The path to the folder that contains all the TDI results.
#' @param fname.output The output file name of the TSD triplet.
#'
#' @export
#'
#' @return none
#'
#' @examples \dontrun{
#'
#' }
#'

generate.TumorSGADEG.triplet <- function(fpath.TDIresults, fname.output) {
  # Tumor files
  tumorfiles = list.files(fpath.TDIresults)
  patients = gsub(".csv", "", tumorfiles)
  numofpats = length(patients)

  # Go through each patient
  for (i in 1:numofpats) {
    tumorfile.i = tumorfiles[i] # tumor file i
    patient.i = patients[i] # patient i
    cat("Processing patient ", i, " ", patient.i, "...\n")

    # TDI result, i.e. posterior probability matrix
    postprob.i = read.data(paste(fpath.TDIresults, tumorfile.i, sep=""))
    SGAs.i = rownames(postprob.i)
    DEGs.i = colnames(postprob.i)
    numofDEGs.i = length(DEGs.i)

    # Find the highest posterior probability for each DEG
    indx.SGAs.toppostprob = apply(postprob.i, 2, which.max)
    SGAs.toppostprob = SGAs.i[indx.SGAs.toppostprob]

    # TDI result table for patient i
    TDI.results.table.i = matrix("null", numofDEGs.i, 4)
    TDI.results.table.i[, 1] = patient.i
    TDI.results.table.i[, 2] = SGAs.toppostprob
    TDI.results.table.i[, 3] = DEGs.i
    TDI.results.table.i[, 4] = postprob.i[cbind(as.matrix(SGAs.toppostprob), t(DEGs.i))]

    # Output the TDI results table for patient i
    if (i==1) {
      colnames(TDI.results.table.i) = c("patient_name", "cause_gene_name", "result_gene_name", "posterior")
      write.table(TDI.results.table.i, file=fname.output, quote=F, sep=",", row.names=F, col.names=T)
    } else {
      write.table(TDI.results.table.i, file=fname.output, append=T, quote=F, sep=",", row.names=F, col.names=F)
    }
  }
}



