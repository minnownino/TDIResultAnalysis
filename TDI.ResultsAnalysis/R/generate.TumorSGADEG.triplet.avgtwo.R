#' @title
#'
#' @description
#'
#' @details
#'
#'
#'
#'
#'
#' @param fpath.TDIresults.1
#' @param fpath.TDIresults.2
#' @param fname.output
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

generate.TumorSGADEG.triplet.avgtwo <- function(fpath.TDIresults.1, fpath.TDIresults.2, fname.output) {
  # Tumor files for file path 1
  tumorfiles.1 = list.files(fpath.TDIresults.1)
  patients.1 = gsub(".csv", "", tumorfiles.1)

  # Tumor files for file path 2
  tumorfiles.2 = list.files(fpath.TDIresults.2)
  patients.2 = gsub(".csv", "", tumorfiles.2)

  # Check whether two folders contain exactly the same set of patients
  if (length(union(patients.1, patients.2))>length(intersect(patients.1, patients.2))) {
    stop("There are mismatch between two patient sets!\n")
  }
  numofpats = length(patients.1)

  # Go through each patient
  for (i in 1:numofpats) {
    # patient i
    tumorfile.i = tumorfiles.1[i]
    patient.i = patients.1[i] 
    cat("Processing patient ", i, " ", patient.i, "...\n")

    # TDI result from folder 1, i.e. posterior probability matrix
    postprob1.i = read.data(paste(fpath.TDIresults.1, tumorfile.i, sep=""))
    SGAs1.i = rownames(postprob1.i)
    DEGs1.i = colnames(postprob1.i)
    
    # TDI result from folder 2, i.e. posterior probability matrix
    postprob2.i = read.data(paste(fpath.TDIresults.2, tumorfile.i, sep=""))
    SGAs2.i = rownames(postprob2.i)
    DEGs2.i = colnames(postprob2.i)
    
    # Check if two matrices contain the exact same set of SGAs and DEGs
    if (length(union(SGAs1.i, SGAs2.i))>length(intersect(SGAs1.i, SGAs2.i))) {
      stop("There are mismatch between two SGAs sets in ", patient.i, "!")
    }
    if (length(union(DEGs1.i, DEGs2.i))>length(intersect(DEGs1.i, DEGs2.i))) {
      stop("There are mismatch between two DEGs sets in ", patient.i, "!")
    }

    # Calculate the average between two datasets
    SGAs.i = intersect(SGAs1.i, SGAs2.i)
    DEGs.i = intersect(DEGs1.i, DEGs2.i)
    numofDEGs.i = length(DEGs.i)
    postprob1.i = postprob1.i[SGAs.i, ]
    postprob1.i = postprob1.i[, DEGs.i]
    postprob2.i = postprob2.i[SGAs.i, ]
    postprob2.i = postprob2.i[, DEGs.i]
    postprob.i = (postprob1.i + postprob2.i)/2

    # Find the highest posterior probability for each DEG
    indx.SGAs.toppostprob = apply(postprob.i, 2, which.max)
    SGAs.toppostprob = SGAs.i[indx.SGAs.toppostprob]

    # TDI result table for patient i
    TDI.results.table.i = matrix("null", numofDEGs.i, 4)
    TDI.results.table.i[, 1] = patient.i
    TDI.results.table.i[, 2] = SGAs.toppostprob
    TDI.results.table.i[, 3] = DEGs.i
    TDI.results.table.i[, 4] = postprob.i[cbind(as.matrix(SGAs.toppostprob), DEGs.i)]

    # Output the TDI results table for patient i
    if (i==1) {
      colnames(TDI.results.table.i) = c("patient_name", "cause_gene_name", "result_gene_name", "posterior")
      write.table(TDI.results.table.i, file=fname.output, quote=F, sep=",", row.names=F, col.names=T)
    } else {
      write.table(TDI.results.table.i, file=fname.output, append=T, quote=F, sep=",", row.names=F, col.names=F)
    }
  }
}



