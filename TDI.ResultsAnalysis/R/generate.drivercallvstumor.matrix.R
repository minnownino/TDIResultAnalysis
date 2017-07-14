#' @title Generate Driver Calls per Tumor
#'
#' @description Identify an SGA as a driver for a given tumor, if the number of its target DEGs in the tumor is larger or equal to the given cutoff.
#'
#' @param TDIresults TDI Tumor-SGA-DEG triplets
#' @param postprobnoiselv posterior probability cutoff for each SGA
#' @param numoftargenes.cutoff minimal number of target DEGs for an SGA to be called a driver in the tumor
#'
#'
#' @export
#'
#'
#' @return drivercall.vs.tumor.matrix
#'
#'
#' @examples \dontrun{
#'
#' }
#'

generate.drivercallvstumor.matrix <- function(TDIresults, postprobnoiselv, numoftargenes.cutoff=5) {
  # Get the tumor names
  tumors = unique(TDIresults[, "patient_name"])
  numoftumors = length(tumors)

  # Get the SGA names
  SGAs = unique(TDIresults[, "cause_gene_name"])
  numofSGAs = length(SGAs)

  # Create drivercall versus tumor matrix
  drivercall.vs.tumor.matrix = matrix(0, numoftumors, numofSGAs)
  rownames(drivercall.vs.tumor.matrix) = tumors
  colnames(drivercall.vs.tumor.matrix) = SGAs

  # Go through each tumor
  for (i in 1:numoftumors) {
    # TDI results for tumor i
    tumor.i = tumors[i]
    cat("Processing tumor", i, tumor.i, "...\n")
    idx.i = which(TDIresults[, "patient_name"]==tumor.i)
    TDIresults.i = TDIresults[idx.i, ]

    # SGAs in tumor i
    SGAs.i = unique(TDIresults.i[, "cause_gene_name"])
    SGAs.i = setdiff(SGAs.i, "A0") # Exlcude A0
    SGAs.i = intersect(SGAs.i, rownames(postprobnoiselv))
    numofSGAs.i = length(SGAs.i)

    if (numofSGAs.i==0) {
      next # no driver found
    }

    # Go through all SGAs in tumor i
    for (j in 1:numofSGAs.i) {
      SGA.j = SGAs.i[j]
      idx.j = which(TDIresults.i[, "cause_gene_name"]==SGA.j)
      if (length(idx.j)<numoftargenes.cutoff) {
        next
      }

      TDIresults.ij = TDIresults.i[idx.j, ]
      postprobnoiselv.j = postprobnoiselv[SGA.j, "PostProbcutoff"]
      idx.ij = which(TDIresults.ij[, "posterior"]>=postprobnoiselv.j)
      if (length(idx.ij)>=numoftargenes.cutoff) {
        # Update the matrix
        drivercall.vs.tumor.matrix[tumor.i, SGA.j] = 1
      }
    }
  }

  # Save the matrix
  return(drivercall.vs.tumor.matrix)
}

