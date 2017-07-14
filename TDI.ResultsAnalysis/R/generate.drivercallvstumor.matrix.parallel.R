#' @title Generate Driver Calls per Tumor (parallel computing)
#'
#' @description Identify an SGA as a driver for a given tumor, if the number of its target DEGs in the tumor is larger or equal to the given cutoff.
#'
#' @param TDIresults TDI Tumor-SGA-DEG triplets
#' @param postprobnoiselv posterior probability cutoff for each SGA
#' @param numoftargenes.cutoff minimal number of target DEGs for an SGA to be called a driver in the tumor
#' @param numofcores Number of CPU cores to be requested for parallel computing
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

generate.drivercallvstumor.matrix.parallel <- function(TDIresults, postprobnoiselv, numoftargenes.cutoff=5, numofcores=6) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # Get the tumor names
  tumors = unique(TDIresults[, "patient_name"])
  numoftumors = length(tumors)

  # Get the SGA names
  SGAs = unique(TDIresults[, "cause_gene_name"])
  numofSGAs = length(SGAs)

  # Go through each tumor
  Driverspertumor.list <- foreach(tumor.i=tumors, .inorder=FALSE) %dopar% {
    # TDI results for tumor i
    idx.i = which(TDIresults[, "patient_name"]==tumor.i)
    TDIresults.i = TDIresults[idx.i, ]

    # SGAs in tumor i
    SGAs.i = unique(TDIresults.i[, "cause_gene_name"])
    SGAs.i = setdiff(SGAs.i, "A0") # Exlcude A0
    SGAs.i = intersect(SGAs.i, rownames(postprobnoiselv))
    numofSGAs.i = length(SGAs.i)

    # Go through all SGAs in tumor i
    drivers.i = c()
    if (numofSGAs.i==0) {
      # no potential drivers found
    } else {
      # potential drivers found
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
          drivers.i = c(drivers.i, SGA.j)
        }
      }
    }

    # drivers4tumor.i is the returned value for each loop
    drivers4tumor.i = list()
    drivers4tumor.i[[tumor.i]] = drivers.i
    drivers4tumor.i
  }

  #stop cluster
  stopCluster(cl)

  # Create drivercall versus tumor matrix
  drivercall.vs.tumor.matrix = matrix(0, numoftumors, numofSGAs)
  rownames(drivercall.vs.tumor.matrix) = tumors
  colnames(drivercall.vs.tumor.matrix) = SGAs
  for (i in 1:numoftumors) {
    drivers4tumor.i = Driverspertumor.list[[i]]

    # No potential drivers
    if (length(drivers4tumor.i)==0) {
      next
    }

    # At least one potential driver(s)
    tumor.i = names(drivers4tumor.i)
    drivers.i = drivers4tumor.i[[tumor.i]]
    if (length(drivers.i)>0) {
      drivercall.vs.tumor.matrix[tumor.i, drivers.i] = 1
    }
  }

  # Return the matrix
  return(drivercall.vs.tumor.matrix)
}

