#' @title Generate Driver Calls per Tumor (parallel computing) Using the Filtered TSD triplets
#'
#' @description Identify an SGA as a driver for a given tumor, if the number of its target DEGs in the tumor is larger or equal to the given cutoff.
#'
#' @param TDIresults.filtered TDI Tumor-SGA-DEG triplets filtered by posterior probability threshold
#' @param numoftargenes.cutoff minimal number of target DEGs for an SGA to be called a driver in the tumor
#' @param numofcores Number of CPU cores to be requested for parallel computing
#'
#' @export
#'
#' @return drivercall.vs.tumor.matrix
#'
#'
#' @examples \dontrun{
#'
#' }
#'

generate.drivercallvstumor.matrix.parallel.v2 <- function(TDIresults.filtered, numoftargenes.cutoff=5, numofcores=6) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # Get the tumor names
  tumors = unique(TDIresults.filtered[, "patient_name"])
  numoftumors = length(tumors)

  # Get the SGA names
  SGAs = unique(TDIresults.filtered[, "cause_gene_name"])
  numofSGAs = length(SGAs)

  if (numoftargenes.cutoff==1) { # only 1 target DEG
    tumors_SGAs = unique(TSDtriplet.filtered[, c("patient_name", "cause_gene_name")])
    tumors = unique(TSDtriplet.filtered[, "patient_name"])
    SGAs = unique(TSDtriplet.filtered[, "cause_gene_name"])
    drivercall.vs.tumor.matrix = matrix(0, length(tumors), length(SGAs))
    rownames(drivercall.vs.tumor.matrix) = tumors
    colnames(drivercall.vs.tumor.matrix) = SGAs
    drivercall.vs.tumor.matrix[tumors_SGAs] = 1
  } else { # More than just 1 target DEG
    # Go through each tumor
    Driverspertumor.list <- foreach(tumor.i=tumors, .inorder=FALSE) %dopar% {
      # TDI results for tumor i
      drivers.i = c()
      idx.i = which(TDIresults.filtered[, "patient_name"]==tumor.i)
      if (length(idx.i)==1) {
        # No potential drivers identified. Do nothing.
      } else {
        TDIresults.i = TDIresults.filtered[idx.i, ]

        # SGAs in tumor i
        SGAs.i = unique(TDIresults.i[, "cause_gene_name"])
        numofSGAs.i = length(SGAs.i)

        # Go through all SGAs in tumor i
      
        if (numofSGAs.i==0) {
          # no potential drivers found
        } else {
          # potential drivers found
          for (j in 1:numofSGAs.i) {
            SGA.j = SGAs.i[j]
            idx.j = which(TDIresults.i[, "cause_gene_name"]==SGA.j)
            if (length(idx.j)>=numoftargenes.cutoff) {
              drivers.i = c(drivers.i, SGA.j)
            }
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
  }

  # Return the matrix
  return(drivercall.vs.tumor.matrix)
}

