#' @title Identify Target DEGs Associated with Each Significant Driver
#'
#' @description If a DEG is a target of a driver, the causal links between the driver and the DEG has to present in at least 20% (by default) of the tumors where the driver is present.
#'
#' @param driverlist identify the target DEGs for the given list of drivers
#' @param drivercallmatrix drivers per tumor
#' @param TDIresults TDI Tumor-SGA-DEG triplets
#' @param postprobnoiselv posterior probability cutoff per SGA
#' @param numoftargenes.cutoff minimal number of target DEGs for an SGA to be called a driver in the tumor
#' @param DEGcallrate.thd minimal call rate for a DEG to be associated with a driver
#'
#' @export
#'
#' @return Driver.DEGlist
#'
#' @examples \dontrun{
#'
#' }
#'

identify.driverAssociatedDEGlist <- function(driverlist, drivercallmatrix, TDIresults, postprobnoiselv, numoftargenes.cutoff = 5, DEGcallrate.thd = 0.2) {
  # Get the tumor names
  tumors = unique(TDIresults[, "patient_name"])
  numoftumors = length(tumors)

  # Get DEG names
  DEGs = unique(TDIresults[, "result_gene_name"])
  numofDEGs = length(DEGs)

  # Driver call matrix for the given drivers
  drivercallmatrix = drivercallmatrix[, driverlist]
  drivercallfreq = colSums(drivercallmatrix)
  numofdrivers = length(driverlist)
  
  # Create DriverDEGpair.countmatrix
  DriverDEGpair.countmatrix = matrix(0, numofdrivers, numofDEGs)
  rownames(DriverDEGpair.countmatrix) = driverlist
  colnames(DriverDEGpair.countmatrix) = DEGs

  # Get all the tumor TDI result files and update the counts in DriverDEGpair.countmatrix
  for (i in 1:numoftumors) {
    # TDI results for tumor i
    tumor.i = tumors[i]
    cat("Processing tumor", i, tumor.i, "...\n")
    idx.i = which(TDIresults[, "patient_name"]==tumor.i)
    TDIresults.i = TDIresults[idx.i, ]

    # SGAs in tumor i
    SGAs.i = unique(TDIresults.i[, "cause_gene_name"])
    SGAs.i = intersect(SGAs.i, driverlist)
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

      # Target DEGs for SGA j in tumor i
      idx.ij = which(TDIresults.ij[, "posterior"]>=postprobnoiselv.j)
      if (length(idx.ij)>=numoftargenes.cutoff) {
        DEGs.ij = TDIresults.ij[idx.ij, "result_gene_name"]
        # Update the matrix
        DriverDEGpair.countmatrix[SGA.j, DEGs.ij] = DriverDEGpair.countmatrix[SGA.j, DEGs.ij] + 1
      }
    }
  }

  # Find the target DEGs for each driver
  Driver.DEGlist = list()
  for (i in 1:numofdrivers) {
    # Name, number of driver call events, driver deg pair total counts for driver i
    driver.i = driverlist[i]
    drivercallfreq.i = drivercallfreq[driver.i]
    driverdegpair.count.i = DriverDEGpair.countmatrix[driver.i, ]

    # Caculate the DEG call rate and get the DEGs with call rate higher than the threshold
    degcallrate.i = driverdegpair.count.i/drivercallfreq.i
    idx.tardegs.i = which(degcallrate.i>=DEGcallrate.thd)
    tardegs.i = degcallrate.i[idx.tardegs.i]
    if (length(tardegs.i)>0) {
      Driver.DEGlist[[driver.i]] = sort(tardegs.i, decreasing=T)
    }
  }

  # Return the Driver.DEGlist
  return(Driver.DEGlist)
}



