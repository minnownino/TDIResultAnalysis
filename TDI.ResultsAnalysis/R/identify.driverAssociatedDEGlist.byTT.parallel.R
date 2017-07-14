#' @title Identify Target DEGs Associated with Each Significant Driver
#'
#' @description If a DEG is a target of a driver, the causal links between the driver and the DEG has to present in at least 20% (by default) of the tumors where the driver is present.
#'
#' @param driverlist identify the target DEGs for the given list of drivers
#' @param drivercallmatrixbyTT drivers per tumor in tissue type specific manner
#' @param TDIresults TDI Tumor-SGA-DEG triplets
#' @param postprobnoiselvbyTT posterior probability cutoff per SGA by tissue types
#' @param numoftargenes.cutoff minimal number of target DEGs for an SGA to be called a driver in the tumor
#' @param DEGcallrate.thd minimal call rate for a DEG to be associated with a driver
#' @param numofcores Number of CPU cores to be requested for parallel computing (6 by default)
#'
#' @export
#'
#' @return Driver.DEGlist
#'
#' @examples \dontrun{
#'
#' }
#'

identify.driverAssociatedDEGlist.byTT.parallel <- function(driverlist, drivercallmatrixbyTT, TDIresults, postprobnoiselvbyTT, numoftargenes.cutoff = 5, DEGcallrate.thd = 0.2, numofcores = 6) {
  # Set parallel computing enviroment
  cl <- makeCluster(numofcores)
  registerDoParallel(cl)

  # Driver call matrix for the given drivers
  drivercallmatrixbyTT = drivercallmatrixbyTT[, driverlist]
  drivercallfreq = colSums(drivercallmatrixbyTT)
  numofdrivers = length(driverlist)

  # Posterior probability threshold for the given drivers
  postprobnoiselvbyTT = postprobnoiselvbyTT[driverlist, ]

  # Tissue types
  tts = names(tumorsgroupedbyTT)
  numoftts = length(tts)

  # Get the tumor names and their corresponding tissue types
  tumors = unique(TDIresults[, "patient_name"])
  numoftumors = length(tumors)
  tissuetypes = matrix(NA, numoftumors, 1)
  rownames(tissuetypes) = tumors
  colnames(tissuetypes) = "tissuetype"
  for (n in 1:numoftts) {
    tt.n = tts[n]
    tumors.n = intersect(tumorsgroupedbyTT[[tt.n]], tumors)
    if (length(tumors.n)>0) {
      tissuetypes[tumors.n, "tissuetype"] = tt.n
    }
  }
  # Check if we find the tissue types for all tumors
  if (sum(is.na(tissuetypes))>0) {
    tumorsmissingtissuetypes = names(which(is.na(tissuetypes)))
    stop("Tissue types are missing for ", tumorsmissingtissuetypes, "...")
  }

  # Get DEG names
  DEGs = unique(TDIresults[, "result_gene_name"])
  numofDEGs = length(DEGs)

  # Get all the tumor TDI result files and extract a list of Driver-DEG pair counts for each tumor
  DriverDEGpair.countmatrix.list <- foreach(tumor.i=tumors, .inorder=F) %dopar% {
    # Get the tissue specific posterior probability threshold
    tt.i = tissuetypes[tumor.i, "tissuetype"]
    postprobnoiselv = as.matrix(postprobnoiselvbyTT[, tt.i])
    colnames(postprobnoiselv) = "PostProbcutoff"

    # Initialization of the returned matrix for tumor i
    DriverDEGpair.countmatrix.i = matrix(0, numofdrivers, numofDEGs)
    rownames(DriverDEGpair.countmatrix.i) = driverlist
    colnames(DriverDEGpair.countmatrix.i) = DEGs

    # TDI results for tumor i
    idx.i = which(TDIresults[, "patient_name"]==tumor.i)
    TDIresults.i = TDIresults[idx.i, ]

    # SGAs in tumor i
    SGAs.i = unique(TDIresults.i[, "cause_gene_name"])
    SGAs.i = intersect(SGAs.i, driverlist)
    numofSGAs.i = length(SGAs.i)

    if (numofSGAs.i==0) {
      # no driver found
    } else {
      # Go through all drivers in tumor i
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
          DriverDEGpair.countmatrix.i[SGA.j, DEGs.ij] = 1
        }
      }
    }

    # Update DriverDEGpair.countmatrix.i matrix
    if (sum(DriverDEGpair.countmatrix.i)==0) {
      DriverDEGpair.countmatrix.i = c()
    } else {
      idx2keep.degs = which(colSums(DriverDEGpair.countmatrix.i)>0)
      if (length(idx2keep.degs)==1) {
        DriverDEGpair.countmatrix.i = c()
      } else {
        DriverDEGpair.countmatrix.i = DriverDEGpair.countmatrix.i[, idx2keep.degs]
        idx2keep.sgas = which(rowSums(DriverDEGpair.countmatrix.i)>0)
        if (length(idx2keep.sgas)==1) {
          DriverDEGpair.countmatrix.i = t(as.matrix(DriverDEGpair.countmatrix.i[idx2keep.sgas, ]))
          rownames(DriverDEGpair.countmatrix.i) = driverlist[idx2keep.sgas]
        } else {
          DriverDEGpair.countmatrix.i = DriverDEGpair.countmatrix.i[idx2keep.sgas, ]
        }
      }
    }
    DriverDEGpair.countmatrix.i # Returned value for each loop
  }

  # Driver DEG pair count matrix for all tumors
  DriverDEGpair.countmatrix = matrix(0, numofdrivers, numofDEGs)
  rownames(DriverDEGpair.countmatrix) = driverlist
  colnames(DriverDEGpair.countmatrix) = DEGs
  for (i in 1:length(DriverDEGpair.countmatrix.list)) {
    DriverDEGpair.countmatrix.i = DriverDEGpair.countmatrix.list[[i]]
    if (length(DriverDEGpair.countmatrix.i)==0) {
      next
    } else {
      drivers.i = rownames(DriverDEGpair.countmatrix.i)
      degs.i = colnames(DriverDEGpair.countmatrix.i)
      DriverDEGpair.countmatrix[drivers.i, degs.i] = DriverDEGpair.countmatrix[drivers.i, degs.i] + DriverDEGpair.countmatrix.i
    }
  }
  rm(DriverDEGpair.countmatrix.list)
  gc()

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

  #stop cluster
  stopCluster(cl)

  # Return the Driver.DEGlist
  return(Driver.DEGlist)
}











