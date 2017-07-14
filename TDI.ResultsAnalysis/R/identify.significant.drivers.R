#' @title Identify Significant Drivers
#'
#' @description A driver is called as a significant driver if it satisfy the following conditions: 1. number of driver events is larger or equal to the given threshold, and 2. driver call rate is larger or equal to the given threshold.
#'
#' @param drivercallmatrix driver call per tumor
#' @param SGAmatrix SGA events per tumor
#' @param numofdrivercalls.thd Minimal number of driver events
#' @param drivercallrate.thd driver call rate threshold
#'
#' @export
#'
#' @return sigdrivers
#'
#' @examples \dontrun{
#'
#' }
#'

identify.significant.drivers <- function(drivercallmatrix, SGAmatrix, numofdrivercalls.thd=30, drivercallrate.thd=0.5) {
  ### Identify the significant drivers
  # Driver call matrix
  tumors = rownames(drivercallmatrix)
  numoftumors = length(tumors)
  SGAmatrix = SGAmatrix[tumors, ]
  genes.olp = intersect(colnames(drivercallmatrix), colnames(SGAmatrix))
  numofgenes.olp = length(genes.olp)
  drivercallmatrix = drivercallmatrix[, genes.olp]
  SGAmatrix = SGAmatrix[, genes.olp]

  # Count the number of tumors per SGA/driver.
  numoftumors.perSGA = as.matrix(colSums(SGAmatrix))
  numoftumors.perdriver = as.matrix(colSums(drivercallmatrix))
  drivercall.perSGA = matrix(0, numofgenes.olp, 3)
  rownames(drivercall.perSGA) = genes.olp
  colnames(drivercall.perSGA) = c("#SGAEvents", "#DriverEvents", "Drivercallrate")
  drivercall.perSGA[, "#SGAEvents"] = numoftumors.perSGA
  drivercall.perSGA[, "#DriverEvents"] = numoftumors.perdriver
  drivercall.perSGA[, "Drivercallrate"] = drivercall.perSGA[, "#DriverEvents"]/drivercall.perSGA[, "#SGAEvents"]

  # Keep the drivers that at least cover certain number/percentage of all tumors
  idx2keep = intersect(which(drivercall.perSGA[, "#DriverEvents"]>=numofdrivercalls.thd), which(drivercall.perSGA[, "Drivercallrate"]>=drivercallrate.thd))
  if (length(idx2keep)==0) {
    stop("No drivers have been identified under the current parameter set for", cancer.i, "\n")
  }

  # Significant drivers
  sigdrivers = drivercall.perSGA[idx2keep, ]

  # Return the significant driver list
  return(sigdrivers)
}

