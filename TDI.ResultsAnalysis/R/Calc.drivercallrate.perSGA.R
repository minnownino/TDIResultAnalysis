#' @title
#'
#' @description
#'
#' @details
#'
#' @param drivercallmatrix
#' @param SGAmatrix
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

Calc.drivercallrate.perSGA <- function(drivercallmatrix, SGAmatrix) {
  # Driver call matrix and SGA matrix
  tumors = rownames(drivercallmatrix)
  numoftumors = length(tumors)
  drivers = colnames(drivercallmatrix)
  SGAmatrix = SGAmatrix[tumors, ]
  SGAnames = colnames(SGAmatrix)
  numofSGAs = length(SGAnames)

  # Calculate the Driver call rate for each SGA
  drivercallrate = matrix(0, numofSGAs, 3)
  rownames(drivercallrate) = SGAnames
  colnames(drivercallrate) = c("#SGAEvents", "#DriverEvents", "Drivercallrate")

  # Count number of SGA events for each gene
  numoftumors.perSGA = colSums(SGAmatrix)
  drivercallrate[, "#SGAEvents"] = numoftumors.perSGA

  # Count number of driver events for each gene
  numoftumors.perdriver = colSums(drivercallmatrix)
  drivercallrate[intersect(drivers, SGAnames), "#DriverEvents"] = numoftumors.perdriver[intersect(drivers, SGAnames)]

  # Calculate driver call rate for each gene
  drivercallrate[, "Drivercallrate"] = drivercallrate[, "#DriverEvents"]/drivercallrate[, "#SGAEvents"]

  # Return the driver call rate matrix
  return(drivercallrate)
}






