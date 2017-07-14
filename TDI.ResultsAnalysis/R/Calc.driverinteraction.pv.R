#' @title
#'
#' @description
#'
#'
#' @param sigdrivers.tarDEGs target DEGs of the significant drivers
#' @param numofgenes.total
#' @param log10pv.cutoff
#' @param numoftarDEGoverlaps.cutoff
#'
#' @export
#'
#' @return driverlinks
#'
#' @examples \dontrun{
#'
#' }
#'

Calc.driverinteraction.pv <- function(sigdriver.tarDEGs, numofgenes.total = 20000, log10pv.cutoff = -3, numoftarDEGoverlaps.cutoff = 3) {
  # Driver names
  driverlist = names(sigdrivers.tarDEGs)
  numofdrivers = length(driverlist)

  # Hyper-geometric test for each pair of driver based on their gene targets
  drivers.pairs = t(combn(driverlist, 2))
  numofdriverpairs = nrow(drivers.pairs)
  driverpairs.links = matrix(NA, numofdriverpairs, 6)
  colnames(driverpairs.links) = c("Driver1", "Driver2", "# driver1 DEG targets", "# driver2 DEG targets", "# DEG target overlaps", "log10PV (hypergeometric test)")
  driverpairs.links[, c("Driver1", "Driver2")] = drivers.pairs
  for (i in 1:numofdriverpairs) {
    # Number of target genes for driver 1 and driver 2, and their overlaps
    driver1.i = driverpairs.links[i, "Driver1"]
    driver2.i = driverpairs.links[i, "Driver2"]
    driver1tars.i = unique(names(sigdrivers.tarDEGs[[driver1.i]]))
    driver2tars.i = unique(names(sigdrivers.tarDEGs[[driver2.i]]))
    numofdriver1tars.i = length(driver1tars.i)
    numofdriver2tars.i = length(driver2tars.i)
    numofoverlaptars.i = length(intersect(driver1tars.i, driver2tars.i))

    # Hyper-geometric test
    log10.pv.i = log10(phyper(numofoverlaptars.i, numofdriver1tars.i, numofgenes.total-numofdriver1tars.i, numofdriver2tars.i, lower.tail=F) + dhyper(numofoverlaptars.i, numofdriver1tars.i, numofgenes.total-numofdriver1tars.i, numofdriver2tars.i))

    # Update the driverpairs.links matrix
    driverpairs.links[i, "# driver1 DEG targets"] = numofdriver1tars.i
    driverpairs.links[i, "# driver2 DEG targets"] = numofdriver2tars.i
    driverpairs.links[i, "# DEG target overlaps"] = numofoverlaptars.i
    driverpairs.links[i, "log10PV (hypergeometric test)"] = log10.pv.i
  }

  # Select the most significant driver-driver interactions
  indxsigdriverpairs = intersect(which(as.numeric(driverpairs.links[, "log10PV (hypergeometric test)"])<=log10pv.cutoff), which(as.numeric(driverpairs.links[, "# DEG target overlaps"])>=numoftarDEGoverlaps.cutoff))
  seldriverpair.links = driverpairs.links[indxsigdriverpairs, ]
  indx.order = sort.int(as.numeric(seldriverpair.links[, "log10PV (hypergeometric test)"]), index.return=T)$ix
  seldriverpair.links = seldriverpair.links[indx.order, ]

  # Return the significant interacting driver pairs
  return(seldriverpair.links)
}









