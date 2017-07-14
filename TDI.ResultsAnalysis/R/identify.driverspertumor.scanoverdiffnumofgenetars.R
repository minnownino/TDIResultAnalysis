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
#' @param fpath.TDIresults
#' @param fname.postprobcutoff
#' @param fname.output
#' @param numofDEGtars.2scan
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

identify.driverspertumor.scanoverdiffnumofgenetars <- function(fpath.TDIresults, fname.postprobcutoff, fname.output, numofDEGtars.2scan=1:15) {
  ### Go through all tumor samples
  # Posterior probability cutoff per SGA
  PostProb.cutoff.perSGA = as.matrix(read.table(fname.postprobcutoff, header=T, row.names=1))
  SGAnames = rownames(PostProb.cutoff.perSGA)

  # Get the list of tumor samples
  files = list.files(fpath.TDIresults)
  files = files[which(regexpr("TCGA", files)>0)]
  numoftumors = length(files)

  # Initialization of list to store the drivers per tumor
  driverspertumor = list()
  for (j in 1:length(numofDEGtars.2scan)) {
    driverspertumor[[j]] = list()
  }
  names(driverspertumor) = paste("numofDEGtars.cutoff=", numofDEGtars.2scan, sep="")

  # Go through all tumor samples
  cat("Start processing tumor samples from ", fpath.TDIresults, "\n")
  for (j in 1:numoftumors) {
    # tumor j
    tumor.j = files[j]
    cat("Processing sample.", j, ": ", tumor.j, "...\n", sep="")

    # Posterior probability for cancer type i and tumor j
    file2read.j = paste(fpath.TDIresults, "/", tumor.j, sep="")
    PostProb.j = read.data(file2read.j)
    PostProb.j[which(is.na(PostProb.j))] = 0 # Set NA values to 0
    SGAset.j = as.matrix(rownames(PostProb.j))
    numofSGAs.j = length(SGAset.j)

    # Find the SGA with maximal posterior probability for each DEG
    numofDEGs.j = ncol(PostProb.j)
    maxPostProbperDEG.j = matrix(0, numofDEGs.j, 3)
    colnames(maxPostProbperDEG.j) = c("idx.SGA", "idx.DEG", "PostProb")
    maxPostProbperDEG.j[, "idx.DEG"] = seq(1,numofDEGs.j,by=1)
    maxPostProbperDEG.j[, "idx.SGA"] = apply(PostProb.j, 2, which.max)
    maxPostProbperDEG.j[, "PostProb"] = apply(PostProb.j, 2, max)

    # Only consider SGAs with posterior probability cutoff, A0 and cancer type labels not included as well.
    indx2keep = which(is.element(SGAset.j[maxPostProbperDEG.j[, "idx.SGA"]], SGAnames)>0)
    if (length(indx2keep)==1) {
      maxPostProbperDEG.j = t(as.matrix(maxPostProbperDEG.j[indx2keep, ]))
    } else {
      maxPostProbperDEG.j = maxPostProbperDEG.j[indx2keep, ]
    }

    # Use the posterior probability cutoff to screen out the SGA-DEG pairs with smaller posterior probabilities
    PostProb.cutoff.j = as.matrix(PostProb.cutoff.perSGA[SGAset.j[maxPostProbperDEG.j[, "idx.SGA"]], 1])
    indx2keep = which(maxPostProbperDEG.j[, "PostProb"]>=PostProb.cutoff.j)
    if (length(indx2keep)==1) {
      maxPostProbperDEG.j = t(as.matrix(maxPostProbperDEG.j[indx2keep, ]))
    } else {
      maxPostProbperDEG.j = maxPostProbperDEG.j[indx2keep, ]
    }

    # Count the number of target DEGs for each SGA
    numoftarDEGs.perSGA.j = matrix(0, numofSGAs.j, 1)
    rownames(numoftarDEGs.perSGA.j) = SGAset.j
    colnames(numoftarDEGs.perSGA.j) = "numoftarDEGs"
    for (k in 1:numofSGAs.j) {
      numoftarDEGs.perSGA.j[SGAset.j[k],] = sum(is.element(SGAset.j[maxPostProbperDEG.j[,"idx.SGA"]], SGAset.j[k]))
    }

    # Find the drivers for different number of target DEGs cutoff
    for (numofDEGtars.k in numofDEGtars.2scan) {
      driver.k = rownames(numoftarDEGs.perSGA.j)[which(numoftarDEGs.perSGA.j>=numofDEGtars.k)]
      driverspertumor[[paste("numofDEGtars.cutoff=", numofDEGtars.k, sep="")]][[gsub(".csv", "", tumor.j)]] = driver.k
    }
  }

  # Save the driverspertumor list to .RData for cancer i
  save(driverspertumor, file=fname.output)
}

