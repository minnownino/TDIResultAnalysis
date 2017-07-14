#' @title Convert a List of Gene Symbols to Standard NCBI Gene symbols
#'
#' @description Convert the given gene symbols to standard NCBI gene symbols by comparing with the NCBI gene id, symbol and synonyms.
#'
#' @details Each NCBI gene id has a corresponding gene symbol and may have one or multiple gene synonyms. Given the gene symbol, the function first check whether it is a standard NCBI gene symbol. If the gene symbol is not found in NCBI standard gene symbol list, the function further searches the gene synonyms list and return the standard NCBI gene symbol. A gene symbol will be set to NA, if it has no match or multiple match with the NCBI gene symbols.
#'
#'
#' @param gnames2convert A vector of gene symbols to be converted to standard NCBI symbols.
#'
#' @param geneidsymbolsynonyms.list Mapping between gene NCBI ids and corresponding gene symbols and symnonms. It consists of three lists, i.e. 1. gene NCBI id, 2. gene symbol, and 3. gene symnonms.
#'
#' @export
#'
#'
#' @return A vector which contains the standard NCBI gene symbols.
#'
#'
#' @examples \dontrun{
#'
#' }
#'
geneNCBIsymbol.conversion <- function(gnames2convert, geneidsymbolsynonyms.list) {
  # Standard NCBI gene symbols and synonyms
  genesymbols = geneidsymbolsynonyms.list$symbols
  genesynonyms = geneidsymbolsynonyms.list$synonyms
  numofgenesymbols = length(genesymbols)

  # Go through all genes
  numofgenes = length(gnames2convert)
  gnamesconverted = matrix(NA, numofgenes, 1)
  # Genes that are standard NCBI symbols
  idx.symbols = which(is.element(gnames2convert, genesymbols))
  gnamesconverted[idx.symbols] = gnames2convert[idx.symbols]
  for (i in 1:numofgenes) {
    if (!is.na(gnamesconverted[i])) {
      next # Standard gene symbol
    }

    gname2convert.i = gnames2convert[i]
    for (j in 1:numofgenesymbols) {
      genesynonyms.j = genesynonyms[[j]]
      if (is.element(gname2convert.i, genesynonyms.j)) {
        gnamesconverted[i] = genesymbols[j]
        break
      }
    }
  }

  return(gnamesconverted)
}

