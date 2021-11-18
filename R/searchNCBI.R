##' NCBI Database API - Search single NCBI locus ID
##'
##' Search single NCBI locus ID and retrieve its.
##'
##' @title Search NCBI locus ID
##' @param NCBILocusID A \code{string} of NCBI locus ID.
##' @param type A \code{string} of database.
##' @return A \code{string} of NCBI gene ID.
##' @examples
##' searchNCBILocus('G4B84_011174')
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom xml2 read_xml xml_find_all xml_text
##' @export
##'
searchNCBILocus <- function(NCBILocusID, type = 'gene') {

  ##~~~~~~~~~~~~~~~~~~~~~~~~ESearch~~~~~~~~~~~~~~~~~~~~~~~~
  ## fetch url base
  fetchUrlBase <- EUrl('esearch')

  searchXML <- postForm(uri = fetchUrlBase,
                        db = type,
                        term = NCBILocusID,
                        retmode = 'xml') %>%
    read_xml

  geneID <- xml_find_all(searchXML, 'IdList/Id') %>%
    xml_text

  return(geneID)
}



