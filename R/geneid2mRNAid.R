##' NCBI Database API - Get single NCBI gene information
##'
##' Get mRNA accession number.
##' @title Get mRNA accession number from gene id for one gene
##' @return NA or mRNA accession number
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @inheritParams singleGeneInfo
##' @importFrom xml2 xml_find_all
##' @keywords internal
##'
ExtractmRNAid <- function(geneXml) {
  ## 1: Gene-commentary_type value="mRNA
  ## 2. parent nodes
  ## 3. Gene-commentary_accession
  pat <- './/Gene-commentary_type[@value="mRNA"]/parent::*/Gene-commentary_accession'
  mRNAid <- unique(xml_text(xml_find_all(geneXml, pat)))

  return(mRNAid)
}



##' NCBI Database API - Get NCBI mRNA accession number from given NCBI gene IDs
##'
##' Get NCBI mRNA accession number
##' @title Get mRNA accession number
##' @return A list containing retrieve mRNA accession number
##' @examples
##' geneid2mRNAid(c('100286922', '948242', '15486644'), n = 2)
##'
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' ##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children
##' @importFrom foreach foreach %do% %dopar%
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom ParaMisc CutSeqEqu
##' @inheritParams getNCBIGenesInfo
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
geneid2mRNAid <- function(NCBIGeneIDs, n = 1, maxEach = 10000) {

  ## register multiple core
  registerDoParallel(cores = n)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~EPost~~~~~~~~~~~~~~~~~~~~~~~
  ## compress gene IDs
  geneIDs <- paste(NCBIGeneIDs, collapse = ',')
  infoPostPara <- list(db = 'gene', id = geneIDs)
  infoPost <- EPostNCBI(infoPostPara)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~~~~~~~EFetch~~~~~~~~~~~~~~~~~~~~~~~~~
  cutMat <- CutSeqEqu(length(NCBIGeneIDs), maxEach)
  ## The start number is from 0.
  cutMat <- cutMat - 1

  ## fetch url base
  fetchUrlBase <- EUrl('efetch')
  key = infoPost$QueryKey
  webEnv = infoPost$WebEnv

  geneInfo <- foreach (i = 1:ncol(cutMat), .combine = c) %do% {

    eachFetchStr <- postForm(uri = fetchUrlBase,
                             db = 'gene',
                             query_key = key,
                             WebEnv = webEnv,
                             retstart = cutMat[1, i],
                             retmax = maxEach,
                             retmode = 'xml')
    eachFetchXml <- read_xml(eachFetchStr)

    childXml <- xml_children(eachFetchXml)

    eachInfo <- foreach(j = 1 : length(childXml)) %dopar% {

      singleInfo <- ExtractmRNAid(childXml[[j]])

      return(singleInfo)
    }

    names(eachInfo) <- xml_text(xml_find_all(eachFetchXml, '//Gene-track_geneid/text()'))

    return(eachInfo)
  }

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## stop multiple core
  stopImplicitCluster()

  return(geneInfo)
}
