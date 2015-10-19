##' NCBI Database API - Get NCBI taxonomy information from given NCBI taxonomy IDs
##'
##' Get NCBI taxonomy information. 
##' @title Get NCBI taxonomy information
##' @param NCBITaxoIDs A vector of NCBI taxonomy IDs.
##' @return A list containing taxonomy information for each ID.
##' @examples
##' threeTax <- getNCBITaxo(c('9606', '511145', '797302'))
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children xml_text
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##' 
getNCBITaxo <- function(NCBITaxoIDs) {

  ## compress taxonomy IDs
  taxoIDs <- paste(NCBITaxoIDs, collapse = ',')

  ## read in xml content with HTTP POST method
  urlBase <- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
  xmlStr <- postForm(uri = urlBase,
                     db = 'taxonomy',
                     id = taxoIDs,
                     retmode = 'xml')
  taxXML <- read_xml(xmlStr)

  ## extract child node of each tax
  taxChildren <- xml_children(taxXML)
  taxList <- lapply(taxChildren, function(x){

    ## all TaxId, all ScientificName, and all Rank
    xPathVec <- c('.//TaxId', './/ScientificName', './/Rank')
    allTax <- sapply(xPathVec, function(eachXPath) {
      eachNodeSet <- xml_find_all(x, eachXPath)
      eachContent <- xml_text(eachNodeSet)
    })

    ## move input tax to the last
    allTax <- rbind(allTax[-1, ], allTax[1, ])
    colnames(allTax) <- c('TaxId', 'ScientificName', 'Rank')

    return(allTax)
  })

  names(taxList) <- NCBITaxoIDs

  return(taxList)
}


##' NCBI Database API - Get NCBI gene information from given NCBI gene IDs
##'
##' Get NCBI gene information, including gene name, description, genetic source, aliases, gene location. To retrieve thousands of proteins, use EPost to post record into the web server and then retrieve data using ESummary. If the gene ID is not found, return an error information in the list.
##' @title Get NCBI genes information
##' @param NCBIGeneIDs A vector of NCBI gene IDs.
##' @param n The number of CPUs or processors, and the default value is 4.
##' @return A list containing gene information for each ID. A empty character vector (whose length is 0) will be returned for the items if the contents are not found.
##' @examples
##' gene3 <- getNCBIGenesInfo(c('100286922', '948242', '15486644'), n = 2)
##' ## not found
##' ghostInfo <- getNCBIGenesInfo('111111111', n = 1)
##' \dontrun{
##' require(KEGGAPI)
##' ## signle genome
##' smuGenes <- convKEGG('smu', 'ncbi-geneid')
##' smuGeneNames <- sapply(strsplit(smuGenes[, 1], split = ':', fixed = TRUE), '[[', 2)
##' smuInfo <- getNCBIGenesInfo(smuGeneNames)
##' ## multiple genomes
##' draGenes <- convKEGG('dra', 'ncbi-geneid')
##' draGeneNames <- sapply(strsplit(draGenes[, 1], split = ':', fixed = TRUE), '[[', 2)
##' draInfo <- getNCBIGenesInfo(draGeneNames)
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children
##' @importFrom foreach foreach %do% %dopar%
##' @importFrom doMC registerDoMC
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##' 
getNCBIGenesInfo <- function(NCBIGeneIDs, n = 4) {

  registerDoMC(n)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~EPost~~~~~~~~~~~~~~~~~~~~~~~
  ## compress gene IDs
  geneIDs <- paste(NCBIGeneIDs, collapse = ',')
  infoPostPara <- list(db = 'gene', id = geneIDs)
  infoPost <- EPostNCBI(infoPostPara)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~~~~~~~ESummary~~~~~~~~~~~~~~~~~~~~~~~~~
  ## For EFetch, the max number in one query is 10,000. The input gene IDs is splited every 500.
  cutMat <- CutSeqEqu(length(NCBIGeneIDs), 500)
  ## The start number is from 0.
  cutMat <- cutMat - 1

  ## fetch url base
  fetchUrlBase <- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
  key = infoPost$QueryKey
  webEnv = infoPost$WebEnv
  
  geneInfo <- foreach (i = 1:ncol(cutMat), .combine = c) %do% {
    
    eachFetchStr <-  postForm(uri = fetchUrlBase,
                              db = 'gene',
                              query_key = key,
                              WebEnv = webEnv,
                              retstart = cutMat[1, i],
                              retmax = 500,
                              retmode = 'xml')
    eachFetchXml <- read_xml(eachFetchStr)
    childXml <- xml_find_all(eachFetchXml, './/DocumentSummary')

    eachInfo <- foreach(j = 1 : length(childXml)) %dopar% {
      
      singleInfo <- singleGeneInfo(childXml[[j]])

      return(singleInfo)
    }

    return(eachInfo)
  }
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  names(geneInfo) <- NCBIGeneIDs
  
  return(geneInfo)
}


##' NCBI Database API - Get single NCBI gene information
##'
##' Get gene information form single NCBI gene ID.
##' @title Get single NCBI gene information
##' @param geneXml xml data.
##' @return A list of gene information.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 xml_find_all xml_text
##' @keywords internal
##'
##' 
singleGeneInfo <- function(geneXml) {

  ## first check if the no candidate for input gene
  errorChild <- xml_find_all(geneXml, 'error')
  errorNum <- length(errorChild)
  
  if (errorNum == 0) {
    ## gene summary
    docSumPrefix <- ''
    docSumItems <- c('Name', 'Description', 'Chromosome', 'GeneticSource', 'MapLocation', 'OtherAliases')
    geneInfo <- BatchXmlText(geneXml, docSumPrefix, docSumItems)

    ## gene location
    ## LocationHist also includes gene location which is not what we want
    locPrefix <- 'GenomicInfo/GenomicInfoType/'
    locItems <- c('ChrLoc', 'ChrAccVer', 'ChrStart', 'ChrStop', 'ExonCount')
    locText <- BatchXmlText(geneXml, locPrefix, locItems)
    locMat <- t(do.call(rbind, locText))

    ## combine summary and gene location
    geneInfo$GenomicInfo = locMat
  }
  else if (errorNum > 0) {
    ## return error info
    geneInfo <- xml_text(errorChild)
  }

  return(geneInfo)
}
