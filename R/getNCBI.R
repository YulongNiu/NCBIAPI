##' NCBI Database API - Get NCBI taxonomy information from given NCBI taxonomy IDs
##'
##' Get NCBI taxonomy information. 
##' @title Get NCBI taxonomy information
##' @param NCBITaxoIDs A vector of NCBI taxonomy IDs.
##' @inheritParams getNCBIGenesInfo
##' @return A list containing taxonomy information for each ID.
##' @examples
##' tax3 <- getNCBITaxo(c('9606', '511145', '797302'), n = 2)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children xml_text
##' @importFrom foreach foreach %do% %dopar%
##' @importFrom doMC registerDoMC
##' @importFrom ParaMisc CutSeqEqu
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##' 
getNCBITaxo <- function(NCBITaxoIDs, n = 4) {

  ## multiple core
  registerDoMC(n)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~EPost~~~~~~~~~~~~~~~~~~~~~~~
  ## compress taxonomy IDs
  taxoIDs <- paste(NCBITaxoIDs, collapse = ',')
  infoPostPara <- list(db = 'taxonomy', id = taxoIDs)
  infoPost <- EPostNCBI(infoPostPara)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~~~~~~~ESummary~~~~~~~~~~~~~~~~~~~~~~~~~
  ## For EFetch, the max number in one query is 10,000. The input gene IDs is splited every 500.
  cutMat <- CutSeqEqu(length(NCBITaxoIDs), 500)
  ## The start number is from 0.
  cutMat <- cutMat - 1

  ## fetch url base
  fetchUrlBase <- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
  key = infoPost$QueryKey
  webEnv = infoPost$WebEnv

  taxoInfo <- foreach (i = 1:ncol(cutMat), .combine = c) %do% {
    eachFetchStr <-  postForm(uri = fetchUrlBase,
                              db = 'taxonomy',
                              query_key = key,
                              WebEnv = webEnv,
                              retstart = cutMat[1, i],
                              retmax = 500,
                              retmode = 'xml')
    eachFetchXml <- read_xml(eachFetchStr)
    childXml <- xml_find_all(eachFetchXml, 'Taxon')

    eachInfo <- foreach(j = 1 : length(childXml)) %dopar% {
      
      singleInfo <- singleTaxoInfo(childXml[[j]])

      return(singleInfo)
    }
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return(eachInfo)
  }
  
  names(taxoInfo) <- NCBITaxoIDs
  
  return(taxoInfo)
}


##' NCBI Database API - Get NCBI taxonomy information from given NCBI taxonomy IDs
##'
##' Get taxonomy information form single NCBI taxonomy ID.
##' @title Get single NCBI taxonomy information
##' @param taxoXml Taxonomy xml data
##' @return A matrix of taxonomy information
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##'
##' 
singleTaxoInfo <- function(taxoXml) {

  taxoPrefix <- './/'
  taxoItems <- c('TaxId', 'ScientificName', 'Rank')
  taxoInfo <- BatchXmlText(taxoXml, taxoPrefix, taxoItems)

  taxoMat <- do.call(cbind, taxoInfo)

  return(taxoMat)
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
##' @importFrom ParaMisc CutSeqEqu
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##' 
getNCBIGenesInfo <- function(NCBIGeneIDs, n = 4) {

  ## multiple core
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
    childXml <- xml_find_all(eachFetchXml, 'DocumentSummarySet/DocumentSummary')

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
##' @param geneXml Gene xml data.
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
    locMat <- do.call(cbind, locText)

    ## combine summary and gene location
    geneInfo$GenomicInfo = locMat
  }
  else if (errorNum > 0) {
    ## return error info
    geneInfo <- xml_text(errorChild)
  }

  return(geneInfo)
}
