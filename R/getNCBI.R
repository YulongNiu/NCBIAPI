##' NCBI Database API - Get NCBI taxonomy information from given NCBI taxonomy IDs
##'
##' Get NCBI taxonomy information.
##' @title Get NCBI taxonomy information
##' @param NCBITaxoIDs A vector of NCBI taxonomy IDs.
##' @inheritParams getNCBIGenesInfo
##' @return A list containing taxonomy information for each ID.
##' @examples
##' ## with two cores
##' tax3 <- getNCBITaxo(c('9606', '511145', '797302'), n = 2)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children xml_text
##' @importFrom foreach foreach %do% %dopar%
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom ParaMisc CutSeqEqu
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##'
getNCBITaxo <- function(NCBITaxoIDs, n = 1, maxEach = 10000) {

  ## register multiple core
  registerDoParallel(cores = n)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~EPost~~~~~~~~~~~~~~~~~~~~~~~
  ## compress taxonomy IDs
  taxoIDs <- paste(NCBITaxoIDs, collapse = ',')
  infoPostPara <- list(db = 'taxonomy', id = taxoIDs)
  infoPost <- EPostNCBI(infoPostPara)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~~~~~~~ESummary~~~~~~~~~~~~~~~~~~~~~~~~~
  cutMat <- CutSeqEqu(length(NCBITaxoIDs), maxEach)
  ## The start number is from 0.
  cutMat <- cutMat - 1

  ## fetch url base
  fetchUrlBase <- EUrl('efetch')
  key = infoPost$QueryKey
  webEnv = infoPost$WebEnv

  taxoInfo <- foreach (i = 1:ncol(cutMat), .combine = c) %do% {
    eachFetchStr <- postForm(uri = fetchUrlBase,
                             db = 'taxonomy',
                             query_key = key,
                             WebEnv = webEnv,
                             retstart = cutMat[1, i],
                             retmax = maxEach,
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

  ## stop multiple core
  stopImplicitCluster()

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




##' NCBI Database API - Get NCBI gene or protein information from given NCBI gene IDs
##'
##' Get NCBI gene information, including gene name, description, genetic source, aliases, gene location. To retrieve thousands of proteins, use EPost to post record into the web server and then retrieve data using ESummary. If the gene ID is not found, return an error information in the list.
##' @title Get NCBI genes information
##' @param NCBIGeneIDs A vector of NCBI gene IDs.
##' @param n The number of CPUs or processors, and the default value is 1.
##' @param maxEach The maximum retrieve number in each visit. The ESearch, EFetch, and ESummary, the max number in one query is 10,000.
##' @return A list containing gene information for each ID. A empty character vector (whose length is 0) will be returned for the items if the contents are not found.
##' @examples
##' gene3 <- getNCBIGenesInfo(c('100286922', '948242', '15486644'), type = 'gene', n = 2)
##' protein2 <- getNCBIGenesInfo(c('MBF1669179', 'BAI64724'), type = 'protein', n = 2)
##' ## not found
##' ghostInfo <- getNCBIGenesInfo('111111111', n = 1)
##' \dontrun{
##' require(KEGGAPI)
##' ## signle genome with two plasmids
##' smuGenes <- convKEGG('smu', 'ncbi-geneid')
##' smuGeneNames <- sapply(strsplit(smuGenes[, 1], split = ':', fixed = TRUE), '[[', 2)
##' smuInfo <- getNCBIGenesInfo(smuGeneNames, n = 4)
##'
##' ## two genomes with two plasmids
##' draGenes <- convKEGG('dra', 'ncbi-geneid')
##' draGeneNames <- sapply(strsplit(draGenes[, 1], split = ':', fixed = TRUE), '[[', 2)
##' draInfo <- getNCBIGenesInfo(draGeneNames, n = 4)
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml xml_children
##' @importFrom foreach foreach %do% %dopar%
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom ParaMisc CutSeqEqu
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##'
getNCBIGenesInfo <- function(NCBIGeneIDs, type = 'gene', n = 1, maxEach = 10000) {

  ## register multiple core
  registerDoParallel(cores = n)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~EPost~~~~~~~~~~~~~~~~~~~~~~~
  ## compress gene IDs
  geneIDs <- paste(NCBIGeneIDs, collapse = ',')
  infoPostPara <- list(db = type, id = geneIDs)
  infoPost <- EPostNCBI(infoPostPara)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~~~~~~~ESummary~~~~~~~~~~~~~~~~~~~~~~~~~
  cutMat <- CutSeqEqu(length(NCBIGeneIDs), maxEach)
  ## The start number is from 0.
  cutMat <- cutMat - 1

  ## fetch url base
  fetchUrlBase <- EUrl('esummary')
  key = infoPost$QueryKey
  webEnv = infoPost$WebEnv

  geneInfo <- foreach (i = 1:ncol(cutMat), .combine = c) %do% {

    eachFetchStr <- postForm(uri = fetchUrlBase,
                             db = type,
                             query_key = key,
                             WebEnv = webEnv,
                             retstart = cutMat[1, i],
                             retmax = maxEach,
                             retmode = 'xml')

    eachFetchXml <- read_xml(eachFetchStr)
    topNode <- ifelse(type == 'gene', 'DocumentSummarySet/DocumentSummary', 'DocSum')
    childXml <- xml_find_all(eachFetchXml, topNode)

    eachInfo <- foreach(j = 1 : length(childXml)) %dopar% {

      singleInfo <- ifelse(type == 'gene',
                           singleGeneInfo(childXml[[j]]),
                           singleProteinInfo(childXml[[j]]))

      return(singleInfo)
    }

    return(eachInfo)
  }
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  names(geneInfo) <- NCBIGeneIDs

  ## stop multiple core
  stopImplicitCluster()

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


##' NCBI Database API - Get single NCBI protein information
##'
##' Get gene information form single NCBI protein ID.
##' @title Get single NCBI protein information
##' @param geneXml Gene xml data.
##' @return A list of protein information.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 xml_find_all xml_text
##' @importFrom magrittr %>%
##' @importFrom stringr str_replace
##' @keywords internal
##'
##'
singleProteinInfo <- function(proteinXml) {

  ## first check if the no candidate for input gene
  errorText <- proteinXml %>%
    xml_find_all('error') %>%
    xml_text

  if (length(errorText) > 0) {
    proteinInfo <- errorText
    return(proteinInfo)
  } else {}

  ## protein summary
  itemsAtts <- c('Caption', 'Title', 'Extra', 'Gi', 'CreateDate', 'UpdateDate', 'Flags', 'TaxId', 'Length', 'Status')

  proteinInfo <- sapply(itemsAtts, function(eachAttr) {
    eachAttr %>%
      str_replace('Item[@Name="Attrs"]', 'Attrs', .) %>%
      xml_find_all(proteinXml, .) %>%
      xml_text
  })

  return(proteinInfo)
}


##' NCBI Database API - Get single NCBI whole genomic gene annotation
##'
##' Get whole gene annotation form single NCBI genome ID. The locus tag is used as names for each gene. If one of the gene feature value is missed, a "" (with length of 1) will return. If the genome has no gene featurs, "NULL" will be returned.
##' This function now supports two feature types, "gene" or "CDS" (coding sequence). Other features such as RNAs ("ncRNA", "rRNA", "tRNA", "tmRNA"), "misc_feature", "rep_origin", "repeat_region" are not supported yet. It is shown in E. coli K-12 MG1655 "genes" features not only includes all "CDS" and RNAs, but some sRNA ("b4714"). "misc_feature" are mainly about cryptic prophage genes, and "repeat_region" are repetitive extragentic palindromic (REP) elements.
##' @title Get single NCBI whole genomic gene annotation
##' @param genomeID Single NCBI genome ID.
##' @param type "gene" or "CDS". The KEGG database use "CDS" as the protein gene count.
##' @inheritParams getNCBIGenesInfo
##' @return A list of annotation.
##' @examples
##' ## no gene features
##' nofeature <- singleGenomeAnno('BA000048')
##'
##' \dontrun{
##' aeuGenome <- singleGenomeAnno('CP007715', n = 4)
##'
##' ## missed value is replaced by ''
##' pseudo <- singleGenomeAnno('AE001826')[54]}
##' @importFrom RCurl postForm
##' @importFrom xml2 read_xml
##' @importFrom foreach foreach %dopar%
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
singleGenomeAnno <- function(genomeID, type = 'gene', n = 1) {


  getEachAnno <- function(featureNode) {
    ## USE: extact annotation from each node
    ## INPUT: `featureNode` is the child node in xml format
    ## OUTPUT: A list of gene annotation

    locNode <- xml_find_all(featureNode, 'GBFeature_intervals/GBInterval')
    loc <- BatchXmlText(locNode, './/', c('GBInterval_from', 'GBInterval_to'))
    loc <- do.call(cbind, loc)

    GBNodes <- xml_find_all(featureNode, 'GBFeature_quals/GBQualifier')
    GBf <- lapply(GBNodes, function(x) {
      ## value may be missed, for example a <psedudo>
      ## example: singleGenomeAnno('AE001826')[50]
      eachGB <- BatchXmlText(x, '', c('GBQualifier_name', 'GBQualifier_value'))
      ## '' may be assighed to multiple elements
      eachGB[which(sapply(eachGB, length) == 0)] <- ''
      eachGB <- unlist(eachGB)
      return(eachGB)
    })
    
    GBf <- do.call(rbind, GBf)
    GBf[, 2] <- gsub('\\n', '', GBf[, 2])
    
    geneAnno <- list(GBInterval = loc,
                     GBFeature_quals = GBf)

    return(geneAnno)
  }

  ## register multiple core
  registerDoParallel(cores = n)

  ##~~~~~~~~~~~~~~~~~~~load in whole genomic annotation~~~~~~~~~~~~
  urlBase <- EUrl('efetch')
  postList <- list(db = 'nuccore',
                   id = genomeID,
                   retmode = 'xml')

  annoStr <- postForm(urlBase, .params = postList)
  annoXml <- read_xml(annoStr)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## extract annotation node and keys
  annoNode <- xml_find_all(annoXml, 'GBSeq/GBSeq_feature-table')
  keys <- xml_text(xml_find_all(annoNode, './/GBFeature_key'))

  ## may be no gene features
  if (length(keys) == 1) {
    ## only one key that is "source", and return NULL
    annoList <- NULL
  } else {
    
    rightKeys <- which(keys == type)
    annoChild <- xml_children(annoNode)[rightKeys]
  
    ##~~~~~~~~~~~~~~~~~~extract features~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    annoList <- foreach(i = 1:length(annoChild)) %dopar% {
      eachAnno <- getEachAnno(annoChild[[i]])
    }
    
    locusName <- sapply(annoList, function(x) {
      eachLocus <- x[[2]]
      eachLocus <- eachLocus[eachLocus[, 1] == 'locus_tag', 2]
      return(eachLocus)
    })

    names(annoList) <- locusName
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  }

  ## stop multiple core
  stopImplicitCluster()

  return(annoList)
}
