##' NCBI Database API - Directly use NCBI EPost API
##'
##' NCBI EPost provide thousands of queries in one HTTP POST.
##' @title NCBI EPost API
##' @param postPara A named list of HTTP POST terms. For example, "db" is E-utility Database Name, see \url{http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly}.
##' @param ... Parameters inherited from the "postForm()" function in the "RCurl" package.
##' @return A names list containing "QueryKey", "Count", "WebEnv".
##' @examples
##' genePara <- list(db = "gene", id="948242,15486644")
##' genePost <- EPostNCBI(genePara)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl postForm
##' @references Entrez Programming Utilities Help \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/}
##' @export
##'
##' 
EPostNCBI <- function(postPara, ...) {
  
  ## NCBI EPost url
  urlBase <- EUrl('epost')
  
  ## HTTP POST with EPost
  postStr <- postForm(uri = urlBase, .params = postPara)
  postXml <- read_xml(postStr)

  ## retrieve key and webenv
  postPrefix <- c('/ePostResult/')
  postItems <- c('QueryKey', 'WebEnv')
  postInfo <- BatchXmlText(postXml, postPrefix, postItems)

  return(postInfo)
}
