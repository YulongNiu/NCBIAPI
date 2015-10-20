##' Generate E-utilities url
##'
##' A wrap function to generate E-utilities url. For now, it supports "einfo", "esearch", "epost", "esummary", "efetch", "elink", "egquery", "espell" , and "ecitmatch".
##' @title E-utilities url
##' @param eutils The nine E-utilities tools in lower cases.
##' @return url
##' @examples
##' efetchUrl <- EUrl("efetch")
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
EUrl <- function(eutils) {
  
  wrapUrl <- paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/', eutils, '.fcgi')

  return(wrapUrl)
}
