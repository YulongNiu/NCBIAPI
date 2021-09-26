##' Batch retrieve xml context
##'
##' Retrieve xml context by xPath.
##' @title Retrieve xml context
##' @param xmlObj xml object.
##' @param xPrefix xPath prefix.
##' @param xItems Node names.
##' @return A named list.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 xml_find_all xml_text
##' @importFrom magrittr %>%
##' @keywords internal
##'
##'
BatchXmlText <- function(xmlObj, xPrefix, xItems) {

  xPathBatch <- paste0(xPrefix, xItems)

  batchText <- lapply(xPathBatch, function(eachXPath) {
    xmlObj %>%
      xml_find_all(eachXPath) %>%
      xml_text
  })

  names(batchText) <- xItems

  return(batchText)
}
