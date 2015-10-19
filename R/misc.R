##' Cut vectors with invervals
##'
##' CutSeq() is used to cut a vector with different invervals. CutSeqEqu() is used to cut a vector with same invervals.
##' @title Cut vectors
##' @param cutSeq The inverals vector. The length of cutSeq could be more than 1, and 0 will be automatically excluded. 
##' @return A cut matrix in which the first row is the start point and second row is the end point.
##' @examples
##' # with one interval
##' CutSeq(10)
##' # with multiple interval
##' CutSeq(c(2, 3, 5))
##' # exclude 0
##' @rdname CutSeqInterval
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
##' 
CutSeq <- function(cutSeq){

  # remove 0, because we cannot cut a sequence by the internal of 0.
  cutSeq <- cutSeq[cutSeq != 0]
  vecCutseq <- length(cutSeq)

  if (vecCutseq == 1) {
    headCut <- 1
    endCut <- cutSeq
  } else {
    # loopCutSeq is the circle of vecCutseq
    loopCutSeq <- list()
    for(i in 1:vecCutseq) {
      loopCutSeq[[i]] <- cutSeq[1:i]
    }

    loopSumCutSeq <-  sapply(loopCutSeq, sum)

    # the head and tail sequence
    headCut <- c(1,loopSumCutSeq[1:(vecCutseq-1)]+1)
    endCut <- loopSumCutSeq
  }

  cutMat <- matrix(c(headCut, endCut), 2, byrow=TRUE)

  return(cutMat)

}


##' @param vecLen The length of vector used to cut.
##' @param equNum The equal internal.
##' @examples
##' # equal interval is the same as the length of vector
##' CutSeqEqu(10, equNum = 10)
##' CutSeqEqu(21, equNum = 10)
##' # euqal interval is larger than the length of vector
##' CutSeqEqu(10, equNum = 20)
##' @rdname CutSeqInterval
##' @export
##'
##' 
CutSeqEqu <- function(vecLen, equNum){
  
  if (equNum > vecLen){
    # the internal is bigger than the length of vecLen. So we use the full vecLen.
    cutMat <- matrix(c(1, vecLen))
  } else {
    timeNum <- vecLen %/% equNum
    remainer <- vecLen %% equNum
    cutSeq <- c(rep(equNum, timeNum), remainer)
    cutMat <- CutSeq(cutSeq)
  }

  return(cutMat)
}

 

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
##' @keywords internal
##'
##' 
BatchXmlText <- function(xmlObj, xPrefix, xItems) {
  
  xPathBatch <- paste0(xPrefix, xItems)

  batchText <- lapply(xPathBatch, function(eachXPath) {
    eachNodeSet <- xml_find_all(xmlObj, eachXPath)
    eachContent <- xml_text(eachNodeSet)
  })
  
  names(batchText) <- xItems

  return(batchText)
}
