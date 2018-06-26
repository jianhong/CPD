#' merge the results from multiple folders
#' @description merge multiple on.target.summary.csv into one file
#' @param output output file of summary file. Default is topNgRNA.csv.
#' @param path Path to be searched. Default is current folder.
#' @param filename file name to be searched. Default is on.target.summary.csv
#' @param n top N gRNAs per file. Default is 2.
#' @param cutoffEfficacy cutoff value of gRNA efficacy value. Default is 0
#' @param cutofftop5OfftargetTotalScore cutoff value of top5OfftargetTotalScore. Default is 5.
#' @param ... not used.
#' @return a dataframe
#' @export
#' 
topNgRNA <- function(output="topNgRNA.csv", path=".", filename="on.target.summary.csv",
                     n = 2, cutoffEfficacy = 0, cutofftop5OfftargetTotalScore=5,
                     ...){
  if(length(output)>0) stopifnot(is.character(output))
  stopifnot(is.character(filename))
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(cutoffEfficacy))
  stopifnot(is.numeric(cutofftop5OfftargetTotalScore))
  files <- dir(path = path, pattern = filename, recursive = TRUE, full.names = TRUE)
  features <- basename(dirname(files))
  names(files) <- features
  topTbl <- lapply(files, read.csv, nrows=n)
  topTbl <- lapply(topTbl, function(.ele){
    .ele[.ele$gRNAefficacy>=cutoffEfficacy & 
           .ele$top5OfftargetTotalScore<=cutofftop5OfftargetTotalScore, , drop=FALSE]
  })
  feas <- rep(features, sapply(topTbl, nrow))
  topTbl <- do.call(rbind, topTbl)
  topTbl <- cbind(features=feas, topTbl)
  if(length(output)>0) {
    tryCatch(write.csv(topTbl, output),
             error = function(e){ message(e) })
  }
  topTbl
}
