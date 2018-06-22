#' filter offTargetAnalysis results
#' @description filter offTargetAnalysis results by sorting the results with best efficacy and shortest off target list.
#' @param results output of offTargetAnalysis
#' @param outputDir output folder of oofTargetAnalysis
#' @param ... not used
#' @return a sorted gRNAs table with efficay, off-target number, and unique restriction enzyme cut sites.
#' 

filterRes <- function(results, outputDir, ...){
  stopifnot(all(c("on.target", "summary") %in% names(results)))
  on.target <- results$on.target
  summary <- results$summary
  on.target.summary <- merge(on.target, summary, all.x=TRUE)
  on.target.summary$top10OfftargetTotalScore[is.na(on.target.summary$top10OfftargetTotalScore)] <- 0
  on.target.summary$top5OfftargetTotalScore[is.na(on.target.summary$top5OfftargetTotalScore)] <- 0
  on.target.summary$gRNAefficacy <- as.numeric(as.character(on.target.summary$gRNAefficacy))
  on.target.summary <- on.target.summary[with(on.target.summary, 
                                              order(top5OfftargetTotalScore, 
                                                    top10OfftargetTotalScore,
                                                    -gRNAefficacy)), ]
  strand <- ifelse(grepl("f$", as.character(on.target.summary$names)), "+", "-")
  on.target.summary <- cbind(on.target.summary[, 1:2], strand, on.target.summary[, -(1:2)])
  write.csv(on.target.summary, file.path(outputDir, "on.target.summary.csv"), row.names = FALSE)
  results$on.target <- on.target.summary
  results
}
