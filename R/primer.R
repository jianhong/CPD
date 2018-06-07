#' primer designer
#' @description CRISPR/Cas9 primer designer.
#' @param gr A \link[GenomicRanges:GRanges-class]{GRanges} object or coordinates files could be converted into GRanges.
#' @param upstream,downstream numeric(1).upstream or downstream of the coordinates.
#' @param anchor Anchor point of upstream and downstream. TSS: search primter for promoter of given coordinates.
#' @param species character(1). Avaliable values: "mm9", "mm10", "hg19", "hg38", "danRer10"
#' @param ... Parameters could be passed to \link[CRISPRseek:offTargetAnalysis]{offTargetAnalysis}
#' @return NULL
#' @importFrom ChIPpeakAnno toGRanges 
#' @importFrom CRISPRseek offTargetAnalysis
#' @importFrom rtracklayer export.bed
#' @import GenomicRanges
#' @export
#' @examples 
#' gr <- system.file('extdata', 'Lamp3.bed', package = 'CPD')
#' primer(gr, anchor="onlyUpstream&Downstream", species="mm10")
#' 
primer <- function(gr, upstream=2000, downstream=2000, anchor=c("TSS", "TES", "ALL", "onlyUpstream&Downstream"),
                   species=c("mm10", "mm9", "hg19", "hg38", "danRer10"), ...){
  if(missing(gr)){
    stop("gr is required.")
  }
  species <- match.arg(species)
  anchor <- match.arg(anchor)
  if(!is(gr, "GRanges")){
    gr <- toGRanges(gr)
  }
  stopifnot(is(gr, "GRanges"))
  gr <- switch(anchor,
               TSS=promoters(x = gr, upstream = upstream, downstream = downstream),
               TES={
                 gr1 <- gr
                 gr.neg <- as.character(strand(gr))=="-"
                 strand(gr1[!gr.neg]) <- rep("-", length(gr[!gr.neg]))
                 strand(gr1[gr.neg]) <- rep("+", length(gr[gr.neg]))
                 gr1 <- promoters(x = gr, upstream = upstream, downstream = downstream)
                 stopifnot(length(gr1)==length(gr))
                 strand(gr1) <- strand(gr)
                 gr1
               },
               ALL={
                 suppressWarnings({
                   start(gr) <- start(gr) - upstream
                   end(gr) <- end(gr) + downstream
                   trim(gr)
                 })
               },
               "onlyUpstream&Downstream"={
                 suppressWarnings({
                   left <- right <- gr
                   strand(left) <- "+"
                   strand(right) <- "-"
                   left <- promoters(x = left, upstream = upstream, downstream = 0)
                   right <- promoters(x = right, upstream = downstream, downstream = 0)
                   strand(left) <- strand(right) <- strand(gr)
                   reduce(c(left, right))
                 })
               })
  inputFilePath <- tempfile()
  export.bed(object = gr, con = inputFilePath, format = "BED")
  bgn <- cbind(genome=c("BSgenome.Mmusculus.UCSC.mm10", 
                        "BSgenome.Mmusculus.UCSC.mm9",
                        "BSgenome.Hsapiens.UCSC.hg19",
                        "BSgenome.Hsapiens.UCSC.hg38",
                        "BSgenome.Drerio.UCSC.danRer10"),
               txdb=c("TxDb.Mmusculus.UCSC.mm10.knownGene",
                      "TxDb.Mmusculus.UCSC.mm9.knownGene",
                      "TxDb.Hsapiens.UCSC.hg19.knownGene",
                      "TxDb.Hsapiens.UCSC.hg38.knownGene",
                      "TxDb.Drerio.UCSC.danRer10.refGene"),
               orgAnn=c("org.Mm.egSYMBOL",
                     "org.Mm.egSYMBOL",
                     "org.Hs.egSYMBOL",
                     "org.Hs.egSYMBOL",
                     "org.Dr.egSYMBOL"))
  rownames(bgn) <- c("mm10", "mm9", "hg19", "hg38", "danRer10")
  args <- list(...)
  dots <- substitute(list(...))[-1]
  names <- make.names(unlist(sapply(dots, deparse)))
  names(args) <- names
  if(is.null(args$BSgenomeName)){
    require(bgn[species, "genome"], character.only = TRUE)
    args$BSgenomeName <- get(bgn[species, "genome"])
  }
  if(is.null(args$txdb)){
    require(bgn[species, "txdb"], character.only = TRUE)
    args$txdb <- get(bgn[species, "txdb"])
  }
  if(is.null(args$orgAnn)){
    require(sub("SYMBOL", ".db", bgn[species, "orgAnn"]), character.only = TRUE)
    args$orgAnn <- get(bgn[species, "orgAnn"])
  }
  if(is.null(args$outputDir)){
    args$outputDir <- getwd()
  }
  args$format <- "bed"
  args$inputFilePath <- inputFilePath
  do.call(what = offTargetAnalysis, args = args)
}