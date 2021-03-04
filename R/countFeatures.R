#' countFeatures
#'
#' @param bamfiles A vector of paths to bam files
#' @param bins A GRanges of bins in which to count reads (or path to a rds file
#' containing such an object
#' @param strandSpecific Passed to `Rsubread::featureCounts`
#' @param readLength Used as a minimum width to estimate read density.
#' @param allowMultiOverlap Passed to `Rsubread::featureCounts`
#' @param inclNormalized Logical; whether to include normalized assays (needed
#' for plotting)
#' @param ... Passed to `Rsubread::featureCounts`
#'
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment-class}
#' @import SummarizedExperiment GenomicRanges
#' @importFrom Rsubread featureCounts
#' @importFrom edgeR DGEList cpm calcNormFactors
#' @export
countFeatures <- function(bamfiles, bins, strandSpecific=1, readLength=50L,
                          allowMultiOverlap=TRUE, inclNormalized=TRUE, ...){
  if(is.character(bins) && length(bins)==1 && grepl("\\.rds$", bins)){
    bins <- readRDS(bins)
  }
  if(!is(bins, "GRanges")) stop("Bins have to be provided as a GRanges")

  if(is.null(names(bins))) names(bins) <- paste(bins$gene,bins$bin_id,sep=".")
  binsframe<-as.data.frame(bins)
  names(binsframe)[names(binsframe) == "seqnames"] <- "Chr"
  binsframe$GeneID <- names(bins)

  if(is.null(names(bamfiles))) names(bamfiles) <- .cleanNames(bamfiles)

  hits <- Rsubread::featureCounts( bamfiles, ...,
                                   annot.ext=binsframe,
                                   isGTFAnnotationFile=FALSE,
                                   strandSpecific=strandSpecific,
                                   allowMultiOverlap=allowMultiOverlap,
                                   useMetaFeatures=FALSE )
  se <- SummarizedExperiment(list(counts=as.matrix(hits$counts)),rowRanges=bins)
  colnames(se) <- names(bamfiles)
  wi <- pmax(width(se),readLength)
  assays(se)$logcpm <- log1p(cpm(calcNormFactors(DGEList(assay(se)))))
  assays(se)$logNormDensity <- log1p(exp(assays(se)$logcpm)/wi)
  rowData(se)$meanLogCPM <- rowMeans(assays(se)$logcpm)
  rowData(se)$logWidth <- log1p(width(se))
  rowData(se)$meanLogDensity <- rowMeans(assays(se)$logNormDensity)
  colData(se)$assigned <- as.numeric(hits$stat[1,-1])
  colData(se)$unassigned <- colSums(hits$stat[-1,-1])
  if(!inclNormalized) assays(se) <- assays(se)[1]
  se
}
