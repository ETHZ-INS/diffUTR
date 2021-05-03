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
#' @param tmpDir Passed to `Rsubread::featureCounts`
#' @param ... Passed to `Rsubread::featureCounts`
#'
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment-class}
#' @import SummarizedExperiment GenomicRanges
#' @importFrom Rsubread featureCounts
#' @importFrom edgeR DGEList cpm calcNormFactors
#' @export
#' @examples
#' data("example_gene_annotation", package="diffUTR")
#' bins <- prepareBins(example_gene_annotation)
#' bam_files <- list.files(system.file("extdata", package="diffUTR"),
#'                         pattern="bam$", full=TRUE)
#' se <- countFeatures(bam_files, bins, verbose=FALSE)
#' se
countFeatures <- function(bamfiles, bins, strandSpecific=0, readLength=50L,
                          allowMultiOverlap=TRUE, inclNormalized=TRUE,
                          tmpDir=tempdir(), ...){
  if(is.character(bins) && length(bins)==1 && grepl("\\.rds$", bins)){
    bins <- readRDS(bins)
  }
  if(!is(bins, "GRanges")) stop("Bins have to be provided as a GRanges")

  if(is.null(names(bins))) names(bins) <- paste(bins$gene,bins$bin_id,sep=".")
  binsframe<-as.data.frame(bins)
  names(binsframe)[names(binsframe) == "seqnames"] <- "Chr"
  binsframe$GeneID <- names(bins)

  if(is.null(names(bamfiles))) names(bamfiles) <- .cleanNames(bamfiles)

  hits <- Rsubread::featureCounts( bamfiles, ..., tmpDir=tmpDir,
                                   annot.ext=binsframe,
                                   isGTFAnnotationFile=FALSE,
                                   strandSpecific=strandSpecific,
                                   allowMultiOverlap=allowMultiOverlap,
                                   useMetaFeatures=FALSE )
  se <- SummarizedExperiment(list(counts=as.matrix(hits$counts)),
                             rowRanges=bins)
  colnames(se) <- names(bamfiles)
  wi <- pmax(width(se),readLength)
  if(ncol(se)==1){
    assays(se)$logcpm <- log1p(assay(se)*10^6/sum(as.numeric(assay(se))))
  }else{
    assays(se)$logcpm <- log1p(cpm(calcNormFactors(DGEList(assay(se)))))
  }
  assays(se)$logNormDensity <- log1p(exp(assays(se)$logcpm)/wi)
  rowData(se)$meanLogCPM <- rowMeans(assays(se)$logcpm)
  rowData(se)$logWidth <- log1p(width(se))
  rowData(se)$meanLogDensity <- rowMeans(assays(se)$logNormDensity)
  rowData(se)$logDensityRatio <- 
    log(.getDensityRatio(rowData(se)$meanLogDensity, rowData(se)$gene))
  colData(se)$assigned <- as.numeric(hits$stat[1,-1,drop=FALSE])
  colData(se)$unassigned <- colSums(hits$stat[-1,-1,drop=FALSE])
  if(!inclNormalized) assays(se) <- assays(se)[1]
  se
}
