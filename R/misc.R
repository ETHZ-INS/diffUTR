#' addNormalizedAssays
#'
#' @param se A bin-wise `SummarizedExperiment` as produced by
#' \code{\link{countFeatures}}
#'
#' @return The `se` object with populated `logcpm` and `logNormDensity` assays.
#' @param readLength Used as a minimum width to estimate read density
#' (default 50).
#' @export
#'
#' @importFrom edgeR DGEList cpm calcNormFactors
#' @examples
#' data(example_bin_se)
#' example_bin_se <- addNormalizedAssays(example_bin_se)
addNormalizedAssays <- function(se, readLength=50L){
  cpm <- cpm(calcNormFactors(DGEList(assay(se))))
  assays(se)$logcpm <- log1p(cpm)
  assays(se)$logNormDensity <- log1p(1000*cpm/pmax(width(se),readLength))
  se
}

#' @importFrom stringi stri_reverse
.cleanNames <- function(x){
  x <- gsub(" ","_",gsub("-","_",gsub("/",".",x,fixed=TRUE),fixed=TRUE),
            fixed=TRUE)
  x <- rmLCS(rmLCS(x,"."),"_")
  stri_reverse(rmLCS(rmLCS(stri_reverse(x),"."),"_"))
}


.matchGene <- function(se, x){
  w <- which(rowData(se)$gene == x)
  if(length(w)==0) w <- which(any(rowData(se)$gene_name == x))
  w
}


# remove the longest common string at the beginning of all elements of a
# character vector
rmLCS <- function(x, delim=""){
  tmp <- strsplit(as.character(x),delim,fixed=TRUE)
  if(any(lengths(tmp)==1)) return(x)
  i <- 1
  while(length(unique(vapply(tmp, FUN.VALUE=character(1),
                             FUN=function(x) x[i])))==1){
    i <- i+1
  }
  if(i==1) return(x)
  vapply(tmp, FUN.VALUE=character(1), FUN=function(x){
    x <- x[!is.null(x)]
    paste(x[-seq_len(i-1)], collapse=delim)
  })
}

.typeColors <- function(){
  c("3UTR"="#117733", "CDS"="#332288", "CDS/3UTR"="#44AA99",
    "CDS/UTR"="#44AA99", "CDS/UTR/3UTR"="#44AA99", "UTR"="#DDCC77",
    "UTR/3UTR"="#999933", "non-coding"="#CC6677")
}

#' @importFrom methods is
.checkSE <- function(se, checkNorm=FALSE, requireStats=FALSE){
  stopifnot(is(se,"RangedSummarizedExperiment"))
  stopifnot(all(c("type","meanLogCPM","logWidth","meanLogDensity","gene") %in%
                  colnames(rowData(se))))
  if(requireStats){
    if( !all(c("bin.p.value","bin.FDR") %in% colnames(rowData(se))) ||
        is.null(gl <- metadata(se)$geneLevel) ||
        !(is.data.frame(gl) || is(gl, "DFrame")))
      stop("The object does not contain differential bin usage statistics. ",
           "You should run a DEU wrapper first (see `?DEUwrappers`).")
  }
  if(checkNorm){
    if(!all(c("logcpm","logNormDensity") %in% assayNames(se))){
      message("Computing normalized assays...")
      message("To avoid doing this again for every plot, run:\n",
              "se <- addNormalizedAssays(se)")
      se <- addNormalizedAssays(se)
    }
  }
  se
}