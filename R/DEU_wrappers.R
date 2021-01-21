#' DEUwrappers
#'
#' Wrappers around commonly-used DEU methods (\code{\link[edgeR]{diffSpliceDGE}},
#' \code{\link[DEXSeq]{DEXSeq}} and an improved version of \code{\link[limma]{diffSplice}}
#'
#' @param se A bin-wise SummarizedExperiment as produced by \code{\link{countFeatures}}
#' @param design A formula (using columns of `colData(se)`) or (for `diffSplice.wrapper`
#' or `diffSpliceDGE.wrapper` only) a model.matrix.
#' @param reducedModel A reduced formula (applicable only to `DEXSeq.wrapper`).
#' @param coef The coefficient to be tested (ignored for `DEXSeq.wrapper`).
#' @param QLF Logical; whether to use edgeR's quasi-likelihood negative binomial
#' (applicable only to `diffSpliceDGE.wrapper`).
#' @param robust Logical; whether to use robust fitting for the dispersion trend
#' (ignored for `DEXSeq.wrapper`).
#' @param countFilter Logical; whether to filter out low-count bins (ignored for
#' `DEXSeq.wrapper`).
#' @param excludeTypes A vector of bin types to ignore for testing. To test for any kind
#' of differential usage, leave empty. To test for differential UTR usage, use
#' `excludeTypes=c("CDS","non-coding")`.
#'
#' @return The `se` object with additional rowData columns contain bin (i.e. exon) -level
#' statistics, and a metadata slot containing gene level p-values.
#'
#' @importFrom edgeR DGEList calcNormFactors glmQLFit glmFit diffSpliceDGE filterByExpr
#' @aliases DEUwrappers
#' @export
#' @rdname DEUwrappers
diffSpliceDGE.wrapper <- function(se, design, coef=NULL, QLF=TRUE, robust=TRUE,
                                  countFilter=TRUE, excludeTypes=NULL){
  if(is(design, "formula"))
    design <- model.matrix(design, data=as.data.frame(colData(se)))
  if(is.null(coef)){
    coef <- ifelse(is.null(colnames(design)), ncol(design), rev(colnames(design))[1])
    message("Testing coefficient ", coef)
  }
  if(countFilter) se <- se[filterByExpr(assays(se)$counts, design=design),]
  dds <- calcNormFactors(DGEList(assays(se)$counts))
  dds <- estimateDisp(dds,design)
  if(QLF){
    fit <- glmQLFit(dds, design, robust=robust)
  }else{
    fit <- glmFit(dds, design, robust=robust)
  }
  res <- diffSpliceDGE(fit, coef=coef, geneid=rowData(se)$gene, exonid=row.names(se))
  se <- se[names(res$exon.p.value),]
  rowData(se)$coefficient <- res$coefficients
  rowData(se)$bin.p.value <- res$exon.p.value

  ep <- res$exon.p.value
  if(!is.null(excludeTypes)) ep[rowData(se)$type %in% excludeTypes] <- 1
  metadata(se)$gene.q.values <- simes.aggregation(ep, res$genes[,2])
  se
}

#' @importFrom limma lmFit voom diffSplice
#' @importFrom edgeR DGEList calcNormFactors filterByExpr
#' @export
#' @rdname DEUwrappers
diffSplice.wrapper <- function(se, design, coef, robust=TRUE, filter=TRUE, improved=TRUE,
                               countFilter=TRUE, excludeTypes=NULL){
  if(is(design, "formula"))
    design <- model.matrix(design, data=as.data.frame(colData(se)))
  if(is.null(coef)){
    coef <- ifelse(is.null(colnames(design)), ncol(design), rev(colnames(design))[1])
    message("Testing coefficient ", coef)
  }

  if(countFilter) se <- se[filterByExpr(assays(se)$counts, design=design),]
  dds <- calcNormFactors(DGEList(assays(se)$counts))
  dds <- voom(dds, design)
  dds <- lmFit(dds, design)

  if(improved){
    res <- diffSplice2(dds, geneid=rowData(se)$gene, exonid=row.names(se), robust=robust)
  } else{
    res <- diffSplice(dds, geneid=rowData(se)$gene, exonid=row.names(se), robust=robust)
  }

  se <- se[row.names(res$p.value),]
  rowData(se) <- cbind(rowData(se), res$coefficients, bin.p.value=res$p.value)

  ep <- res$p.value
  if(!is.null(excludeTypes)) ep[rowData(se)$type %in% excludeTypes] <- 1
  metadata(se)$gene.q.values <- simes.aggregation(ep, res$genes[,2])
  se
}


#' @importFrom DEXSeq DEXSeqDataSet estimateSizeFactors estimateDispersions
#' @importFrom DEXSeq estimateExonFoldChanges DEXSeqResults perGeneQValue
#' @export
#' @rdname DEUwrappers
DEXSeq.wrapper <- function(se, design=~sample+exon+condition:exon,
                           reducedModel=~sample+exon, excludeTypes=NULL, ...){
  e <- floor(as.matrix(assays(se)$counts))
  e <- matrix(as.integer(e), nrow=nrow(e))
  dds <- DEXSeqDataSet(e, sampleData=as.data.frame(colData(se)),
                         design=design, featureRanges=rowRanges(se),
                         featureID=row.names(se), groupID=rowData(se)$gene)
  dds <- estimateSizeFactors( dds )
  dds <- estimateDispersions( dds )
  dds <- testForDEU( dds, reducedModel=reducedModel,...)
  vars <- setdiff(labels(terms(design)), labels(terms(reducedModel)))
  vars <- gsub(":exon|exon:", "", grep(":exon|exon:", vars, value=TRUE))
  dds <- estimateExonFoldChanges( dds, fitExpToVar=vars, ... )
  res <- DEXSeqResults( dds )

  exonBaseMean <- res$exonBaseMean
  exon.p.value <- res$pvalue
  log2fc <- res[[grep("log2fold",names(res))]]
  rowData(se) <- cbind(rowData(se),exonBaseMean,log2fc,bin.p.value=exon.p.value)

  ep <- exon.p.value
  if(!is.null(excludeTypes)) res$pvalue[rowData(se)$type %in% excludeTypes] <- 1
  metadata(se)$gene.q.values <- perGeneQValue(res)

  se
}
