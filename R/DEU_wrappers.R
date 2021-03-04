#' DEUwrappers
#'
#' Wrappers around commonly-used DEU methods
#' (\code{\link[edgeR]{diffSpliceDGE}}, \code{\link[DEXSeq]{DEXSeq}} and an
#' improved version of \code{\link[limma]{diffSplice}}
#'
#' @param se A bin-wise SummarizedExperiment as produced by
#' \code{\link{countFeatures}}
#' @param design A formula (using columns of `colData(se)`) or (for
#' `diffSpliceWrapper` or `diffSpliceDGEWrapper` only) a model.matrix.
#' @param reducedModel A reduced formula (applicable only to `DEXSeqWrapper`).
#' @param coef The coefficient to be tested (ignored for `DEXSeqWrapper`).
#' @param QLF Logical; whether to use edgeR's quasi-likelihood negative binomial
#' (applicable only to `diffSpliceDGEWrapper`).
#' @param robust Logical; whether to use robust fitting for the dispersion trend
#' (ignored for `DEXSeqWrapper`).
#' @param countFilter Logical; whether to filter out low-count bins (ignored for
#' `DEXSeqWrapper`).
#' @param excludeTypes A vector of bin types to ignore for testing. To test for
#'  any kind of differential usage, leave empty. To test for differential UTR
#'  usage, use `excludeTypes=c("CDS","non-coding")` (or see
#'  \code{\link{geneLevelStats}} for more options).
#'
#' @return The `se` object with additional rowData columns contain bin (i.e.
#' exon) -level statistics, and a metadata slot containing gene level p-values.
#'
#' @importFrom edgeR DGEList calcNormFactors glmQLFit glmFit diffSpliceDGE
#' @importFrom edgeR filterByExpr estimateDisp
#' @importFrom stats p.adjust setNames
#' @importFrom matrixStats rowMins
#' @aliases DEUwrappers
#' @export
#' @rdname DEUwrappers
#' @examples
#' library(SummarizedExperiment)
#' data(example_bin_se)
#' se <- diffSpliceWrapper(example_bin_se, ~condition)
#' head(rowData(se))
diffSpliceDGEWrapper <- function(se, design, coef=NULL, QLF=TRUE, robust=TRUE,
                                  countFilter=TRUE, excludeTypes=NULL){
  se <- .checkSE(se)
  if(is(design, "formula"))
    design <- model.matrix(design, data=as.data.frame(colData(se)))
  if(is.null(coef)){
    coef <- ifelse(is.null(colnames(design)), ncol(design),
                   rev(colnames(design))[1])
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
  res <- diffSpliceDGE(fit, coef=coef, geneid=rowData(se)$gene,
                       exonid=row.names(se))
  se <- se[names(res$exon.p.value),]
  co <- res$coefficients
  if(!is.null(dim(co))) co <- co[,coef]
  if(length(coef)>1 && !is.null(dim(res$exon.p.value))){
    rowData(se) <- cbind(rowData(se), co)
    pv <- res$exon.p.value[,coef]
    colnames(pv) <- paste0(colnames(pv),".p.value")
    rowData(se)$coefficient <- mapply(i=seq_len(nrow(co)),
                                      j=apply(pv,1,which.min),
                                      FUN=function(i,j) co[i,j])
    rowData(se) <- cbind(rowData(se), pv)
    rowData(se)$bin.p.value <- ep <- rowMins(pv)
  }else{
    rowData(se)$coefficient <- co
    rowData(se)$bin.p.value <- ep <- res$exon.p.value
  }

  if(!is.null(excludeTypes)) ep[rowData(se)$type %in% excludeTypes] <- 1
  rowData(se)$bin.FDR <- p.adjust(ep)

  d <- DataFrame(bin.pval=ep,coef=rowData(se)$coefficient,gene=rowData(se)$gene,
                 width=width(se), meanLogDensity=rowData(se)$meanLogDensity)
  if("gene_name" %in% colnames(rowData(se)))
    d$gene_name <- rowData(se)$gene_name
  metadata(se)$geneLevel <- .geneLevelStats(d=d)
  se
}

#' @param improved Logical; whether to use \code{\link{diffSplice2}} instead of
#'  the original \code{\link[limma]{diffSplice}} (default TRUE).´
#' @importFrom stats model.matrix
#' @importFrom limma lmFit voom diffSplice
#' @importFrom edgeR DGEList calcNormFactors filterByExpr
#' @export
#' @rdname DEUwrappers
diffSpliceWrapper <- function(se, design, coef=NULL, robust=TRUE,
                               improved=TRUE, countFilter=TRUE,
                               excludeTypes=NULL){
  se <- .checkSE(se)
  if(is(design, "formula"))
    design <- model.matrix(design, data=as.data.frame(colData(se)))
  if(is.null(coef)){
    coef <- ifelse(is.null(colnames(design)), ncol(design),
                   rev(colnames(design))[1])
    message("Testing coefficient ", coef)
  }

  if(countFilter) se <- se[filterByExpr(assays(se)$counts, design=design),]
  dds <- calcNormFactors(DGEList(assays(se)$counts))
  dds <- voom(dds, design)
  dds <- lmFit(dds, design)

  if(improved){
    res <- diffSplice2(dds, geneid=rowData(se)$gene, exonid=row.names(se),
                       robust=robust)
  } else{
    res <- diffSplice(dds, geneid=rowData(se)$gene, exonid=row.names(se),
                      robust=robust)
  }

  se <- se[row.names(res$p.value),]

  tmp <- rowData(se)[,setdiff(colnames(rowData(se)),
                    c(colnames(res$coefficients),"bin.p.value","coefficient"))]
  if(length(coef)>1){
    co <- res$coefficients[,coef]
    pv <- res$p.value[,coef]
    colnames(pv) <- paste0(colnames(pv),"p.value")
    tmp$coefficient <- mapply(i=seq_len(nrow(co)), j=apply(pv,1,which.min),
                                      FUN=function(i,j) co[i,j])
    tmp <- cbind(tmp, co, pv)
    ep <- rowMins(pv)
  }else{
    tmp$coefficient <- res$coefficients[,coef]
    ep <- res$p.value[,coef]
  }
  tmp$bin.p.value <- ep
  rowData(se) <- tmp
  rowData(se)$bin.FDR <- p.adjust(ep)

  if(!is.null(excludeTypes)) ep[rowData(se)$type %in% excludeTypes] <- 1
  d <- DataFrame(bin.pval=ep, coef=rowData(se)$coefficient,
                 gene=rowData(se)$gene, width=width(se),
                 meanLogDensity=rowData(se)$meanLogDensity)
  if("gene_name" %in% colnames(rowData(se)))
    d$gene_name <- rowData(se)$gene_name
  metadata(se)$geneLevel <- .geneLevelStats(d=d)
  se
}

#' @param ... Further arguments (passed to `testForDEU` and
#' `estimateExonFoldChanges`)
#' @importFrom DEXSeq DEXSeqDataSet estimateSizeFactors estimateDispersions
#' @importFrom DEXSeq estimateExonFoldChanges DEXSeqResults perGeneQValue
#' @importFrom DEXSeq testForDEU
#' @export
#' @rdname DEUwrappers
DEXSeqWrapper <- function(se, design=~sample+exon+condition:exon,
                           reducedModel=~sample+exon, excludeTypes=NULL, ...){
  if(!("exon" %in% labels(terms(design))))
    stop("For DEXSeq, the formula should include the extra 'sample' and 'exon'",
         " terms.\nFor instance, if you wanted to test for an effect of the ",
         "variable 'condition' on bin usage, you would use:\n",
         "~sample+exon+condition:exon")
  se <- .checkSE(se)
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

  rowData(se)$log2fc <- res[[grep("log2fold",names(res))]]
  rowData(se)$exonBaseMean <- res$exonBaseMean
  rowData(se)$bin.p.value <- res$pvalue

  if(!is.null(excludeTypes)) res$pvalue[rowData(se)$type %in% excludeTypes] <- 1
  rowData(se)$bin.FDR <- p.adjust(res$pvalue)

  d <- DataFrame(bin.pval=res$pvalue, coef=rowData(se)$log2fc,
                 gene=rowData(se)$gene, width=width(se),
                 meanLogDensity=rowData(se)$meanLogDensity)
  if("gene_name" %in% colnames(rowData(se)))
    d$gene_name <- rowData(se)$gene_name
  message("Generating gene-level stats...")
  metadata(se)$geneLevel <- .geneLevelStats(d=d, gene.qval=perGeneQValue(res))

  se
}
