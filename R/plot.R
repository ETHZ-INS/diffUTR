#' geneBinHeatmap
#'
#' A wrapper around `ComplexHeatmap`.
#'
#' @param se A bin-wise SummarizedExperiment as produced by
#' \code{\link{countFeatures}}
#' @param gene The gene of interest
#' @param what Type of values (i.e. assay) to plot
#' @param anno_rows Row annotation columns (i.e. columns of `rowData(se)`) to
#' plot
#' @param anno_columns Column annotation columns (i.e. columns of
#' `colData(se)`) to plot
#' @param anno_colors Annotation colors, as a list named with the row/column
#' annotations, see `\code{\link[ComplexHeatmap]{SingleAnnotation}}` for
#' details. Ignored if `left_annotation` and/or `top_annotation` are given
#' directly.
#' @param removeAmbiguous Logical; whether to remove bins that are
#' gene-ambiguous (i.e. overlap multiple genes).
#' @param merge_legends Logical; whether to merge legends. This effectively
#' calls `draw(..., merge_legends=TRUE)` around the heatmap.
#' @param cluster_columns Logical; whether to cluster columns (passed to
#' \code{\link[ComplexHeatmap]{Heatmap}})
#' @param left_annotation Passed to \code{\link[ComplexHeatmap]{Heatmap}},
#' overrides `anno_rows`.
#' @param top_annotation Passed to \code{\link[ComplexHeatmap]{Heatmap}},
#' overrides `anno_columns`.
#' @param ... Passed to `ComplexHeatmap` (see
#' \code{\link[ComplexHeatmap]{Heatmap}})
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap}}
#' @export
#' @import SummarizedExperiment
#' @importFrom ComplexHeatmap draw Heatmap HeatmapAnnotation rowAnnotation
#'
#' @examples
#' data(example_bin_se)
#' se <- diffSpliceWrapper(example_bin_se, ~condition)
#' geneBinHeatmap(se, "Jund")
geneBinHeatmap <- function(se, gene,
                           what=c("logNormDensity", "logCPM", "scaledLogCPM"),
                           anno_rows=c("type","logWidth","meanLogDensity",
                                       "log10PValue","geneAmbiguous"),
                           anno_columns=c(), anno_colors=list(),
                           removeAmbiguous=FALSE, merge_legends=TRUE,
                           cluster_columns=FALSE,
                           left_annotation=NULL, top_annotation=NULL, ...){
  se <- .checkSE(se, checkNorm=TRUE)
  if(length(w <- .matchGene(se, gene))==0) stop("Gene not found!")
  se <- sort(se[w,])
  if(removeAmbiguous) se <- se[!rowData(se)$geneAmbiguous,]
  what <- match.arg(what)
  if("log10PValue" %in% anno_rows &&
     !("log10PValue" %in% colnames(rowData(se))))
    rowData(se)[["log10PValue"]] <- -log10(rowData(se)$bin.p.value)
  if("type" %in% anno_rows && is.null(anno_colors$type))
    anno_colors$type <- .typeColors()
  if(removeAmbiguous | !any(rowData(se)$geneAmbiguous))
    anno_rows <- intersect(setdiff(anno_rows, "geneAmbiguous"),
                           colnames(rowData(se)))
  if(is.null(left_annotation) && !is.null(anno_rows)){
    left_annotation <- rowAnnotation(
      df=as.data.frame(rowData(se)[,anno_rows,drop=FALSE]),
      col=anno_colors[intersect(names(anno_colors),anno_rows)])
  }
  if(is.null(top_annotation) && !is.null(anno_columns)){
    top_annotation <- HeatmapAnnotation(
      df=as.data.frame(colData(se)[,anno_columns,drop=FALSE]),
      col=anno_colors[intersect(names(anno_colors),anno_rows)])
  }
  x <- assays(se)[[ifelse(what=="logNormDensity",what,"logcpm")]]
  if(what=="scaledLogCPM") x <- t(scale(t(x)))
  h <- Heatmap(x, name=what, show_row_names=FALSE, cluster_rows=FALSE,
          cluster_columns=cluster_columns, row_title=paste(gene, "bins"),
          top_annotation=top_annotation, left_annotation=left_annotation, ...)
  if(merge_legends) return(draw(h, merge_legends=TRUE))
  h
}

#' deuBinPlot
#'
#' @param se A bin-wise SummarizedExperiment as produced by
#' \code{\link{countFeatures}} and including bin-level tests (i.e. having been
#' passed through one of the DEU wrappers such as
#' \code{\link{diffSpliceWrapper}} or \code{\link{DEXSeqWrapper}})
#' @param gene The gene of interest
#' @param type Either 'summary' (plot DEU summary), 'sample' (plot sample-wise
#' data), or 'condition' (plot data aggregate by condition)
#' @param intronSize Intron plot size. If <=3, intron size will be this
#' fraction of the mean exon size. If >3, each intron will have the given size.
#' @param exonSize Scaling for exon sizes, either 'sqrt', 'log', or 'linear'.
#' @param y Value to plot on the y-axis. If `type="summary"`, this should be a
#' column of `rowData(se)`, otherwise should be an assay name of `se`.
#' @param condition The colData column containing the samples' condition.
#' @param size rowData variable to use to determine the thickness of the bins.
#' @param lineSize Size of the line connecting the bins. Use `lineSize=0` to
#' omit the line.
#' @param colour rowData variable to use to determine the colour of the bins.
#' If `type="condition"`, can also be "condition"; if `type="sample"` can be
#' any colData column.
#' @param alpha Alpha level, passed to ggplot.
#' @param removeAmbiguous Logical; whether to remove bins that are
#' gene-ambiguous (i.e. overlap multiple genes).
#'
#' @return A ggplot object
#' @importFrom matrixStats colMedians
#' @importFrom dplyr bind_rows
#' @export
#' @import ggplot2
#'
#' @examples
#' data(example_bin_se)
#' se <- diffSpliceWrapper(example_bin_se, ~condition)
#' deuBinPlot(se, "Jund")
deuBinPlot <- function(se, gene, type=c("summary","condition","sample"),
                       intronSize=2, exonSize=c("sqrt","linear","log"), y=NULL,
                       condition=NULL, size="type", lineSize=1,
                       colour=NULL, alpha=NULL, removeAmbiguous=TRUE ){
  se <- .checkSE(se)
  type <- match.arg(type)
  exonSize <- match.arg(exonSize)
  stopifnot(length(gene)==1)
  if(is.null(colour)){
    if(type!="summary" && !is.null(condition)){
      colour <- "condition"
    }else{
      colour <- "log10PValue"
    }
  }
  if(is.null(alpha)) alpha <- ifelse(type=="sample", 0.6, 1)
  w <- .matchGene(se, gene)
  if(length(w)==0) stop("Gene not found in the data!")
  if(type=="condition"){
    if(is.null(condition) && "condition" %in% colnames(colData(se))){
      condition <- "condition"
    }else{
      condition <- "group"
    }
    stopifnot(length(condition)==1 && condition %in% colnames(colData(se)))
  }
  if(is.null(y)){
    if(type=="summary"){
      ff <- colnames(rowData(se))[grep("bin.p.value",colnames(rowData(se)))-1]
      for(f in c("coefficient","coefficients","log2fc","logFC",ff))
        if(is.null(y) && f %in% colnames(rowData(se))) y <- f
      if(is.null(y)) stop("Please specify `y` among the rowData columns.")
    }else{
      if("logNormDensity" %in% assayNames(se)){
        y <- "logNormDensity"
      }else{
        y <- assayNames(se)[1]
        message("Using assay '",y, "' (use `y` to specify another)")
      }
    }
  }

  se <- sort(se[w,])
  if(removeAmbiguous) se <- se[!rowData(se)$geneAmbiguous,]
  de <- rowData(se)
  if(y=="geneNormDensity" && !("geneNormDensity" %in% assayNames(se))){
    si <- vapply(split(seq_len(ncol(se)), se[[condition]]),
                 FUN.VALUE=numeric(1), FUN=function(i){
                   mean(matrixStats::colMedians(assays(se)$logNormDensity[,i]))
                 })
    assays(se)$geneNormDensity <-
      t(t(assays(se)$logNormDensity)-si[se[[condition]]])
  }

  if(size=="type"){
    levels(de$type)[grepl("CDS/.*UTR", levels(de$type))] <- "CDS/UTR"
    levels(de$type)[grep("CDS/.*UTR", levels(de$type))] <- "CDS/UTR"
    levels(de$type)[grep("CDS|non-coding", levels(de$type), invert=TRUE)] <-
      "UTR"
    de$type <- factor(de$type, c("non-coding","UTR","CDS/UTR","CDS"))
  }

  gr <- ranges(se)
  de$log10PValue <- -log10(de$bin.p.value)
  de$order <- seq_len(nrow(de))
  de$width <- width(gr)
  de$width2 <- switch(exonSize, sqrt=sqrt(de$width), log=log(de$width),
                      de$width)
  de$nextIsIntron <- c(start(gr)[-1],rev(end(gr))[1]) > end(gr)+1
  if(intronSize <= 3) intronSize <- intronSize*mean(de$width2)
  de$x_start <- c(0,cumsum(de$width2)+
                    cumsum(de$nextIsIntron*intronSize))[-nrow(de)]
  de$x_end <- de$x_start + de$width2
  getd2 <- function(de)
    data.frame( x_start=de$x_end[-nrow(de)], x_end=de$x_start[-1],
                y_start=de[[y]][-nrow(de)], y_end=de[[y]][-1],
                type=de$nextIsIntron[-nrow(de)])
  p <- ggplot(as.data.frame(de))
  if(type=="summary"){
    if(lineSize>0){
      d2 <- getd2(de)
      p <- p + geom_segment(data=d2,
                            aes(x=x_start, xend=x_end, y=y_start, yend=y_end,
                                linetype=type),
                            colour="grey", size=lineSize, alpha=alpha)
    }
  }else{
    y <- match.arg(y, assayNames(se))
    if(type=="condition"){
      e <- vapply(split(seq_len(ncol(se)), colData(se)[[condition]]),
                  FUN.VALUE=numeric(nrow(de)),
                  FUN=function(i) rowMeans(assays(se)[[y]][,i,drop=FALSE]))
    }else{
      e <- assays(se)[[y]][,]
    }
    de <- lapply(setNames(colnames(e),colnames(e)), FUN=function(x){
      de2 <- as.data.frame(de)
      de2[[y]] <- as.numeric(e[,x])
      de2
    })
    if(lineSize>0) d2 <- dplyr::bind_rows(lapply(de, FUN=getd2), .id=type)
    de <- dplyr::bind_rows(de, .id=type)
    if(type=="sample"){
      de <- cbind(de, colData(se)[de$sample,,drop=FALSE])
      if(lineSize>0) d2 <- cbind(d2, colData(se)[d2$sample,,drop=FALSE])
    }
    if(type=="sample" && lineSize>0 && colour %in% colnames(d2)){
      p <- p + geom_segment(data=d2,
                          aes_string(x="x_start", xend="x_end", y="y_start",
                                     yend="y_end", linetype="type",
                                     colour=colour),
                          size=lineSize, alpha=min(0.6,alpha))
    }else{
      p <- p + geom_segment(data=d2,
                        aes_string(x="x_start", xend="x_end", y="y_start",
                                   yend="y_end", linetype="type", group=type),
                        colour="grey", size=lineSize, alpha=min(0.6,alpha))
    }
  }
  if(colour=="type" && type!="sample")
    p <- p + scale_colour_manual(values=.typeColors())
  if(size=="type")
    p <- p + scale_size_manual(values=c("non-coding"=2, "UTR"=3, "CDS/UTR"=4,
                                        CDS=6))
  xlab <- switch(exonSize,
                 linear="Genomic location",
                 paste0(exonSize,"-scaled genomic location"))
  p + scale_linetype_manual(values=c("TRUE"="dotted", "FALSE"="solid"),
                            guide=FALSE) +
    geom_segment(data=as.data.frame(de), alpha=alpha,
                 aes_string(x="x_start", xend="x_end", y=y, yend=y, size=size,
                            colour=colour)) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle(gene) + xlab(xlab) +
    ylab(ifelse(y=="geneNormDensity","Gene-normalized logDensity",y))
}


#' plotTopGenes
#'
#' @param se A bin-wise SummarizedExperiment as produced by
#' \code{\link{countFeatures}} and including bin-level tests (i.e. having been
#' passed through one of the DEU wrappers such as
#' \code{\link{diffSpliceWrapper}} or \code{\link{DEXSeqWrapper}})
#' @param n The maximum number of genes for which to plot labels
#' @param FDR The FDR threshold above which to plot labels
#' @param diffUTR Logical; if FALSE, uses absolute coefficients (appropriate
#' for normal differential exon usage); if TRUE, uses non-absolute (ie changes
#' should be in the same direction across significant bins) and width-weighted
#' scores (i.e. larger bins have more weight) -- this is relevant only when
#' testing UTR usage.
#' @param alpha Points transparency
#'
#' @return A ggplot
#' @export
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' data(example_bin_se)
#' se <- diffSpliceWrapper(example_bin_se, ~condition)
#' plotTopGenes(se)
plotTopGenes <- function(se, n=25, FDR=0.05, diffUTR=FALSE, alpha=NULL){
  if(is(se, "SummarizedExperiment")){
    se <- .checkSE(se, requireStats=TRUE)
    se <- metadata(se)$geneLevel
  }
  stopifnot(!is.null(se) && (is.data.frame(se) || is(se,"DFrame")))
  se <- as.data.frame(se)
  if(length(w <- which(se$q.value<=0))>0){
    minnzq <- -log10(min(se$q.value[-w],na.rm=TRUE))
    score <- abs(se[[ifelse(diffUTR,"sizeScore","w.abs.coef")]])
    o <- order(order(score[w]))
    se$q.value[w] <- 10^-(minnzq+o*minnzq/(4*length(w)))
  }
  if(diffUTR){
    if(is.null(alpha)) alpha <- 1
    p <- ggplot(se, aes(sizeScore, -log10(q.value))) +
      geom_point(aes(colour=w.coef, size=geneMeanDensity), alpha=alpha) +
      scale_colour_gradient2() + xlab("Weighted size score")
    de <- se[se$q.value<FDR,,drop=FALSE]
    de$tmp <- abs(-log10(de$q.value)*de$sizeScore)
  }else{
    se$tmp <- -log10(se$q.value)*se$w.abs.coef
    if(is.null(alpha)){
      p <- ggplot(se, aes(w.abs.coef, -log10(q.value))) +
        geom_point(aes(colour=density.ratio, size=geneMeanDensity,
                       alpha=abs(tmp))) + scale_alpha(guide=FALSE)

    }else{
      p <- ggplot(se, aes(w.abs.coef, -log10(q.value))) +
        geom_point(aes(colour=density.ratio, size=geneMeanDensity),
                   alpha=alpha)
    }
    p <- p + xlab("Weighted absolute coefficient")
    de <- se[se$q.value<FDR,,drop=FALSE]
  }
  if(n>0){
    de <- de[head(order(de$tmp, decreasing=TRUE),n),,drop=FALSE]
    if(is.null(de$name)) de$name <- row.names(de)
    p <- p + geom_text_repel(data=de, aes(label=name))
  }
  p
}
