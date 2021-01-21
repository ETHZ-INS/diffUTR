#' geneBinHeatmap
#'
#' @param se A bin-wise SummarizedExperiment as produced by \code{\link{countFeatures}}
#' @param gene The gene of interest
#' @param what Type of values (i.e. assay) to plot
#' @param anno_rows Row annotation columns (i.e. columns of `rowData(se)`) to plot
#' @param ... Passed to `sechm` (see \code{\link[SEtools]{SE-heatmap}}).
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap}}
#' @export
#' @import SummarizedExperiment
#' @importFrom SEtools sechm
#'
#' @examples
geneBinHeatmap <- function(se, gene, what=c("logNormDensity", "logCPM", "scaledLogCPM"),
                           anno_rows=c("logWidth","meanLogCPM","meanLogDensity","log10PValue"), ...){
  what <- match.arg(what, unique(c(what, assayNames(se))))
  if("log10PValue" %in% anno_rows && !("log10PValue" %in% colnames(rowData(se))))
    rowData(se)[["log10PValue"]] <- -log10(rowData(se)$bin.p.value)
  sechm(se, row.names(de), do.scale=what=="scaledLogCPM", show_rownames=FALSE,
        sortRowsOn=NULL, cluster_rows=FALSE, row_title=paste(gene, "bins"),
        anno_rows=anno_rows, assayName=ifelse(what=="logNormDensity",what,"logcpm"), ...)
}

#' deuBinPlot
#'
#' @param se A bin-wise SummarizedExperiment as produced by \code{\link{countFeatures}}
#' and including bin-level tests (i.e. having been passed through one of the DEU wrappers
#' such as \code{\link{diffSplice.wrapper}} or \code{\link{DEXSeq.wrapper}})
#' @param gene The gene of interest
#' @param type Either 'summary' (plot DEU summary), 'sample' (plot sample-wise data), or
#' 'condition' (plot data aggregate by condition)
#' @param intronSize Intron plot size. If <=1, total intron size will be this fraction of
#' the total exon size. If >1, each intron will have the given size.
#' @param exonSize Scaling for exon sizes, either 'sqrt', 'log', or 'linear'.
#' @param y Value to plot on the y-axis. If `type="summary"`, this should be a column of
#' `rowData(se)`, otherwise should be an assay name of `se`.
#' @param condition The colData column containing the samples' condition.
#' @param size rowData variable to use to determine the thickness of the bins.
#' @param lineSize Size of the line connecting the bins. Use `lineSize=0` to omit the line.
#' @param colour rowData variable to use to determine the colour of the bins (can also
#' be `condition` or `sample` when using the corresponding plot `type`)
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#'
#' @examples
deuBinPlot <- function(se, gene, type=c("summary","condition","sample"), intronSize=1,
                    exonSize=c("sqrt","linear","log"), y=NULL, condition=NULL,
                    size="meanLogDensity", lineSize=1, colour=NULL ){
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
  w <- rowData(se)$gene==gene
  if(length(w)==0) stop("Gene not found in the data!")
  if(type=="condition"){
    if(is.null(condition)) condition <- "condition"
    stopifnot(length(condition)==1 && condition %in% colnames(colData(se)))
  }
  if(is.null(y)){
    if(type=="summary"){
      for(f in c("coefficient","coefficients","log2fc","logFC"))
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

  de <- rowData(se)[w,]
  gr <- ranges(se)[w]
  de$log10PValue <- -log10(de$bin.p.value)
  de$order <- seq_len(nrow(de))
  de$width2 <- switch(exonSize, sqrt=sqrt(de$width), log=log(de$width), de$width)
  de$nextIsIntron <- c(start(gr)[-1],rev(end(gr))[1]) > end(gr)+1
  if(intronSize <= 1) intronSize <- intronSize*sum(de$width2)/sum(de$nextIsIntron)
  de$x_start <- c(0,cumsum(de$width2)+cumsum(de$nextIsIntron*intronSize))[-nrow(de)]
  de$x_end <- de$x_start + de$width2
  getd2 <- function(de) data.frame(x_start=de$x_end[-nrow(de)], x_end=de$x_start[-1],
                                   y_start=de[[y]][-nrow(de)], y_end=de[[y]][-1],
                                   type=de$nextIsIntron[-nrow(de)])
  p <- ggplot(as.data.frame(de))
  if(type=="summary"){
    if(lineSize>0){
      d2 <- getd2(de)
      p <- p + geom_segment(data=d2, aes(x=x_start, xend=x_end, y=y_start, yend=y_end,
                                         linetype=type), colour="grey", size=lineSize)
    }
  }else{
    y <- match.arg(y, assayNames(se))
    if(type=="condition"){
      e <- vapply(split(seq_len(ncol(se)), colData(se)[[condition]]),
                  FUN.VALUE=numeric(nrow(de)),
                  FUN=function(i) rowMeans(assays(se)[[y]][w,i,drop=FALSE]))
    }else{
      e <- assays(se)[[y]][w,]
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
                                       yend="y_end", linetype="type", colour=colour),
                            size=lineSize)
    }else{
      p <- p + geom_segment(data=d2,
                            aes_string(x="x_start", xend="x_end", y="y_start",
                                       yend="y_end", linetype="type", group=type),
                            colour="grey", size=lineSize)
    }
  }
  xlab <- switch(exonSize,
                 linear="Genomic location",
                 paste0(exonSize,"-scaled genomic location"))
  p + scale_linetype_manual(values=c("TRUE"="dotted", "FALSE"="solid"), guide=FALSE) +
    geom_segment(data=as.data.frame(de), aes_string(x="x_start", xend="x_end", y=y,
                                                    yend=y, size=size, colour=colour)) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle(gene) + xlab(xlab) + ylab(y)
}
