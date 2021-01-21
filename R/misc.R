#' addNormalizedAssays
#'
#' @param se A bin-wise SummarizedExperiment as produced by \code{\link{countFeatures}}
#'
#' @return The `se` object with populated `logcpm` and `logNormDensity` assays.
#' @export
#'
#' @importFrom edgeR DGEList cpm calcNormFactors
#' @examples
addNormalizedAssays <- function(se){
  assays(se)$logcpm <- log1p(cpm(calcNormFactors(DGEList(assay(se)))))
  assays(se)$logNormDensity <- log1p(exp(assays(se)$logcpm)/width(se))
  se
}

#' simes.aggregation
#'
#' Simes p-value correction and aggregation, adapted from \code{link[limma]{diffSplice}}
#'
#' @param p.value A vector of p-values
#' @param geneid A vector of group labels such as gene identifiers
#'
#' @return A named vector of aggregated p-values
#' @export
#'
#' @examples
simes.aggregation <- function(p.value, geneid){
  stopifnot(is.numeric(p.value))
  stopifnot(length(p.value)==length(geneid))
  o <- order(geneid)
  geneid <- geneid[o]
  p.value <- p.value[o]

  ngenes <- length(unique(geneid))

  gene.nexons <- rowsum(rep(1,length(p.value)), geneid, reorder=FALSE)
  g <- rep(1:ngenes, times=gene.nexons)

  gene.lastexon <- cumsum(gene.nexons)
  gene.firstexon <- gene.lastexon-gene.nexons+1

  penalty <- rep_len(1L,length(g))
  names(p.value)<-geneid
  penalty[gene.lastexon] <- 1L-gene.nexons
  penalty <- cumsum(penalty)[-gene.lastexon]
  penalty <- penalty / rep(gene.nexons-1L,gene.nexons-1L)
  g2 <- g[-gene.lastexon]

  gene.simes.p.value <- rep(1, length(gene.nexons))

  o <- order(g,p.value)
  p.adj <- pmin(p.value[o][-gene.lastexon] / penalty, 1)
  o <- order(g2,p.adj)
  gene.simes.p.value <- p.adj[o][gene.firstexon-0L:(ngenes-1L)]

  gene.simes.p.value
}
