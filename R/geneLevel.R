#' simesAggregation
#'
#' Simes p-value correction and aggregation, adapted from
#' \code{link[limma]{diffSplice}}
#'
#' @param p.value A vector of p-values
#' @param geneid A vector of group labels such as gene identifiers
#'
#' @return A named vector of aggregated p-values
#' @export
#'
#' @examples
#' p <- runif(50)
#' genes <- sample(LETTERS,50,replace=TRUE)
#' simesAggregation(p, genes)
simesAggregation <- function(p.value, geneid){
  stopifnot(is.numeric(p.value))
  stopifnot(length(p.value)==length(geneid))
  o <- order(geneid)
  geneid <- geneid[o]
  p.value <- p.value[o]

  ngenes <- length(unique(geneid))

  gene.nexons <- rowsum(rep(1,length(p.value)), geneid, reorder=FALSE)
  g <- rep(seq_len(ngenes), times=gene.nexons)

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

#' geneLevelStats
#'
#' Aggregates bin-level statistics to the gene-level
#'
#' @param se A `RangedSummarizedExperiment` containing the results of one of
#' the DEU wrappers.
#' @param coef The coefficients tested (if the model included more than one
#' term).
#' @param excludeTypes Vector of bin types to exclude.
#' @param includeTypes Vector of bin types to include (overrides
#' `excludeTypes`)
#' @param returnSE Logical; whether to return the updated `se` object
#' (default), or the gene-level table.
#'
#' @return If `returnSE=TRUE` (default), returns the `se` object with an
#' updated `metadata(se)$geneLevel` slot, otherwise returns the gene-level
#' data.frame.
#' @export
#'
#' @import S4Vectors
#' @importFrom stats weighted.mean
#' @examples
#' library(SummarizedExperiment)
#' data(example_bin_se)
#' se <- diffSpliceWrapper(example_bin_se, ~condition)
#' se <- geneLevelStats(se, includeTypes="3UTR")
#' head(metadata(se)$geneLevel)
geneLevelStats <- function(se, coef=NULL, excludeTypes=NULL, includeTypes=NULL,
                           returnSE=TRUE){
  rd <- rowData(se)
  if(!is.null(includeTypes)){
    excludeTypes <- setdiff(levels(rd$type), includeTypes)
  }
  if(!is.null(excludeTypes)){
    rd$bin.p.value[which(rd$type %in% excludeTypes)] <- 1
  }
  if(is.null(coef)){
    tmp <- setdiff(colnames(rd),"exonBaseMean")
    coef <- colnames(rd)[grep("bin.p.value",tmp)-1]
  }

  metadata(se)$geneLevel <- .geneLevelStats(DataFrame(
    bin.pval=rd$bin.p.value, coef=rd[[coef]], width=width(se),
    gene=rd$gene, gene_name=rd$gene_name, meanLogDensity=rd$meanLogDensity))

  if(!returnSE) return(metadata(se)$geneLevel)
  se
}

#' @importFrom methods is as
#' @importFrom BiocParallel bpnworkers
#' @importFrom IRanges LogicalList NumericList IntegerList
.geneLevelStats <- function(d, gene.qval=NULL){
  stopifnot(c("bin.pval","coef","gene","width","meanLogDensity") %in%
              colnames(d))
  d$bin.pval[is.na(d$bin.pval)] <- 1
  if(is.null(gene.qval)) gene.qval <- simesAggregation(d$bin.pval, d$gene)
  d$gene <- droplevels(as.factor(d$gene))
  d$coef[is.na(d$coef)] <- 0
  g0 <- names(which(!any(LogicalList(split(d$bin.pval<1, d$gene)))))
  d$weight <- -log10(d$bin.pval)
  d$weight[is.na(d$weight)] <- 0
  d$weight[d$gene %in% g0] <- 1
  w <- NumericList(split(d$weight, d$gene))
  w <- w/sum(w)
  co <- NumericList(split(d$coef, d$gene))
  wi <- IntegerList(split(d$width, d$gene))
  de <- NumericList(split(d$meanLogDensity, d$gene))
  d2 <- data.frame( row.names = names(co),
                    nb.bins = lengths(co),
                    w.coef=sum(w*co),
                    w.abs.coef=sum(w*abs(co)),
                    w.width=sum(w*wi),
                    w.density=log2(sum(w*de)),
                    sizeScore=sum(wi/sum(wi)*w*abs(co)),
                    abs.sizeScore=sum(wi/sum(wi)*w*co),
                    geneMeanDensity=mean(de,trim=0.05) )
  d2$density.ratio <- sum(w*(de-d2$geneMeanDensity))

  for(f in c(1:2,5:8)) d2[,f] <- round(d2[,f],3)
  d2$w.width <- as.integer(round(d2$w.width))
  d2$q.value <- gene.qval[row.names(d2)]
  if("gene_name" %in% colnames(d)){
    d <- d[order(d$gene, lengths(d$gene_name)),c("gene","gene_name")]
    d <- d[!duplicated(d$gene),]
    row.names(d) <- d$gene
    if(is(d$gene_name,"FactorList"))
      d$gene_name <- as(d$gene_name, "CharacterList")
    if(is(d$gene_name,"CharacteList")){
      d$name <- unstrsplit(d$gene_name, sep = "/")
    }else{
      d$name <- d$gene_name
    }
    d$name <- unstrsplit(d$gene_name, sep = "/")
    d2 <- merge(d[,"name",drop=FALSE], d2, by="row.names")
    row.names(d2) <- d2$Row.names
    d2 <- d2[,-1]
  }
  d2[order(d2$q.value),]
}
