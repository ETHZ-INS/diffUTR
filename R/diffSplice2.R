#' diffSplice2
#'
#' This is a small improvement to the \code{\link[limma]{diffSplice}} function
#' written by Gordon Smyth and Charity Law.
#'
#' @param fit an \code{\link[limma]{MArrayLM-class}} fitted model object
#' produced by \code{\link[limma]{lmFit}} or `contrasts.fit`, with rows
#' corresponding to exons.
#' @param geneid gene identifiers (as in \code{\link[limma]{diffSplice}})
#' @param exonid exon identifiers (as in \code{\link[limma]{diffSplice}})
#' @param robust logical, should the estimation of the empirical Bayes prior
#' parameters be robustified against outlier sample variances?
#' @param verbose logical, if TRUE will output some diagnostic information
#'
#' @return An \code{\link[limma]{MArrayLM-class}} objectcontaining both exon
#' level and gene level tests. Results are sorted by geneid and by exonid
#' within gene.
#' @export
#' @import limma
#' @importFrom methods new
#' @importFrom stats ave pf pt terms
#' @examples
#' library(SummarizedExperiment)
#' library(limma)
#' data(example_bin_se)
#' se <- example_bin_se
#' design <- model.matrix(~condition, data=as.data.frame(colData(se)))
#' dds <- calcNormFactors(DGEList(assays(se)$counts))
#' dds <- voom(dds, design)
#' dds <- lmFit(dds, design)
#' res <- diffSplice2(dds, geneid=rowData(se)$gene, exonid=row.names(se))
#' topSplice(res)
diffSplice2 <- function(fit, geneid, exonid=NULL, robust=FALSE, verbose=TRUE){
  # Exon Level squeeze, exon.s2.post t-test and weighting,
  # sum of exon.df as gene.df after squeeze
  exon.genes <- fit$genes
  if(is.null(exon.genes)) exon.genes <- data.frame(ExonID=seq_len(nrow(fit)))

  # Get ID columns for genes and exons
  if(length(geneid)==1) {
    genecolname <- as.character(geneid)
    geneid <- exon.genes[[genecolname]]
  } else {
    exon.genes$GeneID <- geneid
    genecolname <- "GeneID"
  }
  if(is.null(exonid)) {
    exoncolname <- NULL
  } else {
    if(length(exonid)==1) {
      exoncolname <- as.character(exonid)
      exonid <- exon.genes[[exoncolname]]
    } else {
      exon.genes$ExonID <- exonid
      exoncolname <- "ExonID"
    }
  }

  # Treat NA geneids as genes with one exon
  if(anyNA(geneid)) {
    isna <- which(is.na(geneid))
    geneid[isna] <- paste0("NA",seq_along(isna))
  }

  # Sort by geneid
  if(is.null(exonid))
    o <- order(geneid)
  else
    o <- order(geneid,exonid)
  geneid <- geneid[o]
  exon.genes <- exon.genes[o,,drop=FALSE]
  exon.coefficients <- fit$coefficients[o,,drop=FALSE]
  exon.stdev.unscaled <- fit$stdev.unscaled[o,,drop=FALSE]
  exon.df.residual <- fit$df.residual[o]
  exon.s2 <- fit$sigma[o]^2

  # Count exons by gene and get genewise variances
  exon.stat <- cbind(1,exon.df.residual,exon.s2)
  gene.sum <- rowsum(exon.stat,geneid,reorder=FALSE)
  gene.nexons <- gene.sum[,1]
  gene.df.residual <- gene.sum[,2]
  gene.s2 <- gene.sum[,3] / gene.sum[,1]
  if(verbose) {
    cat("Total number of exons: ", length(geneid), "\n")
    cat("Total number of genes: ", length(gene.nexons), "\n")
    cat("Number of genes with 1 exon: ", sum(gene.nexons==1), "\n")
    cat("Mean number of exons in a gene: ", round(mean(gene.nexons),0), "\n")
    cat("Max number of exons in a gene: ", max(gene.nexons), "\n")
  }

  # Posterior genewise variances
  squeeze <- squeezeVar(var=exon.s2, df=exon.df.residual, robust=robust)


  # Remove genes with only 1 exon
  gene.keep <- gene.nexons>1
  ngenes <- sum(gene.keep)
  if(ngenes==0) stop("No genes with more than one exon")

  exon.keep <- rep(gene.keep,gene.nexons)
  geneid <- geneid[exon.keep]
  exon.genes <- exon.genes[exon.keep,,drop=FALSE]
  exon.coefficients <- exon.coefficients[exon.keep,,drop=FALSE]
  exon.stdev.unscaled <- exon.stdev.unscaled[exon.keep,,drop=FALSE]
  exon.df.residual <- exon.df.residual[exon.keep]

  gene.nexons <- gene.nexons[gene.keep]
  gene.df.test <- gene.nexons-1
  gene.df.residual <- gene.df.residual[gene.keep]
  if(robust) squeeze$df.prior <- squeeze$df.prior[gene.keep]
  exon.df.total <- exon.df.residual
  exon.df.total <- pmin(exon.df.total,sum(exon.df.residual))
  gene.df.total <- rowsum(exon.df.total,geneid,reorder=FALSE)
  exon.s2.post <- squeeze$var.post[exon.keep]

  # Genewise betas
  u2 <- 1/(exon.stdev.unscaled^2*exon.s2.post)
  u2.rowsum <- rowsum(u2,geneid,reorder=FALSE)
  gene.betabar <- rowsum(exon.coefficients*u2,geneid,reorder=FALSE) / u2.rowsum

  # T-statistics for exon-level tests
  g <- rep(seq_len(ngenes), times=gene.nexons)
  exon.coefficients <- exon.coefficients-gene.betabar[g,,drop=FALSE]
  exon.1mleverage <- 1 - (u2 / u2.rowsum[g,,drop=FALSE])
  exon.coefficients <- exon.coefficients / exon.1mleverage
  exon.t <- exon.coefficients / exon.stdev.unscaled / sqrt(exon.s2.post)
  gene.F <- rowsum(exon.t^2,geneid,reorder=FALSE) / gene.df.test

  exon.p.value <- 2 * pt(abs(exon.t), df=gene.df.total[g], lower.tail=FALSE)
  gene.F.p.value <- pf(gene.F, df1=gene.df.test, df2=gene.df.total,
                       lower.tail=FALSE)

  # Exon level output
  out <- new("MArrayLM",list())
  out$genes <- exon.genes
  out$genecolname <- genecolname
  out$exoncolname <- exoncolname
  out$coefficients <- exon.coefficients
  out$t <- exon.t
  out$p.value <- exon.p.value

  # Gene level output
  out$gene.df.prior <- squeeze$df.prior
  out$gene.df.residual <- gene.df.residual
  out$gene.df.total <- gene.df.total
  out$gene.s2 <- gene.s2[gene.keep]
  out$gene.F <- gene.F
  out$gene.F.p.value <- gene.F.p.value

  # Which columns of exon.genes contain gene level annotation?
  gene.lastexon <- cumsum(gene.nexons)
  gene.firstexon <- gene.lastexon-gene.nexons+1
  no <- logical(nrow(exon.genes))
  isdup <- vapply(exon.genes,duplicated,no)[-gene.firstexon,,drop=FALSE]
  isgenelevel <- apply(isdup,2,all)
  out$gene.genes <- exon.genes[gene.lastexon,isgenelevel, drop=FALSE]
  out$gene.genes$NExons <- gene.nexons
  out$gene.firstexon <- gene.firstexon
  out$gene.lastexon <- gene.lastexon

  # Simes adjustment of exon level p-values
  penalty <- rep_len(1L,length(g))
  penalty[gene.lastexon] <- 1L-gene.nexons
  penalty <- cumsum(penalty)[-gene.lastexon]
  penalty <- penalty / rep(gene.nexons-1L,gene.nexons-1L)
  g2 <- g[-gene.lastexon]

  out$gene.simes.p.value <- gene.F.p.value
  for (j in seq_len(ncol(fit))) {
    o <- order(g,exon.p.value[,j])
    p.adj <- pmin(exon.p.value[o,j][-gene.lastexon] / penalty, 1)
    o <- order(g2,p.adj)
    out$gene.simes.p.value[,j] <- p.adj[o][gene.firstexon-0L:(ngenes-1L)]
  }

  out
}
