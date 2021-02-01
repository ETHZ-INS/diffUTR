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

.geneLevelStats <- function(d, gene.qval=NULL){
  stopifnot(c("bin.pval","coef","gene","width","meanLogDensity") %in% colnames(d))
  if(is.null(gene.qval)) gene.qval <- simes.aggregation(d$bin.pval, d$gene)
  si <- split(seq_len(nrow(d)), d$gene)
  d2 <- data.frame( row.names=names(si))
  d$log10p <- -log10(d$bin.pval)
  d2 <- vapply( si, FUN.VALUE=numeric(8), FUN=function(i){
    w <- d$log10p[i]
    if(sum(w)==0) w <- rep(1,length(i))
    c( w.coef=weighted.mean(d$coef[i], w),
       w.abs.coef=weighted.mean(abs(d$coef[i]), w),
       w.width=sum(d$width[i]*w)/sum(w),
       w.density=log(weighted.mean(d$meanLogDensity[i], w)),
       sizeScore=weighted.mean(d$coef[i], w*d$width[i]),
       abs.sizeScore=weighted.mean(abs(d$coef[i]), w*d$width[i]),
       geneMeanDensity=mean(d$meanLogDensity[i], trim=0.05),
       density.ratio=weighted.mean(d$meanLogDensity[i]-mean(d$meanLogDensity[i]), w)
    )
  })
  d2 <- as.data.frame(t(d2))
  for(f in c(1:2,5:8)) d2[,f] <- round(d2[,f],3)
  d2$w.width <- as.integer(round(d2$w.width))
  d2$nb.bins <- lengths(si)
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
