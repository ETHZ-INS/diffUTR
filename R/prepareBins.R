#' prepareBins
#'
#' @param g A GRanges (or path to RDS file containing a GRanges) or path to a
#' gtf file or EnsDb object containing the gene annotation.
#' @param APA A GRanges (or path to a GRanges in RDS format) or bed file
#' containing the alternative poly-A site database
#' @param onlyMainChr Logical; whether to keep only main chromosomes
#' @param removeAntisense Logical; whether to remove antisense APA sites
#' @param chrStyle Chromosome notation to convert to, default "UCSC"
#' @param maxUTRbinSize Max width of new alternative UTR bins
#' @param genewise Logical, whether annotation should be flattened genewise
#' @param stranded Logical, whether to perform disjoin in a stranded fashion.
#' @param codingOnly Logical, whether to keep only coding transcripts
#' @param verbose Logical, whether to print run information
#'
#' @details See the vignette for more details.
#'
#' @return A `GRanges` object.
#'
#' @author Stefan Greber
#' @export
#' @import GenomicRanges rtracklayer GenomeInfoDb
#' @importFrom S4Vectors mcols
#' @importFrom IRanges FactorList CharacterList extractList
#' @examples
#' data(example_gene_annotation)
#' bins <- prepareBins(example_gene_annotation)
prepareBins <- function( g, APA=NULL, onlyMainChr=TRUE, removeAntisense=TRUE,
                         chrStyle="UCSC", maxUTRbinSize=15000,
                         codingOnly=FALSE,
                         genewise=FALSE, stranded=FALSE, verbose=TRUE ){

  if(verbose) message("Preparing annotation")
  if(!is.null(APA)) APA <- .prepAPA(APA)
  g <- .prepGeneAnno(g, verbose=verbose)

  if(onlyMainChr){
    #Extract only main chromosomes
    g <- keepStandardChromosomes(g, pruning.mode="coarse")
    if(!is.null(APA))
      APA <- keepStandardChromosomes(APA, pruning.mode="coarse")
  }
  seqlevelsStyle(g) <- chrStyle
  seqlevels(g) <- seqlevelsInUse(g)

  if(!is.null(APA)){
    if(verbose) message("Generating alternative UTR bins")
    seqlevelsStyle(APA) <- chrStyle
    seqlevels(APA) <- seqlevelsInUse(APA)
    g <- .extendWithAPA(g, APA)
  }

  if(codingOnly) g <- g[g$tx_biotype=="protein_coding"]

  if(verbose) message("Merging and disjoining bins")

  if(genewise){
    bins <- .disjoinGeneWise(g)
  } else{
    bins <- disjoin(g, with.revmap=TRUE, ignore.strand=!stranded)
    f <- rep(seq_along(bins), lengths(bins$revmap))
    mcols(bins) <- DataFrame(lapply( mcols(g)[unlist(bins$revmap),],
                                     FUN=function(x){
                                       sort(unique(splitAsList(x,f)))
                                     }))
  }
  rm(g)

  bt <- CharacterList(bins$type)
  if(is.null(APA)){
    bt <- as.factor(ifelse(any(bt=="CDS"), "CDS", "UTR"))
  }else{
    #concatenate the lists with a + as seperator so each combination of
    # transcripts or type has a unique string as identifier
    #bins$type <- unstrsplit(CharacterList(bins$type),sep="+")
    bt <- 100L*any(bt=="CDS") + 10L*any(bt=="UTR") +
      as.integer(any(bt=="3UTR"))
    dict <- c("111"="CDS/UTR/3UTR", "110"="CDS/UTR", "100"="CDS",
              "11"="UTR/3UTR", "10"="UTR", "1"="3UTR")
    bt <- factor(bt, as.integer(names(dict)), as.character(dict))
  }
  bins$type <- bt
  rm(bt)

  levels(bins$type) <- c(levels(bins$type), "non-coding")
  bins$tx <- unstrsplit(CharacterList(bins$tx), sep="+")

  # duplicate ranges with more then one overlapping gene so each gene has a
  # separate range

  #find the number of genes overlapping this bin
  len<-lengths(bins$gene)
  #for each bin repeat its index as many times as genes overlapping it
  repid<-rep(seq_along(bins),len)
  #save gene_id list for after
  gene_id=bins$gene
  #index bins with index created before in order to multiply the bins as many
  #times as genes overlapping it
  bins<-bins[repid]
  #overwrite the gene_id so that each multipl of a bin has only one
  bins$gene<-unlist(gene_id)
  for(f in c("gene_name", "tx_biotype")){
    if(f %in% colnames(mcols(bins)) && !is.vector(mcols(bins)[[f]]))
      mcols(bins)[[f]] <- CharacterList(mcols(bins)[[f]])
  }
  if(is(bins$tx_biotype,"AtomicList")){
    bins$type[!any(bins$tx_biotype=="protein_coding")] <- "non-coding"
  }

  bins <- sort(bins, by=~gene)
  bins$bin_id <- ave(as.character(strand(bins)),bins$gene,FUN=function(x){
    res<-seq_along(x)
    if(x[1]=="-"){
      res<-rev(res)
    }
    return(res)
  })
  names(bins)<-NULL

  for(f in c("gene_biotype", "tx_biotype")){
    if(f %in% colnames(mcols(bins)))
      mcols(bins)[[f]] <- FactorList(mcols(bins)[[f]])
  }
  for(f in c("tx", "gene")){
    if(f %in% colnames(mcols(bins)))
      mcols(bins)[[f]] <- as.factor(mcols(bins)[[f]])
  }

  o <- findOverlaps(bins, ignore.strand=!stranded)
  o <- o[from(o) %in% unique(from(o)[duplicated(from(o))])]
  ga <- rowsum(as.integer(bins$gene[from(o)] != bins$gene[to(o)]), from(o))
  bins$geneAmbiguous <- FALSE
  bins$geneAmbiguous[as.integer(row.names(ga))] <- ga[,1]>0

  bins
}

#' @importFrom ensembldb genes cdsBy exons
.prepGeneAnno <- function(g, verbose=TRUE){
  if(is.character(g) && length(g)==1) {
    if (grepl("\\.gtf$", g)) {
      g <- rtracklayer::import(g)
    }else{
      stop("Annotation needs to be provided as gtf")
    }
  }

  if(isEnsDb <- is(g,"EnsDb")){
    if(verbose) message("Extracting exon information from EnsDb")
    e <- exons(g, c("exon_id","tx_id","gene_id","gene_name",
                        "gene_biotype","tx_biotype"))
    e$type <- "exon"
    cd <- unlist(cdsBy(g,
                       columns=c("tx_id","tx_biotype","gene_id","gene_name")))
    cd$type <- "CDS"
    g <- genes(g, c("gene_id","gene_name"))
    g$type <- "gene"
    g <- sort(c(g,e,cd))
    rm(e,cd)
  }
  if(!is(g,"GRanges")) stop("Annotation is not valid.")

  synonyms <- list( type="type",
                    tx=c("tx_id","tx_name","transcript","transcript_id",
                         "transcript_name","tx"),
                    gene=c("gene_id","gene","gene_name","gene_symbol",
                           "symbol"),
                    tx_biotype=c("tx_biotype","transcript_biotype",
                                 "transcript_type"),
                    gene_biotype=c("gene_biotype","gene_type"))
  for(f in names(synonyms)){
    ff <- intersect(synonyms[[f]], colnames(mcols(g)))
    if(length(ff)==0) stop(f," information not present!")
    mcols(g)[[f]] <- mcols(g)[[ff[1]]]
    if(f!=ff[1]) mcols(g)[[ff[[1]]]] <- NULL
  }

  if(!isEnsDb){
    if(!any(g$type=="gene") && any(g$type=="transcript")){
      if(!is.factor(g$type)) g$type <- as.factor(g$type)
      levels(g$type) <- unique(levels(g$type), "gene")
      # build the gene entries
      g2 <- g[g$type=="transcript"]
      ff <- c("tx","tx_biotype")
      ff <- intersect(unique(c(ff,unlist(synonyms[ff]))), colnames(mcols(g2)))
      for(f in ff) mcols(g2)[[f]] <- NA
      g2$type <- "gene"
      gmc <- mcols(g2)[!duplicated(g2$gene),]
      row.names(gmc) <- gmc$gene
      g3 <- unlist(reduce(split(g2, g2$gene)))
      mcols(g3) <- gmc[names(g3),]
      g <- c(g,g3)
    }
    g <- g[g$type %in% c("gene", "exon", "CDS")]
    g$type <- droplevels(g$type)
  }
  g
}

.prepAPA <- function(apa, removeAntisense=TRUE){
  if(is.character(apa) && length(apa)==1) {
    if(grepl("\\.rds$", apa)) {
      apa <- readRDS(apa)
      if(!is(apa, "GRanges")) stop("If a RDS file, apa should be a GRanges")
    }else if(grepl("\\.bed$|\\.bed.gz$", apa)){
      apa <- tryCatch(rtracklayer::import.bed(apa), error=function(e){
        extraColvect <- c(percentage="numeric", numberofprot="integer",
                          tpm2="numeric",encod="character",addinfo="character")
        rtracklayer::import("apa.mm38.bed.gz", format="bed",
                            extraCols=extraColvect)
      })
    }
  }else if (!is(apa, "GRanges")) {
    stop("Annotation has to be provided as GRanges or bed")
  }
  if("encod" %in% colnames(mcols(apa)))
    apa <- apa[!(apa$encod %in% c("AE","AI","AU"))]
  apa
}

#' @importFrom IRanges IRanges
.extendWithAPA <- function(g, APA, maxUTRbinSize=15000){
  # get regions which could be boundrys for bins
  gtf1<-g[g$type %in% c("gene", "exon"),]
  #resize gene to 1 as only start of a gene would be a boundry
  gtf1[gtf1$type=="gene",]<-resize(gtf1[gtf1$type=="gene",],1)

  # get regions which could be boundrys for bins

  #seperate by strand because of slightly different syntax
  APAplus<-APA[(as.character(strand(APA))=="+")]
  APAminus<-APA[(as.character(strand(APA))=="-")]

  # find boundry region before each APA site...some have none,
  # so ignore these and do it again
  outmin<-follow(APAminus, gtf1)
  APAminus<-APAminus[is.na(outmin)==FALSE]
  outmin<-follow(APAminus, gtf1, select="all")
  APAminus<-APAminus[from(outmin)]

  outplus<-follow(APAplus, gtf1)
  APAplus<-APAplus[is.na(outplus)==FALSE]
  outplus<-follow(APAplus, gtf1, select="all")
  APAplus<-APAplus[from(outplus)]


  minusranges<-gtf1[to(outmin),]
  plusranges<-gtf1[to(outplus),]

  #create regions from boundry to APA and include metadata

  beginplusbin<-end(plusranges)+1
  endplusbin<-end(APAplus)

  beginminusbin<-start(APAminus)
  endminusbin<-start(minusranges)-1

  plusbin<-GRanges(seqnames=seqnames(plusranges),
                   IRanges(start=beginplusbin, end=endplusbin),
                   strand=as.character(strand(plusranges)))

  minusbin<-GRanges(seqnames=seqnames(minusranges),
                    IRanges(start=beginminusbin, end=endminusbin),
                    strand=as.character(strand(minusranges)))
  mcols(plusbin)<-mcols(plusranges)
  mcols(minusbin)<-mcols(minusranges)

  #merge all regions necessary for differential analysis
  a <- c(plusbin,minusbin)
  a$type<-"UTR"
  #threshold
  a <- a[width(a)<maxUTRbinSize,]

  b <- GenomicRanges::setdiff(g[g$type=="exon"], g[g$type=="CDS"])

  o<-as(findOverlaps(b,g[g$type=="exon"]),"List")
  ex<-extractList(mcols(g[g$type=="exon"]), o)
  b <- rep(b, lengths(o))
  mcols(b)<-unlist(ex)
  if(is.factor(b$type)) levels(b$type) <- unique(c(levels(b$type),"UTR"))
  b$type[b$tx_biotype=="protein_coding"] <- "UTR"
  b<-c(b,g[g$type=="CDS"])
  # #b$type[b$type=="CDS"]<-"exon"

  cds<-g[g$type=="CDS"]

  a<-a[!is.na(a$tx)]

  cdsminus<-cds[as.character(strand(cds))=="-",]
  cdsplus<-cds[as.character(strand(cds))=="+",]

  #split by transcripts
  cdsplus<- split(cdsplus, cdsplus$tx)
  cdsminus<- split(cdsminus, cdsminus$tx)

  #find direct index of the last cds of each transcript and get them
  l<-elementNROWS(cdsplus)
  l<-cumsum(l)
  lastcdsplus<-unlist(sort(cdsplus))[l,]

  l<-elementNROWS(cdsminus)
  l<-cumsum(l)-l+1
  lastcdsminus<-unlist(sort(cdsminus))[l,]
  lastCDS<-c(lastcdsminus,lastcdsplus)

  end<-data.frame(tx=lastCDS$tx,end=end(lastCDS))
  end$end[as.character(strand(lastCDS))=="-"] <-
    start(lastCDS)[as.character(strand(lastCDS))=="-"]

  a <- c(a,b)
  o <- order(a$tx)
  a<-a[o,]
  o <- order(end$tx)
  end<-end[o,]

  tx.with.cds <- a$tx %in% end$tx

  tx.nexons <- rowsum(rep(1,length(a$tx[tx.with.cds])),
                      group=a$tx[tx.with.cds], reorder = FALSE)
  ntx <- length(unique(end$tx))
  tidx <- rep(seq_len(ntx), times=tx.nexons)
  a$type <- as.character(a$type)

  # UTR that is 3'

  nocds<-a[!tx.with.cds]
  a<-a[tx.with.cds]
  a$type[a$type=="UTR" & as.character(strand(a))=="+" &
           end$end[tidx]<start(a)]<-"3UTR"
  a$type[a$type=="UTR" & as.character(strand(a))=="-" &
           end$end[tidx]>end(a)]<-"3UTR"
  c(a,nocds)
}

.disjoinGeneWise <- function(a){
  #disjoin them to counting bins
  a<-split(a,a$gene)
  #disjoin them to counting bins
  bins<-disjoin(a,with.revmap=TRUE)

  #reverse mapping indexes are indexes of the corresponding list not of the
  #unlisted bins, so find the length of all sublists
  #for each sublist find the total number of previous features
  g<-cumsum(c(0,head(lengths(a),-1)))

  #this needs to be added to every revmap index of features that were created
  #from this sublist, so expand by the number of features

  adder<-rep(g,lengths(bins))
  bins<-unlist(bins)

  #and as there are list of indexes, expand by the number of indexes
  adder<-rep(adder,lengths(bins$revmap))
  trueidx<-unlist(bins$revmap)+adder

  #get back metadata, tx_name as string with name of each overlapping tx
  #name seperated by + gene_id as List of all genes overlapping this bin
  #type as string: with types overlapping this bin sperated by +
  #booltype->is utr or not


  # for each bin create as many slots as there were overlapping regions and
  # give them the index of the bin
  f <- rep(seq_along(bins), lengths(bins$revmap))

  #put the overlapping regions in to the slots and split them as list by the
  #index (bin they belong to), apply sort and unique in order to only get the
  #metadata once
  mcols(bins) <- DataFrame(lapply( mcols(unlist(a))[trueidx,],
                                FUN=function(x){
                                  sort(unique(splitAsList(x,f)))
                                }))
  bins
}
