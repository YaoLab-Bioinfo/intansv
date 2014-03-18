
geneAnnotation <- function(structuralVariation, genomeAnnotation) 
{
    if (is.null(structuralVariation)|nrow(structuralVariation)<=0) {
        return(NULL)
    } else {
        structuralVariation <- structuralVariation[, 1:3]
        genomeAnnotationGeneIrange <- genomeAnnotation[genomeAnnotation$type!="gene"]
        structuralVariationIrange <- 
            GRanges(seqnames=structuralVariation$chromosome, 
                IRanges(start=structuralVariation$pos1, 
                        end=structuralVariation$pos2))
        structuralVariationIrange  <- reduce(structuralVariationIrange)
        structuralVariationGeneOp <- findOverlaps(genomeAnnotationGeneIrange, 
                                              structuralVariationIrange)
        strand(genomeAnnotationGeneIrange) <- "+"
        strand(structuralVariationIrange) <- "+"
        structuralVariationGeneOpDf <- NULL
        structuralVariationGeneOpDf$query <- queryHits(structuralVariationGeneOp)
        structuralVariationGeneOpDf$subject <- subjectHits(structuralVariationGeneOp)
        structuralVariationGeneOpDf <- as.data.frame(structuralVariationGeneOpDf)
        genomeAnnotationGeneDf <- NULL
        genomeAnnotationGeneDf$chr <- 
            rep(as.character((seqnames(genomeAnnotationGeneIrange))@values), 
                (seqnames(genomeAnnotationGeneIrange))@lengths)
        genomeAnnotationGeneDf$start <- start(genomeAnnotationGeneIrange)
        genomeAnnotationGeneDf$end <- end(genomeAnnotationGeneIrange)
        genomeAnnotationGeneDf$loc <- genomeAnnotationGeneIrange$ID
        genomeAnnotationGeneDf$name <- (genomeAnnotationGeneIrange$Parent)@unlistData

        genomeAnnotationGeneDf <- as.data.frame(genomeAnnotationGeneDf)
        genomeAnnotationGeneDf$query <- 1:nrow(genomeAnnotationGeneDf)
        genomeAnnotationGeneDfRes <- 
            merge(genomeAnnotationGeneDf, 
                  structuralVariationGeneOpDf, 
                  by="query")
        genomeAnnotationGeneDfRes$inter <- 
            width(pintersect(genomeAnnotationGeneIrange[genomeAnnotationGeneDfRes$query],
                         structuralVariationIrange[genomeAnnotationGeneDfRes$subject]))
        genomeAnnotationGeneDfRes$inter.per <- 
            genomeAnnotationGeneDfRes$inter/(genomeAnnotationGeneDfRes$end-
                                         genomeAnnotationGeneDfRes$start+1)
        genomeAnnotationGeneDfMer <- 
            merge(genomeAnnotationGeneDf, 
              genomeAnnotationGeneDfRes, 
              by=c("query", "chr", "start", "end", "loc", "name"), all=T)
        genomeAnnotationExon <- genomeAnnotation[genomeAnnotation$type=="exon"]
        genomeAnnotationGene <- genomeAnnotation[genomeAnnotation$type=="mRNA"]
        genomeAnnotationGeneDf <- NULL
        genomeAnnotationGeneDf$chr <- 
            rep(as.character((seqnames(genomeAnnotationGene))@values), 
                (seqnames(genomeAnnotationGene))@lengths)
        genomeAnnotationGeneDf$start <- (ranges(genomeAnnotationGene))@start
        genomeAnnotationGeneDf$end <- (ranges(genomeAnnotationGene))@start + 
                                  (ranges(genomeAnnotationGene))@width - 1
        genomeAnnotationGeneDf$name <- genomeAnnotationGene$ID
        genomeAnnotationGeneDf <- as.data.frame(genomeAnnotationGeneDf, 
                                            stringAsfactor=F)
        names(genomeAnnotationGeneDf)[2:3] <- c("gene.start", "gene.end")
        genomeAnnotationExonDf <- NULL
        genomeAnnotationExonDf$chr <- 
            rep(as.character((seqnames(genomeAnnotationExon))@values), 
               (seqnames(genomeAnnotationExon))@lengths)
        genomeAnnotationExonDf$start <- (ranges(genomeAnnotationExon))@start
        genomeAnnotationExonDf$end <- (ranges(genomeAnnotationExon))@start +
                                  (ranges(genomeAnnotationExon))@width - 1
        genomeAnnotationExonDf$name <- (genomeAnnotationExon$Parent)@unlistData
        genomeAnnotationExonDf <- as.data.frame(genomeAnnotationExonDf, 
                                            stringAsfactor=F)
        names(genomeAnnotationExonDf)[2:3] <- c("exon.start", "exon.end")
        genomeAnnotationExonIrange <- 
            GRanges(seqnames=genomeAnnotationExonDf$name,
                IRanges(genomeAnnotationExonDf$exon.start, 
                        genomeAnnotationExonDf$exon.end))
        genomeAnnotationGeneIrange <- 
            GRanges(seqnames=genomeAnnotationGeneDf$name,
                IRanges(genomeAnnotationGeneDf$gene.start, 
                        genomeAnnotationGeneDf$gene.end))
        genomeAnnotationIntronIrange <- setdiff(genomeAnnotationGeneIrange, 
            genomeAnnotationExonIrange)
        locus.chr <- as.character(genomeAnnotationGeneDf$chr)
        names(locus.chr) <- as.character(genomeAnnotationGeneDf$name)
        genomeAnnotationGeneDf$strand <- rep(
            as.character((strand(genomeAnnotationGene))@values), 
            (strand(genomeAnnotationGene))@lengths)
        locus.strand <- as.character(genomeAnnotationGeneDf$strand)
        names(locus.strand) <- as.character(genomeAnnotationGeneDf$name)
        genomeAnnotationExonDf$strand <- locus.strand[genomeAnnotationExonDf$name]
        genomeAnnotationIntronDf <- NULL
        genomeAnnotationIntronDf$intron.start <- start(genomeAnnotationIntronIrange)
        genomeAnnotationIntronDf$intron.end <- end(genomeAnnotationIntronIrange)
        genomeAnnotationIntronDf$name <- rep(
            as.character((seqnames(genomeAnnotationIntronIrange))@values), 
            (seqnames(genomeAnnotationIntronIrange))@lengths)
        genomeAnnotationIntronDf <- as.data.frame(genomeAnnotationIntronDf, 
                                              stringAsfactor=F)
        genomeAnnotationIntronDf$strand <- 
            locus.strand[genomeAnnotationIntronDf$name]
        genomeAnnotationIntronDfPlus <- 
            genomeAnnotationIntronDf[genomeAnnotationIntronDf$strand=="+", ]
        genomeAnnotationIntronDfMinus <- 
            genomeAnnotationIntronDf[genomeAnnotationIntronDf$strand=="-", ]
        intronPlusNum <- (rle(as.character(genomeAnnotationIntronDfPlus$name)))$lengths
        intronMinusNum <- (rle(as.character(genomeAnnotationIntronDfMinus$name)))$lengths
        genomeAnnotationIntronDfPlus$id <- 
            unlist(sapply(intronPlusNum, function(x){return(c(1:x))}))
        genomeAnnotationIntronDfMinus$id <- 
            unlist(sapply(intronMinusNum, function(x){return(rev(c(1:x)))}))
        genomeAnnotationIntronDf <- rbind(genomeAnnotationIntronDfMinus, 
                                      genomeAnnotationIntronDfPlus)
        genomeAnnotationIntronDf$id <- 
            paste("intron_", genomeAnnotationIntronDf$id, sep="")
        genomeAnnotationIntronDf$chr <- 
            locus.chr[as.character(genomeAnnotationIntronDf$name)]
        genomeAnnotationIntron <- GRanges(seqnames=genomeAnnotationIntronDf$chr, 
            IRanges(start=genomeAnnotationIntronDf$intron.start, 
            end=genomeAnnotationIntronDf$intron.end), 
            name=genomeAnnotationIntronDf$name, id=genomeAnnotationIntronDf$id)
        strand(genomeAnnotationIntron) <- "+"
        genomeAnnotationIntronOp <- findOverlaps(genomeAnnotationIntron, 
                                             structuralVariationIrange)
        genomeAnnotationIntronOpDf <- NULL
        genomeAnnotationIntronOpDf$query <- queryHits(genomeAnnotationIntronOp)
        genomeAnnotationIntronOpDf$subject <- subjectHits(genomeAnnotationIntronOp)
        genomeAnnotationIntronOpDf <- as.data.frame(genomeAnnotationIntronOpDf)

        genomeAnnotationIntronDf$query <- 1:nrow(genomeAnnotationIntronDf)
        genomeAnnotationIntronOpDfMer <- merge(genomeAnnotationIntronOpDf, 
                                           genomeAnnotationIntronDf, 
                                           by="query")
        genomeAnnotationIntronOpDfMer$inter <- width(pintersect(
            genomeAnnotationIntron[genomeAnnotationIntronOpDfMer$query], 
            structuralVariationIrange[genomeAnnotationIntronOpDfMer$subject]))
        genomeAnnotationIntronOpDfMer$inter.per <- 
            genomeAnnotationIntronOpDfMer$inter/(genomeAnnotationIntronOpDfMer$intron.end-
            genomeAnnotationIntronOpDfMer$intron.start+1)

        genomeAnnotationIntronDfMer <- 
            merge(genomeAnnotationIntronDf, 
              genomeAnnotationIntronOpDfMer, 
              by=c("intron.start", "intron.end", 
                   "name", "strand", "id", "chr", "query"), all=T)

        genomeAnnotationGeneDfMer$query <- NULL
        genomeAnnotationGeneDfMer$subject <- NULL
        names(genomeAnnotationGeneDfMer)[4] <- "id"
        names(genomeAnnotationGeneDfMer)[5] <- "loc"
        genomeAnnotationGeneDfMer <- genomeAnnotationGeneDfMer[, c(1:3, 7, 4:6)]
    
        genomeAnnotationIntronDfMer$query <- NULL
        genomeAnnotationIntronDfMer$subject <- NULL
        genomeAnnotationIntronDfMer$strand <- NULL
        names(genomeAnnotationIntronDfMer) <- c("start", "end", "loc", "id", 
                                            "chr", "inter", "inter.per")
        genomeAnnotationIntronDfMer <- 
        genomeAnnotationIntronDfMer[, c(5, 1:2, 7, 4, 3, 6)]

        GeneAnno <- rbind(genomeAnnotationGeneDfMer, genomeAnnotationIntronDfMer)
        GeneAnno$type <- NA
        GeneAnno[grepl("exon", GeneAnno$id, ignore.case=TRUE), ]$type <- "exon"
        GeneAnno[grepl("intron", GeneAnno$id, ignore.case=TRUE), ]$type <- "intron"
        GeneAnno[grepl("cds", GeneAnno$id, ignore.case=TRUE), ]$type <- "cds"
        GeneAnno[grepl("utr", GeneAnno$id, ignore.case=TRUE), ]$type <- "utr"
        GeneAnno <- GeneAnno[GeneAnno$type%in%c("exon", "intron", "cds", "utr"), ]
        GeneAnno <- GeneAnno[order(GeneAnno$loc, GeneAnno$type, GeneAnno$id), ]

        GeneAnno$inter <- NULL
        GeneAnno$id <- as.character(GeneAnno$id)
        GeneAnno$id[is.na(GeneAnno$id)] <- 0
        GeneAnno$inter.per[is.na(GeneAnno$inter.per)] <- 0
        GeneAnno$inter.per <- round(GeneAnno$inter.per, 2)
        GeneAnnoRes <- ddply(GeneAnno, ("loc"), GeneReport)
        names(GeneAnnoRes) <- c("locus", "exon", "intron", "cds", "utr")
        GeneAnnoRes$exon[GeneAnnoRes$exon==""] <- NA
        GeneAnnoRes$intron[GeneAnnoRes$intron==""] <- NA
        GeneAnnoRes$cds[GeneAnnoRes$cds==""] <- NA
        GeneAnnoRes$utr[GeneAnnoRes$utr==""] <- NA
        return(GeneAnnoRes)
    }
}


GeneReport <- function(df)
{
    if(all(df$inter.per==0)) {
        return(NULL)
    } else {
        df.exon <- df[df$type=="exon", ]
        dfIntron <- df[df$type=="intron", ]
        df.cds <- df[df$type=="cds", ]
        df.utr <- df[df$type=="utr", ]
        return(c(paste(df.exon$inter.per, sep="", collapse=":"), 
            paste(dfIntron$inter.per, sep="", collapse=":"), 
            paste(df.cds$inter.per, sep="", collapse=":"), 
            paste(df.utr$inter.per, sep="", collapse=":")))
    }
}


