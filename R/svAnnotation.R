
svAnnotation <- function(structuralVariation, genomeAnnotation) 
{
  if (!is.data.frame(structuralVariation)) {
    stop("structuralVariation should be a data frame!\n")
  } else if (!identical(names(structuralVariation)[1:3], 
                        c("chromosome", "pos1", "pos2"))) {
    stop("structuralVariation should have appropriate names!\n")
  } else if (class(genomeAnnotation)!="GRanges") {
    stop("genomeAnnotation should be a Genomic Range\n")
  } else if (!(.hasSlot(genomeAnnotation, "elementMetadata") && 
                 .hasSlot(genomeAnnotation@elementMetadata, "listData") && 
                 any(names(genomeAnnotation@elementMetadata)=="type") && 
                 any(names(genomeAnnotation@elementMetadata)=="Parent") )) {
    stop("The content of genomeAnnotation is not appropriate!\n")
  }
  
  if (is.null(structuralVariation)|nrow(structuralVariation)<=0) {
    return(NULL)
  } else {
    structuralVariation <- structuralVariation[, 1:3]
    
    structuralVariation$start <- as.numeric(structuralVariation$pos1)
    structuralVariation$end <- as.numeric(structuralVariation$pos2)
    structuralVariationIrange <- 
      GRanges(seqnames=structuralVariation$chromosome, 
              ranges=IRanges(start=structuralVariation$pos1, 
                             end=structuralVariation$pos2))
    geneAnnoRes <- findOverlaps(structuralVariationIrange, genomeAnnotation)
    structuralVariation$query <- 1:nrow(structuralVariation)
    geneAnnoResDf <- as.data.frame(cbind(queryHits(geneAnnoRes), 
                                         subjectHits(geneAnnoRes)))
    names(geneAnnoResDf) <- c("query", "subject")
    structuralVariationGeneAnno <- merge(structuralVariation, 
                                         geneAnnoResDf, by="query", all=T)
    structuralVariationGeneAnnoNa <- 
      structuralVariationGeneAnno[is.na(structuralVariationGeneAnno$query) |
                                    is.na(structuralVariationGeneAnno$subject), ]
    structuralVariationGeneAnnoNna <- 
      structuralVariationGeneAnno[!(is.na(structuralVariationGeneAnno$query)) &
                                    !(is.na(structuralVariationGeneAnno$subject)), ]
    structuralVariationGeneAnnoNna$overlap <- sprintf("%.3f", 
                                                      width(pintersect(structuralVariationIrange[structuralVariationGeneAnnoNna$query], 
                                                                       genomeAnnotation[structuralVariationGeneAnnoNna$subject], ignore.strand=T))/
                                                        width(genomeAnnotation[structuralVariationGeneAnnoNna$subject]))
    structuralVariationGeneAnnoNna$annotation <- 
      genomeAnnotation[as.numeric(structuralVariationGeneAnnoNna$subject)]$type
    structuralVariationGeneAnnoNna$id <- 
      genomeAnnotation[as.numeric(structuralVariationGeneAnnoNna$subject)]$ID
    structuralVariationGeneAnnoNna$parent <- NA
    structuralVariationGeneAnnoNna[structuralVariationGeneAnnoNna$annotation=="gene", ]$parent <- 
      genomeAnnotation[as.numeric(structuralVariationGeneAnnoNna[
        structuralVariationGeneAnnoNna$annotation=="gene", ]$subject)]$Name
    structuralVariationGeneAnnoNna[structuralVariationGeneAnnoNna$annotation!="gene", ]$parent <- 
      (genomeAnnotation[as.numeric(structuralVariationGeneAnnoNna[
        structuralVariationGeneAnnoNna$annotation!="gene", ]$subject)]$Parent)@unlistData
    
    if (nrow(structuralVariationGeneAnnoNa)>0) {
      structuralVariationGeneAnnoNa$parent <- NA
      structuralVariationGeneAnnoNa$id <- NA
      structuralVariationGeneAnnoNa$annotation <- "intergenic"
      structuralVariationGeneAnnoNa$overlap <- NA
    }
    
    genomeAnnotationExon <- genomeAnnotation[genomeAnnotation$type=="exon"]
    genomeAnnotationGene <- genomeAnnotation[genomeAnnotation$type=="mRNA"]
    genomeAnnotationGeneDf <- NULL
    genomeAnnotationGeneDf$chr <- 
      rep(as.character((seqnames(genomeAnnotationGene))@values), 
          (seqnames(genomeAnnotationGene))@lengths)
    genomeAnnotationGeneDf$start <- (ranges(genomeAnnotationGene))@start
    genomeAnnotationGeneDf$end <- (ranges(genomeAnnotationGene))@start + 
      (ranges(genomeAnnotationGene))@width - 1
    genomeAnnotationGeneDf$name <- as.character(genomeAnnotationGene$ID)
    genomeAnnotationGeneDf <- as.data.frame(genomeAnnotationGeneDf, 
                                            stringsAsFactors=FALSE)
    
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
    genomeAnnotationIntronIrange <- 
      setdiff(genomeAnnotationGeneIrange, genomeAnnotationExonIrange)
    LocusChr <- as.character(genomeAnnotationGeneDf$chr)
    names(LocusChr) <- as.character(genomeAnnotationGeneDf$name)
    
    
    genomeAnnotationGeneDf$strand <- 
      rep(as.character((strand(genomeAnnotationGene))@values), 
          (strand(genomeAnnotationGene))@lengths)
    LocusStrand <- as.character(genomeAnnotationGeneDf$strand)
    names(LocusStrand) <- as.character(genomeAnnotationGeneDf$name)
    genomeAnnotationExonDf$strand <- LocusStrand[genomeAnnotationExonDf$name]
    genomeAnnotationIntronDf <- NULL
    genomeAnnotationIntronDf$intron.start <- start(genomeAnnotationIntronIrange)
    genomeAnnotationIntronDf$intron.end <- end(genomeAnnotationIntronIrange)
    genomeAnnotationIntronDf$name <- 
      rep(as.character((seqnames(genomeAnnotationIntronIrange))@values), 
          (seqnames(genomeAnnotationIntronIrange))@lengths)
    genomeAnnotationIntronDf <- 
      as.data.frame(genomeAnnotationIntronDf, stringAsfactor=F)
    genomeAnnotationIntronDf$strand <- LocusStrand[genomeAnnotationIntronDf$name]
    genomeAnnotationIntronDfPlus <- 
      genomeAnnotationIntronDf[genomeAnnotationIntronDf$strand=="+", ]
    genomeAnnotationIntronDfMinus <- 
      genomeAnnotationIntronDf[genomeAnnotationIntronDf$strand=="-", ]
    
    intronPlusNum <- (rle(as.character(genomeAnnotationIntronDfPlus$name)))$lengths
    intronMinusNum <- (rle(as.character(genomeAnnotationIntronDfMinus$name)))$lengths
    
    genomeAnnotationIntronDfPlus$id <- unlist(sapply(intronPlusNum, 
                                                     function(x){return(c(1:x))}))
    genomeAnnotationIntronDfMinus$id <- unlist(sapply(intronMinusNum, 
                                                      function(x){return(rev(c(1:x)))}))
    genomeAnnotationIntronDf <- rbind(genomeAnnotationIntronDfMinus, 
                                      genomeAnnotationIntronDfPlus)
    genomeAnnotationIntronDf$id <- paste("intron_", 
                                         genomeAnnotationIntronDf$id, sep="")
    
    genomeAnnotationIntronDf$chr <- LocusChr[as.character(genomeAnnotationIntronDf$name)]
    genomeAnnotationIntron <- GRanges(seqnames=genomeAnnotationIntronDf$chr, 
                                      IRanges(start=genomeAnnotationIntronDf$intron.start, 
                                              end=genomeAnnotationIntronDf$intron.end), 
                                      name=genomeAnnotationIntronDf$name, 
                                      id=genomeAnnotationIntronDf$id)
    intronAnnoRes <- findOverlaps(structuralVariationIrange, 
                                  genomeAnnotationIntron)
    intronAnnoResDf <- as.data.frame(cbind(queryHits(intronAnnoRes), 
                                           subjectHits(intronAnnoRes)))
    names(intronAnnoResDf) <- c("query", "subject")
    structuralVariationIntronAnno <- merge(structuralVariation, 
                                           intronAnnoResDf, by="query", all=T)
    structuralVariationIntronAnnoNna <- 
      structuralVariationIntronAnno[!(is.na(structuralVariationIntronAnno$query)) &
                                      !(is.na(structuralVariationIntronAnno$subject)), ]
    structuralVariationIntronAnnoNna$overlap <- 
      sprintf("%.3f", width(pintersect(structuralVariationIrange[structuralVariationIntronAnnoNna$query], 
                                       genomeAnnotationIntron[structuralVariationIntronAnnoNna$subject])) / 
                width(genomeAnnotationIntron[structuralVariationIntronAnnoNna$subject]))
    structuralVariationIntronAnnoNna$annotation <- "intron"
    
    structuralVariationIntronAnnoNna$id <- 
      genomeAnnotationIntron[as.numeric(structuralVariationIntronAnnoNna$subject)]$id
    structuralVariationIntronAnnoNna$parent <- 
      genomeAnnotationIntron[as.numeric(structuralVariationIntronAnnoNna$subject)]$name
    structuralVariationAnno <- 
      rbind(structuralVariationGeneAnnoNna, 
            structuralVariationGeneAnnoNa, structuralVariationIntronAnnoNna)
    structuralVariationAnnoOrdered <- 
      structuralVariationAnno[order(structuralVariationAnno$chr, 
                                    structuralVariationAnno$start, 
                                    structuralVariationAnno$parent), ]
    rownames(structuralVariationAnnoOrdered) <- 1:nrow(structuralVariationAnnoOrdered)
    structuralVariationAnnoOrdered$query <- NULL
    structuralVariationAnnoOrdered$subject <- NULL
    structuralVariationAnnoOrdered$start <- NULL
    structuralVariationAnnoOrdered$end <- NULL
    return(structuralVariationAnnoOrdered)
  }
}

