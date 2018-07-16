
## Displaying the distribution of SVs in the whole genome
plotChromosome <- function(genome, structuralVariation, 
                           windowSize=1000000) 
{
  genome.gr <- GRanges(genome$chr, IRanges(genome$start, genome$end))
  
  genomeRes <- ddply(genome, c("chr", "end"), function(df){
    if (df$end[1]>windowSize) {
      pos1 <- seq(1, df$end[1], by=as.numeric(windowSize))
      pos2 <- seq(as.numeric(windowSize), df$end[1], by=as.numeric(windowSize))
      if(length(pos1)!=length(pos2)){
        pos2<-c(pos2, df$end[1])
      }
    } else {
      pos1 <- 1; pos2 <- df$end[1];
    }
    dfRes <- as.data.frame(cbind(pos1, pos2))
    return(dfRes)
  })
  genomeResIrange <- GRanges(seqnames=genomeRes$chr, 
                             IRanges(start=genomeRes$pos1, 
                                     end=genomeRes$pos2))
  genomeRes$query <- 1:nrow(genomeRes)
  
  if ( !is.null(structuralVariation$del) && (is.data.frame(structuralVariation$del)) &&
       (nrow(structuralVariation$del)>0) ) {
    delDf <- structuralVariation$del[(structuralVariation$del)$chromosome%in%
                                       genome$chr, ]
    delIrange <- GRanges(seqnames=delDf$chromosome, 
                         IRanges(start=as.numeric(delDf$pos1), 
                                 end=as.numeric(delDf$pos2)))
    delOverlap <- findOverlaps(genomeResIrange, delIrange)
    delOverlapRes <- as.data.frame(cbind(queryHits(delOverlap), 
                                         subjectHits(delOverlap)))
    names(delOverlapRes) <- c("query", "delTag")
    delOverlapCount <- ddply(delOverlapRes, ("query"), nrow)
    names(delOverlapCount)[2] <- "delTag"
    genomeDel <- merge(genomeRes, delOverlapCount, by="query", all=T)
  } else {
    genomeDel <- genomeRes
    delOverlapCount <- NULL
  }
  
  if ( !is.null(structuralVariation$dup) && (is.data.frame(structuralVariation$dup)) && 
       (nrow(structuralVariation$dup)>0) ) {
    dupDf <- structuralVariation$dup[(structuralVariation$dup)$chromosome %in%
                                       genome$chr, ]
    dupIrange <- GRanges(seqnames=dupDf$chromosome, 
                         IRanges(start=as.numeric(dupDf$pos1), 
                                 end=as.numeric(dupDf$pos2)))
    dupOverlap <- findOverlaps(genomeResIrange, dupIrange)    
    dupOverlapRes <- as.data.frame(cbind(queryHits(dupOverlap), 
                                         subjectHits(dupOverlap)))
    names(dupOverlapRes) <- c("query", "subject")
    dupOverlapCount <- ddply(dupOverlapRes, ("query"), nrow)
    names(dupOverlapCount)[2] <- "dupTag"
  } else {
    dupOverlapCount <- NULL
  }
  
  if ( !is.null(structuralVariation$inv) && (is.data.frame(structuralVariation$inv)) && 
       (nrow(structuralVariation$inv)>0) ) {
    invDf <- structuralVariation$inv[(structuralVariation$inv)$chromosome %in% 
                                       genome$chr, ]
    invIrange <- GRanges(seqnames=invDf$chromosome, 
                         IRanges(start=as.numeric(invDf$pos1), 
                                 end=as.numeric(invDf$pos2)))
    invOverlap <- findOverlaps(genomeResIrange, invIrange)
    invOverlapRes <- as.data.frame(cbind(queryHits(invOverlap), 
                                         subjectHits(invOverlap)))
    names(invOverlapRes) <- c("query", "subject")
    invOverlapCount <- ddply(invOverlapRes, ("query"), nrow)
    names(invOverlapCount)[2] <- "invTag"
  } else {
    invOverlapCount <- NULL
  }
  
  if (!is.null(dupOverlapCount)) {
    genomeDelDup <- merge(genomeDel, dupOverlapCount, by="query", all=T)
  } else {
    genomeDelDup <- genomeDel
  }
  if (!is.null(invOverlapCount)) {
    genomeDelDupInv <- merge(genomeDelDup, invOverlapCount, by="query", all=T)
  } else {
    genomeDelDupInv <- genomeDelDup
  }
  
  genomeDelDupInvIrange <- GRanges(seqnames=genomeDelDupInv$chr, 
                                   IRanges(start=genomeDelDupInv$pos1, 
                                           end=genomeDelDupInv$pos2))
  ## number of deletions in each window
  if (!is.null(delOverlapCount)) {
    genomeDelDupInv$delTag[is.na(genomeDelDupInv$delTag)] <- 0
    genomeDelDupInvIrange$delScore <- genomeDelDupInv$delTag
  }
  ## number of duplications in each window
  if (!is.null(dupOverlapCount)) {
    genomeDelDupInv$dupTag[is.na(genomeDelDupInv$dupTag)] <- 0 
    genomeDelDupInvIrange$dupScore <- genomeDelDupInv$dupTag
  }
  ## number of inversions in each window
  if (!is.null(invOverlapCount)) {
    genomeDelDupInv$invTag[is.na(genomeDelDupInv$invTag)] <- 0
    genomeDelDupInvIrange$invScore <- genomeDelDupInv$invTag
  }
  
  seqlengths(genomeDelDupInvIrange) <- genome$end[
    match(names(seqlengths(genomeDelDupInvIrange)), genome$chr)]
  p <- ggbio() 
  
  ## barplot of inversions
  if (!is.null(invOverlapCount)) {
    p <- p + circle(genomeDelDupInvIrange, geom="bar", aes(y="invScore"), 
                    color="green4", fill="green4")
  }
  
  ## barplot of duplications
  if (!is.null(dupOverlapCount)) {
    p <- p + circle(genomeDelDupInvIrange, geom="bar", aes(y="dupScore"), 
                    color="red", fill="red")
  }
  
  ## barplot of deletions
  if (!is.null(delOverlapCount)) {
    p <- p + circle(genomeDelDupInvIrange, geom="bar", aes(y="delScore"), 
                    color="blue", fill="blue")
  }
  
  p <- p + circle(genomeDelDupInvIrange, geom="ideo", 
                                fill="gray70")
  p <- p + circle(genomeDelDupInvIrange,  geom = "scale",  
                         size = 2)
  ## chromosomes names
  p <- p + circle(genome.gr, geom = "text",  
                         aes(label = seqnames),  vjust = 0)  
  
  return(p)
}
