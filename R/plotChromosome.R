
## Displaying the distribution of SVs in the whole genome
plotChromosome <- function(genomeAnnotation, structuralVariation, 
                           windowSize=1000000) 
{
    genomeDf <- 
        as.data.frame(cbind(
                           rep(as.character(seqnames(genomeAnnotation)@values),
                                            seqnames(genomeAnnotation)@lengths), 
                           start(ranges(genomeAnnotation)), 
                           end(ranges(genomeAnnotation))), 
                           stringsAsFactors=F)
    genomeDf$V2 <- as.numeric(genomeDf$V2)
    genomeDf$V3 <- as.numeric(genomeDf$V3)
    genomeDfRes <- ddply(genomeDf, c("V1", "V3"), function(df){
      pos1 <- seq(1, df$V3[1], by=as.numeric(windowSize))
      pos2 <- seq(as.numeric(windowSize), df$V3[1], by=as.numeric(windowSize))
      if(length(pos1)!=length(pos2)){
        pos2<-c(pos2, df$V3[1])
      }
      dfRes <- as.data.frame(cbind(pos1, pos2))
      return(dfRes)
    })
    genomeDfResIrange <- GRanges(seqnames=genomeDfRes$V1, 
                                 IRanges(start=genomeDfRes$pos1, 
                                         end=genomeDfRes$pos2))
    genomeDfRes$query <- 1:nrow(genomeDfRes)
    
    if ( !is.null(structuralVariation$del) && (is.data.frame(structuralVariation$del)) &&
           (nrow(structuralVariation$del)>0) ) {
        delDf <- structuralVariation$del[(structuralVariation$del)$chromosome%in%
                                      seqnames(genomeAnnotation)@values, ]
        delIrange <- GRanges(seqnames=delDf$chromosome, 
                         IRanges(start=as.numeric(delDf$pos1), 
                                 end=as.numeric(delDf$pos2)))
        delOverlap <- findOverlaps(genomeDfResIrange, delIrange)
        delOverlapRes <- as.data.frame(cbind(queryHits(delOverlap), 
                                         subjectHits(delOverlap)))
        names(delOverlapRes) <- c("query", "delTag")
        delOverlapCount <- ddply(delOverlapRes, ("query"), nrow)
        names(delOverlapCount)[2] <- "delTag"
        genomeDfDel <- merge(genomeDfRes, delOverlapCount, by="query", all=T)
    } else {
        genomeDfDel <- genomeDfRes
        delOverlapCount <- NULL
    }
        
    if ( !is.null(structuralVariation$dup) && (is.data.frame(structuralVariation$dup)) && 
           (nrow(structuralVariation$dup)>0) ) {
        dupDf <- structuralVariation$dup[(structuralVariation$dup)$chromosome%in%
                                      seqnames(genomeAnnotation)@values, ]
        dupIrange <- GRanges(seqnames=dupDf$chromosome, 
                         IRanges(start=as.numeric(dupDf$pos1), 
                                 end=as.numeric(dupDf$pos2)))
        dupOverlap <- findOverlaps(genomeDfResIrange, dupIrange)    
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
        invDf <- structuralVariation$inv[(structuralVariation$inv)$chromosome%in%
                                      seqnames(genomeAnnotation)@values, ]
        invIrange <- GRanges(seqnames=invDf$chromosome, 
                         IRanges(start=as.numeric(invDf$pos1), 
                         end=as.numeric(invDf$pos2)))
        invOverlap <- findOverlaps(genomeDfResIrange, invIrange)
        invOverlapRes <- as.data.frame(cbind(queryHits(invOverlap), 
                                         subjectHits(invOverlap)))
        names(invOverlapRes) <- c("query", "subject")
        invOverlapCount <- ddply(invOverlapRes, ("query"), nrow)
        names(invOverlapCount)[2] <- "invTag"
    } else {
      invOverlapCount <- NULL
    }
    
    if (!is.null(dupOverlapCount)) {
        genomeDfDelDup <- merge(genomeDfDel, dupOverlapCount, by="query", all=T)
    } else {
        genomeDfDelDup <- genomeDfDel
    }
    if (!is.null(invOverlapCount)) {
        genomeDfDelDupInv <- merge(genomeDfDelDup, invOverlapCount, by="query", all=T)
    } else {
        genomeDfDelDupInv <- genomeDfDelDup
    }
    
    genomeDfDelDupInvIrange <- GRanges(seqnames=genomeDfDelDupInv$V1, 
                                       IRanges(start=genomeDfDelDupInv$pos1, 
                                               end=genomeDfDelDupInv$pos2))
    ## number of deletions in each window
    if (!is.null(delOverlapCount)) {
        genomeDfDelDupInv$delTag[is.na(genomeDfDelDupInv$delTag)] <- 0
        genomeDfDelDupInvIrange$delScore <- genomeDfDelDupInv$delTag
    }
    ## number of duplications in each window
    if (!is.null(dupOverlapCount)) {
        genomeDfDelDupInv$dupTag[is.na(genomeDfDelDupInv$dupTag)] <- 0 
        genomeDfDelDupInvIrange$dupScore <- genomeDfDelDupInv$dupTag
    }
    ## number of inversions in each window
    if (!is.null(invOverlapCount)) {
        genomeDfDelDupInv$invTag[is.na(genomeDfDelDupInv$invTag)] <- 0
        genomeDfDelDupInvIrange$invScore <- genomeDfDelDupInv$invTag
    }
    
    seqlengths(genomeDfDelDupInvIrange) <- seqlengths(genomeAnnotation)[
      match(names(seqlengths(genomeDfDelDupInvIrange)), names(seqlengths(genomeAnnotation)))]
    p <- ggplot() + layout_circle(genomeDfDelDupInvIrange, geom="ideo", 
                                  fill="gray70", radius=30, trackwidth=4)
    p <- p + layout_circle(genomeDfDelDupInvIrange,  geom = "scale",  
                           size = 2,  radius = 35,  trackWidth = 2)
    ## chromosomes names
    p <- p + layout_circle(genomeAnnotation, geom = "text",  
                           aes(label = seqnames),  vjust = 0, 
                           radius = 38,  trackWidth = 4)  
    ## barplot of deletions
    if (!is.null(delOverlapCount)) {
        p <- p + layout_circle(genomeDfDelDupInvIrange, geom="bar", radius=25, 
                           trackwidth=4, aes(y="delScore"), 
                           color="blue", fill="blue")
    }
    ## barplot of duplications
    if (!is.null(dupOverlapCount)) {
        p <- p + layout_circle(genomeDfDelDupInvIrange, geom="bar", radius=20, 
                           trackwidth=4, aes(y="dupScore"), 
                           color="red", fill="red")
    }
    ## barplot of inversions
    if (!is.null(invOverlapCount)) {
        p <- p + layout_circle(genomeDfDelDupInvIrange, geom="bar", radius=15, 
                           trackwidth=4, aes(y="invScore"), 
                           color="green4", fill="green4")
    }
    return(p)
}

