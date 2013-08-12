
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
    delDf <- structuralVariation$del[(structuralVariation$del)$chromosome%in%
                                      seqnames(genomeAnnotation)@values, ]
    dupDf <- structuralVariation$dup[(structuralVariation$dup)$chromosome%in%
                                      seqnames(genomeAnnotation)@values, ]
    invDf <- structuralVariation$inv[(structuralVariation$inv)$chromosome%in%
                                      seqnames(genomeAnnotation)@values, ]
    delIrange <- GRanges(seqnames=delDf$chromosome, 
                         IRanges(start=as.numeric(delDf$pos1), 
                         end=as.numeric(delDf$pos2)))
    dupIrange <- GRanges(seqnames=dupDf$chromosome, 
                         IRanges(start=as.numeric(dupDf$pos1), 
                         end=as.numeric(dupDf$pos2)))
    invIrange <- GRanges(seqnames=invDf$chromosome, 
                         IRanges(start=as.numeric(invDf$pos1), 
                         end=as.numeric(invDf$pos2)))

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
    delOverlap <- findOverlaps(genomeDfResIrange, delIrange)
    dupOverlap <- findOverlaps(genomeDfResIrange, dupIrange)
    invOverlap <- findOverlaps(genomeDfResIrange, invIrange)
    delOverlapRes <- as.data.frame(cbind(queryHits(delOverlap), 
                                         subjectHits(delOverlap)))
    names(delOverlapRes) <- c("query", "delTag")
    genomeDfRes$query <- 1:nrow(genomeDfRes)
    delOverlapCount <- ddply(delOverlapRes, ("query"), nrow)
    names(delOverlapCount)[2] <- "delTag"
    genomeDfDel <- merge(genomeDfRes, delOverlapCount, by="query", all=T)
    dupOverlapRes <- as.data.frame(cbind(queryHits(dupOverlap), 
                                         subjectHits(dupOverlap)))
    names(dupOverlapRes) <- c("query", "subject")
    dupOverlapCount <- ddply(dupOverlapRes, ("query"), nrow)
    names(dupOverlapCount)[2] <- "dupTag"
    genomeDfDelDup <- merge(genomeDfDel, dupOverlapCount, by="query", all=T)
    invOverlapRes <- as.data.frame(cbind(queryHits(invOverlap), 
                                         subjectHits(invOverlap)))
    names(invOverlapRes) <- c("query", "subject")
    invOverlapCount <- ddply(invOverlapRes, ("query"), nrow)
    names(invOverlapCount)[2] <- "invTag"
    genomeDfDelDupInv <- merge(genomeDfDelDup, invOverlapCount, by="query", all=T)
    ## number of deletions in each window
    genomeDfDelDupInv$delTag[is.na(genomeDfDelDupInv$delTag)] <- 0   
    ## number of duplications in each window
    genomeDfDelDupInv$dupTag[is.na(genomeDfDelDupInv$dupTag)] <- 0   
    ## number of inversions in each window
    genomeDfDelDupInv$invTag[is.na(genomeDfDelDupInv$invTag)] <- 0   
    genomeDfDelDupInvIrange <- GRanges(seqnames=genomeDfDelDupInv$V1, 
                                       IRanges(start=genomeDfDelDupInv$pos1, 
                                               end=genomeDfDelDupInv$pos2), 
                                       delScore=genomeDfDelDupInv$delTag, 
                                       dupScore=genomeDfDelDupInv$dupTag, 
                                       invScore=genomeDfDelDupInv$invTag)
    seqlengths(genomeDfDelDupInvIrange) <- seqlengths(genomeAnnotation)
    p <- ggplot() + layout_circle(genomeDfDelDupInvIrange, geom="ideo", 
                                  fill="gray70", radius=30, trackwidth=4)
    p <- p + layout_circle(genomeDfDelDupInvIrange,  geom = "scale",  
                           size = 2,  radius = 35,  trackWidth = 2)
    ## chromosomes names
    p <- p + layout_circle(genomeAnnotation, geom = "text",  
                           aes(label = seqnames),  vjust = 0, 
                           radius = 38,  trackWidth = 4)  
    ## barplot of deletions
    p <- p + layout_circle(genomeDfDelDupInvIrange, geom="bar", radius=25, 
                           trackwidth=4, aes(y="delScore"), 
                           color="blue", fill="blue")
    ## barplot of duplications
    p <- p + layout_circle(genomeDfDelDupInvIrange, geom="bar", radius=20, 
                           trackwidth=4, aes(y="dupScore"), 
                           color="red", fill="red")
    ## barplot of inversions
    p <- p + layout_circle(genomeDfDelDupInvIrange, geom="bar", radius=15, 
                           trackwidth=4, aes(y="invScore"), 
                           color="green4", fill="green4")
    return(p)
}

