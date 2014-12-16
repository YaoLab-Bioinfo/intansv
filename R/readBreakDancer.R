
## Merging overlapped SVs predicted by breakDancer
breakDancerCluster <- function(df)
{
    maxScore <- max(df$score)
    maxScoreDf <- df[df$score==maxScore, ]
    maxReadPairSupp <- min(maxScoreDf$ReadPairSupp)
    dfFil <- df[df$score>=(maxScore-20)&df$ReadPairSupp>=(maxReadPairSupp/2), ]
    dfFilIrange <- IRanges(start=dfFil$pos1, end=dfFil$pos2)
    outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
    dfFil$clubd <- subjectHits(outTmp)
    dfFilRes <- ddply(dfFil, ("clubd"), function(x){
        if(nrow(x)==1){
            return(x)
        } else {
            leftMin <- min(x$pos1)
            rightMax <- max(x$pos2)
            rangeLength <- rightMax-leftMin
            x$op <- (x$pos2-x$pos1)/rangeLength
            if(any(x$op<0.8)){
                return(NULL)
            } else {
                return(x[which.max(x$ReadPairSupp), ])
            }
        }
    })
}

## Reading in the predicted SVs given by breakDancer
readBreakDancer <- function(file="", scoreCutoff=60, readsSupport=3, 
                            regSizeLowerCutoff=100, regSizeUpperCutoff=1000000,
                            method="BreakDancer", ...)
{
    bdColClass <- c("character", "numeric", "NULL", "character", "numeric", 
                    "NULL", "character", "numeric", "numeric", "numeric")
    bdPred <- read.table(file, colClasses=bdColClass, fill=T, ...)
    bdPred <- bdPred[,1:8]
    names(bdPred) <- c("chr1", "pos1", "chr2", "pos2", "type", "size", 
                       "score", "ReadPairSupp")
    bdPred <- bdPred[bdPred$score>=scoreCutoff&
                     bdPred$ReadPairSupp>=readsSupport&
                     abs(bdPred$size)>=regSizeLowerCutoff&
                     abs(bdPred$size)<=regSizeUpperCutoff, ]

    ## filtering and merging deletions
    bdDel <- bdPred[bdPred$type=="DEL", ]
    bdDel$mid <- (bdDel$pos2+bdDel$pos1)/2
    bdDel$pos1 <- round(bdDel$mid-bdDel$size/2)
    bdDel$pos2 <- round(bdDel$mid+bdDel$size/2)
    bdDelIrange <- GRanges(seqnames=bdDel$chr1, 
                           ranges=IRanges(start=bdDel$pos1, end=bdDel$pos2))
    bdDelIrangeRes <- findOverlaps(bdDelIrange, reduce(bdDelIrange))
    bdDel$clu <- subjectHits(bdDelIrangeRes)
    bdDelFilMer <- ddply(bdDel, ("clu"), breakDancerCluster)
    if (nrow(bdDelFilMer)==0) {
        bdDelFilMer <- NULL
    } else {
        bdDelFilMer <- bdDelFilMer[, c(1, 2, 4, 6)]
        names(bdDelFilMer) <- c("chromosome", "pos1", "pos2", "size")
    }

    ## filtering and merging inversions
    bdInv <- bdPred[bdPred$type=="INV", ]
    if (nrow(bdInv)==0) {
        bdInvFilMer <- NULL
    } else {
        bdInv$size <- abs(bdInv$size)
        bdInv$mid <- (bdInv$pos1 + bdInv$pos2)/2
        bdInv$pos1 <- round(bdInv$mid - bdInv$size/2)
        bdInv$pos2 <- round(bdInv$mid + bdInv$size/2)
        bdInvIrange <- GRanges(seqnames=bdInv$chr1, 
                           ranges=IRanges(start=bdInv$pos1, 
                           end=bdInv$pos2))
        bdInvIrangeRes <- findOverlaps(bdInvIrange, reduce(bdInvIrange))
        bdInv$clu <- subjectHits(bdInvIrangeRes)
        bdInvFilMer <- ddply(bdInv, ("clu"), breakDancerCluster)
        if (nrow(bdInvFilMer)==0) {
            bdInvFilMer <- NULL
        } else {
            bdInvFilMer <- bdInvFilMer[, c(1, 2, 4, 6)]
            names(bdInvFilMer) <- c("chromosome", "pos1", "pos2", "size")
        }
    }

    retuRes <- list(del=bdDelFilMer, inv=bdInvFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes);
}
