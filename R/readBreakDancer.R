
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
                            regSizeLowerCutoff=100, regSizeUpperCutoff=10000000,
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

    bdPred_clean <- bdPred[, c(1, 2, 4, 5:8)]
    names(bdPred_clean) <- c("chromosome", "start", "end", "type","size", "score", "rp_support")

    bdDel <- mergeOLCNVs(bdPred_clean[bdPred_clean$type=="DEL", ],software= method)
    bdInv <- mergeOLCNVs(bdPred_clean[bdPred_clean$type=="INV", ],software= method)

    retuRes <- list(del=bdDel, inv=bdInv)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes);
}
