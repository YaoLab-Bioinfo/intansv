
## merging overlapped SVs predicted by DELLY
LumpyCluster <- function(df) 
{
    maxReadPairSupp <- max(df$rp_support)
    dfFil <- df[df$rp_support>=(maxReadPairSupp/2), ]
    dfFilIrange <- IRanges(start=dfFil$start, end=dfFil$end)
    outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
    dfFil$cluLumpy <- subjectHits(outTmp)
    dfFilRes <- ddply(dfFil, ("cluLumpy"), function(x){
        if(nrow(x)==1){
            return(x)
        } else {
            LeftMin <- min(x$start)
            RightMax <- max(x$end)
            RangeLength <- RightMax-LeftMin
            x$op <- (x$end - x$start)/RangeLength
            if(any(x$op<0.8)){
                return(NULL)
            } else {
                return(x[which.max(x$rp_support), ])
            }
        }
    })
}

## Reading in the predicted SVs given by DELLY
readLumpy <- function(file="", regSizeLowerCutoff=100, 
                      regSizeUpperCutoff=1000000, readsSupport=3,
                      method="Lumpy", ...) 
{
    ## reading in SV predictions
    LumpyCont <- read.table(file, as.is=T, ...)
    LumpyCont <- LumpyCont[, c(1:2, 8, 10)]
    names(LumpyCont) <- c("chr", "start", "info", "detail")
    
    LumpyCont$type <- gsub("SVTYPE=([A-Z]+);.+", "\\1", LumpyCont$info)
    LumpyCont <- LumpyCont[LumpyCont$type %in% c("DEL", "DUP", "INV"), ]
    LumpyCont$end <- as.numeric(gsub(".+;END=(\\d+);.+", "\\1", LumpyCont$info))
    LumpyCont$rp_support <- as.numeric(gsub(".+;SU=(\\d+).+", "\\1", LumpyCont$info))
    LumpyCont$size <- as.numeric(LumpyCont$end - LumpyCont$start)
    
    LumpyCont <- LumpyCont[LumpyCont$rp_support>=readsSupport&
                            LumpyCont$size>=regSizeLowerCutoff&
                            LumpyCont$size<=regSizeUpperCutoff, ]
    
    LumpyCont <- LumpyCont[, c("chr", "start", "end", "size", "rp_support",
                               "type")]

    LumpyDelDf <- LumpyCont[LumpyCont$type=="DEL", ]

    ## filtering and merging deletions
    if (!is.null(LumpyDelDf) && nrow(LumpyDelDf)>0) {
        LumpyDelIrange <- GRanges(seqnames=LumpyDelDf$chr, 
                              ranges=IRanges(start=LumpyDelDf$start, 
                                             end=LumpyDelDf$end))
        LumpyDelIrangeRes <- findOverlaps(LumpyDelIrange, reduce(LumpyDelIrange))
        LumpyDelDf$clu <- subjectHits(LumpyDelIrangeRes)
        LumpyDelDfFilMer <- ddply(LumpyDelDf, ("clu"), LumpyCluster)
        LumpyDelDfFilMer <- LumpyDelDfFilMer[, c("chr", "start", "end", "size", "rp_support")]
        names(LumpyDelDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport")
	LumpyDelDfFilMer$info <- paste0("SU=", LumpyDelDfFilMer$readsSupport)
	LumpyDelDfFilMer$readsSupport <- NULL
    } else {
        LumpyDelDfFilMer <- NULL
    }

    ## reading duplications
    LumpyDupDf <- LumpyCont[LumpyCont$type=="DUP", ]

    ## filtering and merging duplications
    if (!is.null(LumpyDupDf) && nrow(LumpyDupDf)>0) {
        LumpyDupIrange <- GRanges(seqnames=LumpyDupDf$chr, 
                              ranges=IRanges(start=LumpyDupDf$start, 
                                             end=LumpyDupDf$end))
        LumpyDupIrangeRes <- findOverlaps(LumpyDupIrange, reduce(LumpyDupIrange))
        LumpyDupDf$clu <- subjectHits(LumpyDupIrangeRes)
        LumpyDupDfFilMer <- ddply(LumpyDupDf, ("clu"), LumpyCluster)
        LumpyDupDfFilMer <- LumpyDupDfFilMer[, c("chr", "start", "end", "size", "rp_support")]
        names(LumpyDupDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport")
	LumpyDupDfFilMer$info <- paste0("SU=", LumpyDupDfFilMer$readsSupport)
	LumpyDupDfFilMer$readsSupport <- NULL
    } else {
        LumpyDupDfFilMer <- NULL
    }

    ## reading inversions
    LumpyInvDf <- LumpyCont[LumpyCont$type=="INV", ]
    
    ## filtering and merging inversions
    if (!is.null(LumpyInvDf) && nrow(LumpyInvDf)>0) {
        LumpyInvIrange <- GRanges(seqnames=LumpyInvDf$chr, 
                              ranges=IRanges(start=LumpyInvDf$start, 
                                             end=LumpyInvDf$end))
        LumpyInvIrangeRes <- findOverlaps(LumpyInvIrange, reduce(LumpyInvIrange))
        LumpyInvDf$clu <- subjectHits(LumpyInvIrangeRes)
        LumpyInvDfFilMer <- ddply(LumpyInvDf, ("clu"), LumpyCluster)
        LumpyInvDfFilMer <- LumpyInvDfFilMer[, c("chr", "start", "end", "size", "rp_support")]
        names(LumpyInvDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport")
	LumpyInvDfFilMer$info <- paste0("SU=", LumpyInvDfFilMer$readsSupport)
	LumpyInvDfFilMer$readsSupport <- NULL
    } else {
        LumpyInvDfFilMer <- NULL
    }

    retuRes <- list(del=LumpyDelDfFilMer, dup=LumpyDupDfFilMer, 
                inv=LumpyInvDfFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes)
}
