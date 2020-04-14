
## merging overlapped SVs predicted by DELLY
DellyCluster <- function(df) 
{
    maxReadPairSupp <- max(df$rp_support)
    dfFil <- df[df$rp_support>=(maxReadPairSupp/2), ]
    dfFilIrange <- IRanges(start=dfFil$start, end=dfFil$end)
    outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
    dfFil$cluDelly <- subjectHits(outTmp)
    dfFilRes <- ddply(dfFil, ("cluDelly"), function(x){
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
readDelly <- function(file="", regSizeLowerCutoff=100, 
                      regSizeUpperCutoff=1000000, readsSupport=3,
                      method="Delly", ...) 
{
    ## reading in SV predictions
    DellyCont <- read.table(file, as.is=T, ...)
    DellyCont <- DellyCont[, c(1:2, 8, 10)]
    names(DellyCont) <- c("chr", "start", "info", "detail")
    
    DellyCont$type <- gsub("SVTYPE=([A-Z]+);.+", "\\1", DellyCont$info)
    DellyCont$type <- gsub(".+;", "\\1", DellyCont$type)
    DellyCont <- DellyCont[DellyCont$type %in% c("DEL", "DUP", "INV"), ]
    DellyCont$end <- as.numeric(gsub(".+;END=(\\d+);.+", "\\1", DellyCont$info))
    DellyCont$pe_support <- 0
    DellyCont$pe_support[grepl(".+;PE=(\\d+).+", DellyCont$info)] <- 
      as.numeric(gsub(".+;PE=(\\d+).+", "\\1", DellyCont$info[grepl(".+;PE=(\\d+).+", DellyCont$info)]))
    DellyCont$su_support <- 0
    DellyCont$su_support[grepl(".+;SU=(\\d+).+", DellyCont$info)] <- 
      as.numeric(gsub(".+;SU=(\\d+).+", "\\1", DellyCont$info[grepl(".+;SU=(\\d+).+", DellyCont$info)]))
    DellyCont$rp_support <- DellyCont$pe_support + DellyCont$su_support
    DellyCont$pe_support <- NULL
    DellyCont$su_support <- NULL
    DellyCont$size <- as.numeric(DellyCont$end - DellyCont$start)
    
    DellyCont <- DellyCont[DellyCont$rp_support>=readsSupport&
                            DellyCont$size>=regSizeLowerCutoff&
                            DellyCont$size<=regSizeUpperCutoff, ]
    
    DellyCont <- DellyCont[, c("chr", "start", "end", "size", "rp_support",
                               "type")]

    DellyDelDf <- DellyCont[DellyCont$type=="DEL", ]

    ## filtering and merging deletions
    if (!is.null(DellyDelDf) && nrow(DellyDelDf)>0) {
        DellyDelIrange <- GRanges(seqnames=DellyDelDf$chr, 
                              ranges=IRanges(start=DellyDelDf$start, 
                                             end=DellyDelDf$end))
        DellyDelIrangeRes <- findOverlaps(DellyDelIrange, reduce(DellyDelIrange))
        DellyDelDf$clu <- subjectHits(DellyDelIrangeRes)
        DellyDelDfFilMer <- ddply(DellyDelDf, ("clu"), DellyCluster)
        DellyDelDfFilMer <- DellyDelDfFilMer[, c("chr", "start", "end", "size", "rp_support")]
        names(DellyDelDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport")
	DellyDelDfFilMer$info <- paste0("SU=", DellyDelDfFilMer$readsSupport)
	DellyDelDfFilMer$readsSupport <- NULL
    } else {
        DellyDelDfFilMer <- NULL
    }

    ## reading duplications
    DellyDupDf <- DellyCont[DellyCont$type=="DUP", ]

    ## filtering and merging duplications
    if (!is.null(DellyDupDf) && nrow(DellyDupDf)>0) {
        DellyDupIrange <- GRanges(seqnames=DellyDupDf$chr, 
                              ranges=IRanges(start=DellyDupDf$start, 
                                             end=DellyDupDf$end))
        DellyDupIrangeRes <- findOverlaps(DellyDupIrange, reduce(DellyDupIrange))
        DellyDupDf$clu <- subjectHits(DellyDupIrangeRes)
        DellyDupDfFilMer <- ddply(DellyDupDf, ("clu"), DellyCluster)
        DellyDupDfFilMer <- DellyDupDfFilMer[, c("chr", "start", "end", "size", "rp_support")]
        names(DellyDupDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport")
	DellyDupDfFilMer$info <- paste0("SU=", DellyDupDfFilMer$readsSupport)
	DellyDupDfFilMer$readsSupport <- NULL
    } else {
        DellyDupDfFilMer <- NULL
    }

    ## reading inversions
    DellyInvDf <- DellyCont[DellyCont$type=="INV", ]
    
    ## filtering and merging inversions
    if (!is.null(DellyInvDf) && nrow(DellyInvDf)>0) {
        DellyInvIrange <- GRanges(seqnames=DellyInvDf$chr, 
                              ranges=IRanges(start=DellyInvDf$start, 
                                             end=DellyInvDf$end))
        DellyInvIrangeRes <- findOverlaps(DellyInvIrange, reduce(DellyInvIrange))
        DellyInvDf$clu <- subjectHits(DellyInvIrangeRes)
        DellyInvDfFilMer <- ddply(DellyInvDf, ("clu"), DellyCluster)
        DellyInvDfFilMer <- DellyInvDfFilMer[, c("chr", "start", "end", "size", "rp_support")]
        names(DellyInvDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport")
	DellyInvDfFilMer$info <- paste0("SU=", DellyInvDfFilMer$readsSupport)
	DellyInvDfFilMer$readsSupport <- NULL
    } else {
        DellyInvDfFilMer <- NULL
    }

    retuRes <- list(del=DellyDelDfFilMer, dup=DellyDupDfFilMer, 
                inv=DellyInvDfFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes)
}
