
## Reading in the predicted SVs given by SVseq2
readSvseq <- function(dataDir=".", regSizeLowerCutoff=100, method="SVseq2",
                      regSizeUpperCutoff=1000000, readsSupport=3) 
{
    SvseqDelList <- list.files(dataDir, full.names=T, pattern=".+\\.del$")

    SvseqDel <- lapply(SvseqDelList, function(x){
      SvseqData <- try(read.table(x, comment.char=")", fill=T, as.is=T,
                                  col.names=paste("V", 1:8, sep=""),),silent=T)
        if (is.data.frame(SvseqData)) {
            brk <- which(grepl("^#", SvseqData$V1))
            if (length(brk)<2) {
                return(NULL)
            } else {
              SvseqData <- SvseqData[brk[1]:nrow(SvseqData), ]
              brk <- which(grepl("^#", SvseqData$V1))
              brkTmp <- brk[2:length(brk)] - brk[1:(length(brk)-1)]
              SvseqData <- SvseqData[1:(nrow(SvseqData)-1), ]
              brkRes <- rep(1:length(brkTmp), brkTmp)
              SvseqDataList <- split(SvseqData, brkRes)
    
                SvseqDataListDf <- lapply(SvseqDataList, function(df){
                    dfRes <- df[2, ]
                    dfRes$readsSupport <- nrow(df)-3
                    dfRes$V1 <- NULL
                    return(dfRes)
                })
    
                SvseqDataDf <- do.call(rbind, SvseqDataListDf)
                return(SvseqDataDf);
            }
        } else {
            return(NULL)
        }
    })

    SvseqDelDf <- do.call(rbind, SvseqDel)
    SvseqDelDf$left <- round((as.numeric(SvseqDelDf$V3) + 
                              as.numeric(SvseqDelDf$V4))/2)
    SvseqDelDf$right <- round((as.numeric(SvseqDelDf$V5) + 
                               as.numeric(SvseqDelDf$V6))/2)
    SvseqDelDf$size <- SvseqDelDf$right - SvseqDelDf$left 
    names(SvseqDelDf)[1] <- "chr"

    ## filtering and merging deletions
    SvseqDelDf <- SvseqDelDf[SvseqDelDf$size>=regSizeLowerCutoff &
                             SvseqDelDf$size<=regSizeUpperCutoff &
                             SvseqDelDf$readsSupport>=readsSupport, ]
    SvseqDelIrange <- GRanges(seqnames=SvseqDelDf$chr, 
                              ranges=IRanges(start=SvseqDelDf$left, 
                                             end=SvseqDelDf$right))
    SvseqDelIrangeRes <- findOverlaps(SvseqDelIrange, reduce(SvseqDelIrange))
    SvseqDelDf$clu <- subjectHits(SvseqDelIrangeRes)
    SvseqDelDfFilMer <- ddply(SvseqDelDf, ("clu"), function(df){
        maxReadPairSupp <- max(df$readsSupport)
        dfFil <- df[df$readsSupport>=(maxReadPairSupp/2), ]
        dfFilIrange <- IRanges(start=dfFil$left, end=dfFil$right)
        outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
        dfFil$clusvseq <- subjectHits(outTmp)
        dfFilRes <- ddply(dfFil, ("clusvseq"), function(x){
            if(nrow(x)==1){
                return(x)
            } else {
                LeftMin <- min(x$left)
                RightMax <- max(x$right)
                RangeLength <- RightMax-LeftMin
                x$op <- (x$right-x$left)/RangeLength
                if(any(x$op<0.8)){
                    return(NULL)
                } else {
                    return(x[which.max(x$readsSupport), ])
                }
            }
        })
    })
    
    SvseqDelDfFilMer <- SvseqDelDfFilMer[, c(1, 9:11)]
    names(SvseqDelDfFilMer) <- c("chromosome", "pos1", "pos2", "size")
    retuRes <- list(del=SvseqDelDfFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes);
}

