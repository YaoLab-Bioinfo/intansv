
## Merging overlapped SVs predicted by Lumpy
LumpyCluster <- function(df)
{
  maxReadPairSupp <- max(df$ReadPairSupp)
  dfFil <- df[df$ReadPairSupp>=(maxReadPairSupp/2), ]
  dfFilIrange <- IRanges(start=dfFil$pos1, end=dfFil$pos2)
  outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
  dfFil$clulumpy <- subjectHits(outTmp)
  dfFilRes <- ddply(dfFil, ("clulumpy"), function(x){
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


## Reading in the predicted SVs given by Lumpy
readLumpy <- function(file="", regSizeLowerCutoff=100, readsSupport=3, 
                           method="Lumpy", regSizeUpperCutoff=1000000, 
                      breakpointThres=200, scoreCut=0.1, ...)
{
  LumpyColClass <- c("character", "numeric", "numeric", "character", "numeric", 
                     "numeric", "NULL", "numeric", "NULL", "NULL", "character",
                     "NULL", "character", "NULL", "NULL")
  LumpyPred <- read.table(file, colClasses=LumpyColClass, as.is=T, ...)
  names(LumpyPred) <- c("chr1", "start1", "end1", "chr2", "start2", "end2",
                        "score",  "type", "ReadPairSupp")
  LumpyPred$type <- gsub("TYPE:", "", LumpyPred$type)
  LumpyPred$ReadPairSupp <- as.numeric(gsub(".*,", "", 
                                            LumpyPred$ReadPairSupp))
  LumpyPred <- LumpyPred[LumpyPred$end1<(LumpyPred$start2-100),]
  LumpyPred <- LumpyPred[(LumpyPred$end1-LumpyPred$start1)<=breakpointThres&
                         (LumpyPred$end2-LumpyPred$start2)<=breakpointThres, ]
  LumpyPred <- LumpyPred[LumpyPred$score<=scoreCut, ]
  LumpyPred <- LumpyPred[LumpyPred$type%in%c("DELETION", "DUPLICATION", 
                                             "INVERSION"), ]
  LumpyPred$pos1 <- round((LumpyPred$start1+LumpyPred$end1)/2)
  LumpyPred$pos2 <- round((LumpyPred$start2+LumpyPred$end2)/2)
  LumpyPred$size <- LumpyPred$pos2-LumpyPred$pos1+1
  LumpyPred$chr2 <- NULL;  
  LumpyPred <- LumpyPred[abs(LumpyPred$size)>=regSizeLowerCutoff&
                         abs(LumpyPred$size)<=regSizeUpperCutoff&
                         abs(LumpyPred$ReadPairSupp)>=readsSupport, ]
  LumpyPred <- LumpyPred[, c("chr1", "pos1", "pos2", "size", "type", "ReadPairSupp")]
  names(LumpyPred)[1] <- "chromosome"
  
  ## filtering and merging deletions
  LumpyDel <- LumpyPred[LumpyPred$type=="DELETION", ]
  
  if (nrow(LumpyDel)==0) {
    LumpyDelFilMer <- NULL
  } else {
    LumpyDel$size <- abs(LumpyDel$size)
    LumpyDel$mid <- (LumpyDel$pos1 + LumpyDel$pos2)/2
    LumpyDel$pos1 <- round(LumpyDel$mid - LumpyDel$size/2)
    LumpyDel$pos2 <- round(LumpyDel$mid + LumpyDel$size/2)
    LumpyDelIrange <- GRanges(seqnames=LumpyDel$chromosome, 
                                   ranges=IRanges(start=LumpyDel$pos1, 
                                                  end=LumpyDel$pos2))
    LumpyDelIrangeRes <- findOverlaps(LumpyDelIrange, reduce(LumpyDelIrange))
    LumpyDel$clu <- subjectHits(LumpyDelIrangeRes)
    LumpyDelFilMer <- ddply(LumpyDel, ("clu"), LumpyCluster)
    if (nrow(LumpyDelFilMer)==0) {
      LumpyDelFilMer <- NULL
    } else {
      LumpyDelFilMer <- LumpyDelFilMer[, c("chromosome", "pos1", "pos2", "size")]
    }
  }
  
  ## filtering and merging Inversions
  LumpyInv <- LumpyPred[LumpyPred$type=="INVERSION", ]
  
  if (nrow(LumpyInv)==0) {
    LumpyInvFilMer <- NULL
  } else {
    LumpyInv$size <- abs(LumpyInv$size)
    LumpyInv$mid <- (LumpyInv$pos1 + LumpyInv$pos2)/2
    LumpyInv$pos1 <- round(LumpyInv$mid - LumpyInv$size/2)
    LumpyInv$pos2 <- round(LumpyInv$mid + LumpyInv$size/2)
    LumpyInvIrange <- GRanges(seqnames=LumpyInv$chromosome, 
                                   ranges=IRanges(start=LumpyInv$pos1, 
                                                  end=LumpyInv$pos2))
    LumpyInvIrangeRes <- findOverlaps(LumpyInvIrange, reduce(LumpyInvIrange))
    LumpyInv$clu <- subjectHits(LumpyInvIrangeRes)
    LumpyInvFilMer <- ddply(LumpyInv, ("clu"), LumpyCluster)
    if (nrow(LumpyInvFilMer)==0) {
      LumpyInvFilMer <- NULL
    } else {
      LumpyInvFilMer <- LumpyInvFilMer[, c("chromosome", "pos1", "pos2", "size")]
    }
  }
  
  ## filtering and merging Duplications
  LumpyDup <- LumpyPred[LumpyPred$type=="DUPLICATION", ]
  
  if (nrow(LumpyDup)==0) {
    LumpyDupFilMer <- NULL
  } else {
    LumpyDup$size <- abs(LumpyDup$size)
    LumpyDup$mid <- (LumpyDup$pos1 + LumpyDup$pos2)/2
    LumpyDup$pos1 <- round(LumpyDup$mid - LumpyDup$size/2)
    LumpyDup$pos2 <- round(LumpyDup$mid + LumpyDup$size/2)
    LumpyDupIrange <- GRanges(seqnames=LumpyDup$chromosome, 
                                   ranges=IRanges(start=LumpyDup$pos1, 
                                                  end=LumpyDup$pos2))
    LumpyDupIrangeRes <- findOverlaps(LumpyDupIrange, reduce(LumpyDupIrange))
    LumpyDup$clu <- subjectHits(LumpyDupIrangeRes)
    LumpyDupFilMer <- ddply(LumpyDup, ("clu"), LumpyCluster)
    if (nrow(LumpyDupFilMer)==0) {
      LumpyDupFilMer <- NULL
    } else {
      LumpyDupFilMer <- LumpyDupFilMer[, c("chromosome", "pos1", "pos2", "size")]
    }
  }
  
  retuRes <- list(del=LumpyDelFilMer, dup=LumpyDupFilMer, 
                       inv=LumpyInvFilMer)
  attributes(retuRes) <- c(attributes(retuRes), list(method=method))
  
  return(retuRes);
}



