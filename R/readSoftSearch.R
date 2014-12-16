
## Merging overlapped SVs predicted by softSearch
softSearchCluster <- function(df)
{
  maxReadPairSupp <- max(df$ReadPairSupp)
  dfFil <- df[df$ReadPairSupp>=(maxReadPairSupp/2), ]
  dfFilIrange <- IRanges(start=dfFil$pos1, end=dfFil$pos2)
  outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
  dfFil$clusoft <- subjectHits(outTmp)
  dfFilRes <- ddply(dfFil, ("clusoft"), function(x){
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


## Reading in the predicted SVs given by softSearch
readSoftSearch <- function(file="", regSizeLowerCutoff=100, readsSupport=3, 
                           method="softSearch", regSizeUpperCutoff=1000000, 
                           softClipsSupport=3, ...) 
{
  softSearchColClass <- c("character", "numeric", "NULL", "NULL", "character",
                          "NULL", "character", "character", "character", "character")
  softSearchPred <- read.table(file, colClasses=softSearchColClass, as.is=T, ...)
  names(softSearchPred) <- c("chr1", "start", "alt", "filter", "info", 
                             "format", "detail")
  softSearchPred <- softSearchPred[softSearchPred$filter=="PASS", ]
  softSearchPred$alt <- gsub("]", "", softSearchPred$alt)
  softSearchPred$chr2 <- gsub(":.*", "", softSearchPred$alt)
  softSearchPred$alt <- NULL; softSearchPred$filter <- NULL;
  softSearchPred$type <- gsub(".*EVENT=([A-Z_]+);.*", "\\1", softSearchPred$info)
  softSearchPred$size <- as.numeric(gsub(".*ISIZE=([0-9]+);.*", "\\1", 
                                         softSearchPred$info))
  softSearchPred$end <- as.numeric(gsub(".*END=([0-9]+);.*", "\\1", 
                                         softSearchPred$info))
  softSearchPred$info <- NULL;
  
  softSearchPred <- softSearchPred[abs(softSearchPred$size)>=regSizeLowerCutoff&
                                     abs(softSearchPred$size)<=regSizeUpperCutoff, ]
  
  ## filtering and merging deletions
  softSearchDel <- softSearchPred[softSearchPred$type=="DEL", ]
  delTag <- which(unlist(strsplit(softSearchDel$format[1], ":"))=="DEL")
  softSearchDel$ReadPairSupp <- sapply(strsplit(softSearchDel$detail,":"), 
                                       function(x){x[delTag]})
  softTag <- which(unlist(strsplit(softSearchDel$format[1], ":"))=="nSC")
  softSearchDel$softSupp <- sapply(strsplit(softSearchDel$detail,":"), 
                                   function(x){x[softTag]})
  softSearchDel$softSupp[softSearchDel$softSupp=="NA"] <- 0
  softSearchDel$softSupp <- as.numeric(softSearchDel$softSupp)
  
  softSearchDel$ReadPairSupp <- as.numeric(softSearchDel$ReadPairSupp)
  softSearchDel <- softSearchDel[softSearchDel$ReadPairSupp>=readsSupport|
                                   softSearchDel$softSupp>=softClipsSupport,]
  softSearchDel$pos1 <- pmin(softSearchDel$start, softSearchDel$end)
  softSearchDel$pos2 <- pmax(softSearchDel$start, softSearchDel$end)
  softSearchDel <- softSearchDel[, c("chr1", "pos1", "pos2", "size", "ReadPairSupp")]
  names(softSearchDel)[1] <- "chromosome"
  
  if (nrow(softSearchDel)==0) {
    softSearchDelFilMer <- NULL
  } else {
    softSearchDel$size <- abs(softSearchDel$size)
    softSearchDel$mid <- (softSearchDel$pos1 + softSearchDel$pos2)/2
    softSearchDel$pos1 <- round(softSearchDel$mid - softSearchDel$size/2)
    softSearchDel$pos2 <- round(softSearchDel$mid + softSearchDel$size/2)
    softSearchDelIrange <- GRanges(seqnames=softSearchDel$chromosome, 
                                   ranges=IRanges(start=softSearchDel$pos1, 
                                                  end=softSearchDel$pos2))
    softSearchDelIrangeRes <- findOverlaps(softSearchDelIrange, reduce(softSearchDelIrange))
    softSearchDel$clu <- subjectHits(softSearchDelIrangeRes)
    softSearchDelFilMer <- ddply(softSearchDel, ("clu"), softSearchCluster)
    if (nrow(softSearchDelFilMer)==0) {
      softSearchDelFilMer <- NULL
    } else {
      softSearchDelFilMer <- softSearchDelFilMer[, c("chromosome", "pos1", "pos2", "size")]
    }
  }
  
  ## filtering and merging Inversions
  softSearchInv <- softSearchPred[softSearchPred$type=="INV", ]
  InvTag <- which(unlist(strsplit(softSearchInv$format[1], ":"))=="INV")
  softSearchInv$ReadPairSupp <- sapply(strsplit(softSearchInv$detail,":"), 
                                       function(x){x[InvTag]})
  softTag <- which(unlist(strsplit(softSearchInv$format[1], ":"))=="nSC")
  softSearchInv$softSupp <- sapply(strsplit(softSearchInv$detail,":"), 
                                   function(x){x[softTag]})
  softSearchInv$softSupp[softSearchInv$softSupp=="NA"] <- 0
  softSearchInv$softSupp <- as.numeric(softSearchInv$softSupp)
  softSearchInv$ReadPairSupp <- as.numeric(softSearchInv$ReadPairSupp)
  softSearchInv <- softSearchInv[softSearchInv$ReadPairSupp>=readsSupport|
                                   softSearchInv$softSupp>=softClipsSupport,]
  softSearchInv$pos1 <- pmin(softSearchInv$start, softSearchInv$end)
  softSearchInv$pos2 <- pmax(softSearchInv$start, softSearchInv$end)
  softSearchInv <- softSearchInv[, c("chr1", "pos1", "pos2", "size", "ReadPairSupp")]
  names(softSearchInv)[1] <- "chromosome"
  
  if (nrow(softSearchInv)==0) {
    softSearchInvFilMer <- NULL
  } else {
    softSearchInv$size <- abs(softSearchInv$size)
    softSearchInv$mid <- (softSearchInv$pos1 + softSearchInv$pos2)/2
    softSearchInv$pos1 <- round(softSearchInv$mid - softSearchInv$size/2)
    softSearchInv$pos2 <- round(softSearchInv$mid + softSearchInv$size/2)
    softSearchInvIrange <- GRanges(seqnames=softSearchInv$chromosome, 
                                   ranges=IRanges(start=softSearchInv$pos1, 
                                                  end=softSearchInv$pos2))
    softSearchInvIrangeRes <- findOverlaps(softSearchInvIrange, reduce(softSearchInvIrange))
    softSearchInv$clu <- subjectHits(softSearchInvIrangeRes)
    softSearchInvFilMer <- ddply(softSearchInv, ("clu"), softSearchCluster)
    if (nrow(softSearchInvFilMer)==0) {
      softSearchInvFilMer <- NULL
    } else {
      softSearchInvFilMer <- softSearchInvFilMer[, c("chromosome", "pos1", "pos2", "size")]
    }
  }
  
  ## filtering and merging Duplications
  softSearchDup <- softSearchPred[softSearchPred$type=="TDUP", ]
  DupTag <- which(unlist(strsplit(softSearchDup$format[1], ":"))=="TDUP")
  softSearchDup$ReadPairSupp <- sapply(strsplit(softSearchDup$detail,":"), 
                                       function(x){x[DupTag]})
  softTag <- which(unlist(strsplit(softSearchDup$format[1], ":"))=="nSC")
  softSearchDup$softSupp <- sapply(strsplit(softSearchDup$detail,":"), 
                                   function(x){x[softTag]})
  softSearchDup$softSupp[softSearchDup$softSupp=="NA"] <- 0
  softSearchDup$softSupp <- as.numeric(softSearchDup$softSupp)
  softSearchDup$ReadPairSupp <- as.numeric(softSearchDup$ReadPairSupp)
  softSearchDup <- softSearchDup[softSearchDup$ReadPairSupp>=readsSupport&
                                   softSearchDup$softSupp>=softClipsSupport,]
  softSearchDup$pos1 <- pmin(softSearchDup$start, softSearchDup$end)
  softSearchDup$pos2 <- pmax(softSearchDup$start, softSearchDup$end)
  softSearchDup <- softSearchDup[,c("chr1", "pos1", "pos2", "size", "ReadPairSupp")]
  names(softSearchDup)[1] <- "chromosome"
  
  if (nrow(softSearchDup)==0) {
    softSearchDupFilMer <- NULL
  } else {
    softSearchDup$size <- abs(softSearchDup$size)
    softSearchDup$mid <- (softSearchDup$pos1 + softSearchDup$pos2)/2
    softSearchDup$pos1 <- round(softSearchDup$mid - softSearchDup$size/2)
    softSearchDup$pos2 <- round(softSearchDup$mid + softSearchDup$size/2)
    softSearchDupIrange <- GRanges(seqnames=softSearchDup$chromosome, 
                                   ranges=IRanges(start=softSearchDup$pos1, 
                                                  end=softSearchDup$pos2))
    softSearchDupIrangeRes <- findOverlaps(softSearchDupIrange, reduce(softSearchDupIrange))
    softSearchDup$clu <- subjectHits(softSearchDupIrangeRes)
    softSearchDupFilMer <- ddply(softSearchDup, ("clu"), softSearchCluster)
    if (nrow(softSearchDupFilMer)==0) {
      softSearchDupFilMer <- NULL
    } else {
      softSearchDupFilMer <- softSearchDupFilMer[, c("chromosome", "pos1", "pos2", "size")]
    }
  }
  
  retuRes <- list(del=softSearchDelFilMer, dup=softSearchDupFilMer, 
                       inv=softSearchInvFilMer)
  attributes(retuRes) <- c(attributes(retuRes), list(method=method))
  
  return(retuRes);
}

