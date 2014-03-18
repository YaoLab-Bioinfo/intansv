
## merging overlapped SVs predicted by DELLY
dellyCluster <- function(df) 
{
    maxReadPairSupp <- max(df$rp_support)
    dfFil <- df[df$rp_support>=(maxReadPairSupp/2), ]
    dfFilIrange <- IRanges(start=dfFil$rp_left, end=dfFil$rp_right)
    outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
    dfFil$cludelly <- subjectHits(outTmp)
    dfFilRes <- ddply(dfFil, ("cludelly"), function(x){
        if(nrow(x)==1){
            return(x)
        } else {
            LeftMin <- min(x$rp_left)
            RightMax <- max(x$rp_right)
            RangeLength <- RightMax-LeftMin
            x$op <- (x$rp_right-x$rp_left)/RangeLength
            if(any(x$op<0.8)){
                return(NULL)
            } else {
                return(x[which.max(x$rp_support), ])
            }
        }
    })
}

## Reading in the predicted SVs given by DELLY
readDelly <- function(dataDir=".", regSizeLowerCutoff=100, 
                      regSizeUpperCutoff=1000000, readsSupport=3,
                      method="DELLY") 
{
    ## reading deletions
    dellyDelRpList <- list.files(dataDir, full.names=T, pattern=".+\\.del$")
    dellyDelSrList <- list.files(dataDir, full.names=T, 
                                 pattern=".+\\.del\\.br$")

    dellyDelRp <- lapply(dellyDelRpList, function(x) {
        dellyData <- try(read.table(x, fill=T, as.is=T, comment.char="+"), silent=T)
        if (is.data.frame(dellyData)) {
            dellyData <- dellyData[grepl("Deletion", dellyData$V7), ]
            dellyData <- dellyData[, 1:7]
            names(dellyData) <- c("chr", "rp_left", "rp_right", "rp_size", 
                              "rp_support", "rp_map_quality", "id")
            return(dellyData)
        } else {
            return(NULL)
        }
    })

    dellyDelSr <- lapply(dellyDelSrList, function(x) {
        dellyData <- try(read.table(x, fill=T, col.names=paste("V", 1:13, sep=""), 
                                as.is=T), silent=T)
        if (is.data.frame(dellyData)) {
            dellyData <- dellyData[grepl("Deletion", dellyData$V7), ]
            dellyData <- dellyData[, 1:7]
            names(dellyData) <- c("chr", "sr_left", "sr_right", "sr_size", 
                              "sr_support", "sr_aln_quality", "id")
            return(dellyData)
        } else {
            return(NULL)
        }
    })

    dellyDelRpDf <- do.call(rbind, dellyDelRp)
    dellyDelSrDf <- do.call(rbind, dellyDelSr)

    dellyDelDf <- NULL
    if (is.data.frame(dellyDelRpDf) & is.data.frame(dellyDelSrDf)) {
        dellyDelDf <- merge(dellyDelRpDf, dellyDelSrDf, by=c("id", "chr"), all=T)
        dellyDelDf$rp_left <- as.numeric(dellyDelDf$rp_left)
        dellyDelDf$rp_right <- as.numeric(dellyDelDf$rp_right)
        dellyDelDf$rp_size <- as.numeric(dellyDelDf$rp_size)
        dellyDelDf <- dellyDelDf[dellyDelDf$rp_size>=regSizeLowerCutoff &
                             dellyDelDf$rp_size<=regSizeUpperCutoff &
                             dellyDelDf$rp_support>=readsSupport, ]
        dellyDelDf$name <- NULL
        dellyDelSrNna <- which(!(is.na(dellyDelDf$sr_size)))
        if (length(dellyDelSrNna)>0) {
            dellyDelDf[dellyDelSrNna, ]$rp_left <- dellyDelDf[dellyDelSrNna, ]$sr_left
            dellyDelDf[dellyDelSrNna, ]$rp_right <- dellyDelDf[dellyDelSrNna, ]$sr_right
            dellyDelDf[dellyDelSrNna, ]$rp_size <- dellyDelDf[dellyDelSrNna, ]$sr_size
        }
        dellyDelDf <- dellyDelDf[, 1:6]
        dellyDelDf$rp_left <- as.numeric(dellyDelDf$rp_left)
        dellyDelDf$rp_right <- as.numeric(dellyDelDf$rp_right)
        dellyDelDf$rp_size <- as.numeric(dellyDelDf$rp_size)
    } else if (is.data.frame(dellyDelRpDf)) {
        dellyDelDf <- dellyDelRpDf
        dellyDelDf$rp_left <- as.numeric(dellyDelDf$rp_left)
        dellyDelDf$rp_right <- as.numeric(dellyDelDf$rp_right)
        dellyDelDf$rp_size <- as.numeric(dellyDelDf$rp_size)
        dellyDelDf <- dellyDelDf[dellyDelDf$rp_size>=regSizeLowerCutoff &
                             dellyDelDf$rp_size<=regSizeUpperCutoff &
                             dellyDelDf$rp_support>=readsSupport, ]
        dellyDelDf$name <- NULL
        dellyDelDf <- dellyDelDf[, 1:6]
    } else if (is.data.frame(dellyDelSrDf)) {
        dellyDelDf <- dellyDelSrDf
        names(dellyDelDf)[2:5] <- c("rp_left","rp_right","rp_size","rp_support")
        dellyDelDf$rp_left <- as.numeric(dellyDelDf$rp_left)
        dellyDelDf$rp_right <- as.numeric(dellyDelDf$rp_right)
        dellyDelDf$rp_size <- as.numeric(dellyDelDf$rp_size)
        dellyDelDf <- dellyDelDf[dellyDelDf$rp_size>=regSizeLowerCutoff &
                             dellyDelDf$rp_size<=regSizeUpperCutoff &
                             dellyDelDf$rp_support>=readsSupport, ]
        dellyDelDf$name <- NULL
        dellyDelDf <- dellyDelDf[, 1:6]
    } else {
        dellyDelDf <- NULL
    }

    ## filtering and merging deletions
    if (!is.null(dellyDelDf)) {
        dellyDelIrange <- GRanges(seqnames=dellyDelDf$chr, 
                              ranges=IRanges(start=dellyDelDf$rp_left, 
                                             end=dellyDelDf$rp_right))
        dellyDelIrangeRes <- findOverlaps(dellyDelIrange, reduce(dellyDelIrange))
        dellyDelDf$clu <- subjectHits(dellyDelIrangeRes)
        dellyDelDfFilMer <- ddply(dellyDelDf, ("clu"), dellyCluster)
        dellyDelDfFilMer <- dellyDelDfFilMer[, 2:5]
        names(dellyDelDfFilMer) <- c("chromosome", "pos1", "pos2", "size")
    } else {
        dellyDelDfFilMer <- NULL
    }

    ## reading duplications
    dellyDupRpList <- list.files(dataDir, full.names=T, pattern=".+\\Dup$")
    dellyDupSrList <- list.files(dataDir, full.names=T, pattern=".+\\Dup\\.br$")

    dellyDupRp <- lapply(dellyDupRpList, function(x) {
        dellyData <- try(read.table(x, fill=T, as.is=T, comment.char="+"), silent=T)
        if (is.data.frame(dellyData)) {
            dellyData <- dellyData[grepl("Duplication", dellyData$V7), ]
            dellyData <- dellyData[, 1:7]
            names(dellyData) <- c("chr", "rp_left", "rp_right", "rp_size", 
                              "rp_support", "rp_map_quality", "id")
            return(dellyData)
        } else {
            return(NULL)
        }
    })

    dellyDupSr <- lapply(dellyDupSrList, function(x) {
        dellyData <- try(read.table(x, fill=T, col.names=paste("V", 1:13, sep=""), 
                                as.is=T), silent=T)
        if (is.data.frame(dellyData)) {
            dellyData <- dellyData[grepl("Duplication", dellyData$V7), ]
            dellyData <- dellyData[, 1:7]
            names(dellyData) <- c("chr", "sr_left", "sr_right", "sr_size", 
                              "sr_support", "sr_aln_quality", "id")
            return(dellyData)
        } else {
            return(NULL)
        }
    })

    dellyDupRpDf <- do.call(rbind, dellyDupRp)
    dellyDupSrDf <- do.call(rbind, dellyDupSr)

    dellyDupDf <- NULL
    if (is.data.frame(dellyDupRpDf) & is.data.frame(dellyDupSrDf)) {
        dellyDupDf <- merge(dellyDupRpDf, dellyDupSrDf, by=c("id", "chr"), all=T)
        dellyDupDf$rp_left <- as.numeric(dellyDupDf$rp_left)
        dellyDupDf$rp_right <- as.numeric(dellyDupDf$rp_right)
        dellyDupDf$rp_size <- as.numeric(dellyDupDf$rp_size)
        dellyDupDf <- dellyDupDf[dellyDupDf$rp_size>=regSizeLowerCutoff &
                             dellyDupDf$rp_size<=regSizeUpperCutoff &
                             dellyDupDf$rp_support>=readsSupport, ]
        dellyDupDf$name <- NULL
        dellyDupSrNna <- which(!(is.na(dellyDupDf$sr_size)))
        if (length(dellyDupSrNna)>0) {
            dellyDupDf[dellyDupSrNna, ]$rp_left <- dellyDupDf[dellyDupSrNna, ]$sr_left
            dellyDupDf[dellyDupSrNna, ]$rp_right <- dellyDupDf[dellyDupSrNna, ]$sr_right
            dellyDupDf[dellyDupSrNna, ]$rp_size <- dellyDupDf[dellyDupSrNna, ]$sr_size
        }
        dellyDupDf <- dellyDupDf[, 1:6]
        dellyDupDf$rp_left <- as.numeric(dellyDupDf$rp_left)
        dellyDupDf$rp_right <- as.numeric(dellyDupDf$rp_right)
        dellyDupDf$rp_size <- as.numeric(dellyDupDf$rp_size)
    } else if (is.data.frame(dellyDupRpDf)) {
        dellyDupDf <- dellyDupRpDf
        dellyDupDf$rp_left <- as.numeric(dellyDupDf$rp_left)
        dellyDupDf$rp_right <- as.numeric(dellyDupDf$rp_right)
        dellyDupDf$rp_size <- as.numeric(dellyDupDf$rp_size)
        dellyDupDf <- dellyDupDf[dellyDupDf$rp_size>=regSizeLowerCutoff &
                             dellyDupDf$rp_size<=regSizeUpperCutoff &
                             dellyDupDf$rp_support>=readsSupport, ]
        dellyDupDf$name <- NULL
        dellyDupDf <- dellyDupDf[, 1:6]
    } else if (is.data.frame(dellyDupSrDf)) {
        dellyDupDf <- dellyDupSrDf
        names(dellyDupDf)[2:5] <- c("rp_left","rp_right","rp_size","rp_support")
        dellyDupDf$rp_left <- as.numeric(dellyDupDf$rp_left)
        dellyDupDf$rp_right <- as.numeric(dellyDupDf$rp_right)
        dellyDupDf$rp_size <- as.numeric(dellyDupDf$rp_size)
        dellyDupDf <- dellyDupDf[dellyDupDf$rp_size>=regSizeLowerCutoff &
                             dellyDupDf$rp_size<=regSizeUpperCutoff &
                             dellyDupDf$rp_support>=readsSupport, ]
        dellyDupDf$name <- NULL
        dellyDupDf <- dellyDupDf[, 1:6]
    } else {
        dellyDupDf <- NULL
    }

    ## filtering and merging duplications
    if (!is.null(dellyDupDf)) {
        dellyDupIrange <- GRanges(seqnames=dellyDupDf$chr, 
                              ranges=IRanges(start=dellyDupDf$rp_left, 
                                             end=dellyDupDf$rp_right))
        dellyDupIrangeRes <- findOverlaps(dellyDupIrange, reduce(dellyDupIrange))
        dellyDupDf$clu <- subjectHits(dellyDupIrangeRes)
        dellyDupDfFilMer <- ddply(dellyDupDf, ("clu"), dellyCluster)
        dellyDupDfFilMer <- dellyDupDfFilMer[, 2:5]
        names(dellyDupDfFilMer) <- c("chromosome", "pos1", "pos2", "size")
    } else {
        dellyDupDfFilMer <- NULL
    }

    ## reading inversions
    dellyInvRp.list <- list.files(dataDir, full.names=T, pattern=".+\\.inv$")
    dellyInvSr.list <- list.files(dataDir, full.names=T, 
                                  pattern=".+\\.inv\\.br$")

    dellyInvRp <- lapply(dellyInvRp.list, function(x) {
        dellyData <- try(read.table(x, fill=T, as.is=T, comment.char="+"), silent=T)
        if (is.data.frame(dellyData)) {
            dellyData <- dellyData[grepl("Inversion", dellyData$V7), ]
            dellyData <- dellyData[, 1:7]
            names(dellyData) <- c("chr", "rp_left", "rp_right", 
                              "rp_size", "rp_support", "o_rp_map_qul", "id")
            return(dellyData)
        } else {
            return(NULL)
        }
    })
    dellyInvRpDf <- do.call(rbind, dellyInvRp)
    
    dellyInvSr <- lapply(dellyInvSr.list, function(x) {
        dellyData <- try(read.table(x, fill=T, col.names=paste("V", 1:13, sep=""), 
                                as.is=T), silent=T)
        if (is.data.frame(dellyData)) { 
            dellyData <- dellyData[grepl("Inversion", dellyData$V7), ]
            dellyData <- dellyData[, 1:7]
            names(dellyData) <- c("chr", "sr_left", "sr_right", "sr_size", 
                              "sr_supp", "sr_aln_qul", "id")
            return(dellyData)
        } else {
            return(NULL)
        }
    })
    dellyInvSrDf <- do.call(rbind, dellyInvSr)

    dellyInvDf <- NULL
    if (is.data.frame(dellyInvRpDf) & is.data.frame(dellyInvSrDf)) {
        dellyInvDf <- merge(dellyInvRpDf, dellyInvSrDf, by=c("id", "chr"), all=T)
        dellyInvDf$rp_left <- as.numeric(dellyInvDf$rp_left)
        dellyInvDf$rp_right <- as.numeric(dellyInvDf$rp_right)
        dellyInvDf$rp_size <- as.numeric(dellyInvDf$rp_size)
        dellyInvDf <- dellyInvDf[dellyInvDf$rp_size>=regSizeLowerCutoff &
                             dellyInvDf$rp_size<=regSizeUpperCutoff &
                             dellyInvDf$rp_support>=readsSupport, ]
        dellyInvDf$name <- NULL
        dellyInvSrNna <- which(!(is.na(dellyInvDf$sr_size)))
        if (length(dellyInvSrNna)>0) {
            dellyInvDf[dellyInvSrNna, ]$rp_left <- dellyInvDf[dellyInvSrNna, ]$sr_left
            dellyInvDf[dellyInvSrNna, ]$rp_right <- dellyInvDf[dellyInvSrNna, ]$sr_right
            dellyInvDf[dellyInvSrNna, ]$rp_size <- dellyInvDf[dellyInvSrNna, ]$sr_size
        }
        dellyInvDf <- dellyInvDf[, 1:6]
        dellyInvDf$rp_left <- as.numeric(dellyInvDf$rp_left)
        dellyInvDf$rp_right <- as.numeric(dellyInvDf$rp_right)
        dellyInvDf$rp_size <- as.numeric(dellyInvDf$rp_size)
    } else if (is.data.frame(dellyInvRpDf)) {
        dellyInvDf <- dellyInvRpDf
        dellyInvDf$rp_left <- as.numeric(dellyInvDf$rp_left)
        dellyInvDf$rp_right <- as.numeric(dellyInvDf$rp_right)
        dellyInvDf$rp_size <- as.numeric(dellyInvDf$rp_size)
        dellyInvDf <- dellyInvDf[dellyInvDf$rp_size>=regSizeLowerCutoff &
                             dellyInvDf$rp_size<=regSizeUpperCutoff &
                             dellyInvDf$rp_support>=readsSupport, ]
        dellyInvDf$name <- NULL
        dellyInvDf <- dellyInvDf[, 1:6]
    } else if (is.data.frame(dellyInvSrDf)) {
        dellyInvDf <- dellyInvSrDf
        names(dellyInvDf)[2:5] <- c("rp_left","rp_right","rp_size","rp_support")
        dellyInvDf$rp_left <- as.numeric(dellyInvDf$rp_left)
        dellyInvDf$rp_right <- as.numeric(dellyInvDf$rp_right)
        dellyInvDf$rp_size <- as.numeric(dellyInvDf$rp_size)
        dellyInvDf <- dellyInvDf[dellyInvDf$rp_size>=regSizeLowerCutoff &
                             dellyInvDf$rp_size<=regSizeUpperCutoff &
                             dellyInvDf$rp_support>=readsSupport, ]
        dellyInvDf$name <- NULL
        dellyInvDf <- dellyInvDf[, 1:6]
    } else {
        dellyInvDf <- NULL
    }
    
    ## filtering and merging inversions
    if (!is.null(dellyInvDf)) {
        dellyInvIrange <- GRanges(seqnames=dellyInvDf$chr, 
                              ranges=IRanges(start=dellyInvDf$rp_left, 
                                             end=dellyInvDf$rp_right))
        dellyInvIrangeRes <- findOverlaps(dellyInvIrange, reduce(dellyInvIrange))
        dellyInvDf$clu <- subjectHits(dellyInvIrangeRes)
        dellyInvDfFilMer <- ddply(dellyInvDf, ("clu"), dellyCluster)
        dellyInvDfFilMer <- dellyInvDfFilMer[, 2:5]
        names(dellyInvDfFilMer) <- c("chromosome", "pos1", "pos2", "size")
    } else {
        dellyInvDfFilMer <- NULL
    }

    retuRes <- list(del=dellyDelDfFilMer, dup=dellyDupDfFilMer, 
                inv=dellyInvDfFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes);
}
