## Merging overlapped SVs predicted by Pindel
PindelCluster <- function(df) 
{
    maxReadPairSupp <- max(df$ReadPairSupport)
    dfFil <- df[df$ReadPairSupport>=(maxReadPairSupp/2), ]
    dfFilIrange <- IRanges(start=dfFil$BP_left, end=dfFil$BP_right)
    outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
    dfFil$clupindel <- subjectHits(outTmp)
    dfFilRes <- ddply(dfFil, ("clupindel"), function(x){
        if(nrow(x)==1){
            return(x)
        } else {
            LeftMin <- min(x$BP_left)
            RightMax <- max(x$BP_right)
            RangeLength <- RightMax-LeftMin
            x$op <- (x$BP_right-x$BP_left)/RangeLength
            if(any(x$op<0.8)){
                return(NULL)
            } else {
                return(x[which.max(x$ReadPairSupport), ])
            }
        }
    })
}

## Reading in the predicted SVs given by Pindel
readPindel <- function(dataDir=".", regSizeLowerCutoff=100, 
                       regSizeUpperCutoff=1000000, readsSupport=3,
                       method="Pindel") 
{
    PindelDelList <- list.files(dataDir, full.names=T, pattern=".+_D$")
    PindelInvList <- list.files(dataDir, full.names=T, pattern=".+_INV$")
    PindelTdList <- list.files(dataDir, full.names=T, pattern=".+_TD$")

    ## reading predicted deletions
    PindelDel <- lapply(PindelDelList, function(x){
        PdPredDf <- try(read.table(x, fill=T, as.is=T),silent=T)
        if (is.data.frame(PdPredDf)) {
            PdPredDf <- PdPredDf[grepl("^\\d", PdPredDf$V1), ]
            PdPredDf <- PdPredDf[, c(2, 3, 8, 10, 11, 13, 14, 16, 25)]
            names(PdPredDf) <- c("SV_type", "SV_len", "chromosome", "BP_left", 
                             "BP_right", "BP_range_left", "BP_range_right", 
                             "ReadPairSupport", "score")
        } else {
            PdPredDf <- NULL
        }
        return(PdPredDf)
    })

    ## reading predicted inversions
    PindelInv <- lapply(PindelInvList, function(x){
        PdPredDf <- try(read.table(x, fill=T, quote="", as.is=T),silent=T)
        if (is.data.frame(PdPredDf)) {
            PdPredDf <- PdPredDf[grepl("^\\d", PdPredDf$V1), ]
            PdPredDf <- PdPredDf[, c(2, 3, 8, 10, 11, 13, 14, 16, 25)]
            names(PdPredDf) <- c("SV_type", "SV_len", "chromosome", "BP_left", 
                             "BP_right", "BP_range_left", "BP_range_right", 
                             "ReadPairSupport", "score")
        } else {
            PdPredDf <- NULL
        }
        return(PdPredDf)
    })

    ## reading predicted tandom duplications
    PindelTd <- lapply(PindelTdList, function(x){
        PdPredDf <- try(read.table(x, fill=T, as.is=T),silent=T)
        if (is.data.frame(PdPredDf)) {
            PdPredDf <- PdPredDf[grepl("^\\d", PdPredDf$V1), ]
            PdPredDf <- PdPredDf[, c(2, 3, 8, 10, 11, 13, 14, 16, 25)]
            names(PdPredDf) <- c("SV_type", "SV_len", "chromosome", "BP_left", 
                             "BP_right", "BP_range_left", "BP_range_right", 
                             "ReadPairSupport", "score")
        } else {
            PdPredDf <- NULL
        }
        return(PdPredDf)
    })

    ## merging predictions of different chromosome
    PindelDelDf <- do.call(rbind, PindelDel)
    PindelDelDf$SV_len <- as.numeric(PindelDelDf$SV_len)
    PindelDelDf <- PindelDelDf[PindelDelDf$SV_len>=regSizeLowerCutoff &
                               PindelDelDf$SV_len<=regSizeUpperCutoff &
                               PindelDelDf$ReadPairSupport>=readsSupport, ]
    PindelInvDf <- do.call(rbind, PindelInv)
    PindelInvDf <- PindelInvDf[PindelInvDf$SV_len>=regSizeLowerCutoff &
                               PindelInvDf$SV_len<=regSizeUpperCutoff &
                               PindelInvDf$ReadPairSupport>=readsSupport, ]
    PindelTdDf <- do.call(rbind, PindelTd)
    PindelTdDf <- PindelTdDf[PindelTdDf$SV_len>=regSizeLowerCutoff &
                             PindelTdDf$SV_len<=regSizeUpperCutoff &
                             PindelTdDf$ReadPairSupport>=readsSupport, ]

    ## filtering and merging deletions
    if (is.data.frame(PindelDelDf)) {
        PindelDelIrange <- GRanges(seqnames=PindelDelDf$chromosome, 
                               ranges=IRanges(start=PindelDelDf$BP_left, 
                                              end=PindelDelDf$BP_right))
        PindelDelIrangeRes <- findOverlaps(PindelDelIrange, reduce(PindelDelIrange))
        PindelDelDf$clu <- subjectHits(PindelDelIrangeRes)
        PindelDelDfFilMer <- ddply(PindelDelDf, ("clu"), PindelCluster)
        PindelDelDfFilMer <- PindelDelDfFilMer[, c(3:5, 2)]
        names(PindelDelDfFilMer)[2:4] <- c("pos1", "pos2", "size")
        PindelDelDfFilMer$size <- as.numeric(PindelDelDfFilMer$size)
        PindelDelDfFilMer$pos1 <- as.numeric(PindelDelDfFilMer$pos1)
        PindelDelDfFilMer$pos2 <- as.numeric(PindelDelDfFilMer$pos2)
    } else {
        PindelDelDfFilMer <- NULL
    }

    ## filtering and merging inversions
    if (is.data.frame(PindelInvDf)) {
        PindelInvIrange <- GRanges(seqnames=PindelInvDf$chromosome, 
                               ranges=IRanges(start=PindelInvDf$BP_left, 
                                              end=PindelInvDf$BP_right))
        PindelInvIrangeRes <- findOverlaps(PindelInvIrange, reduce(PindelInvIrange))
        PindelInvDf$clu <- subjectHits(PindelInvIrangeRes)
        PindelInvDfFilMer <- ddply(PindelInvDf, ("clu"), PindelCluster)
        PindelInvDfFilMer <- PindelInvDfFilMer[, c(3:5, 2)]
        names(PindelInvDfFilMer)[2:4] <- c("pos1", "pos2", "size")
        PindelInvDfFilMer$size <- as.numeric(PindelInvDfFilMer$size)
        PindelInvDfFilMer$pos1 <- as.numeric(PindelInvDfFilMer$pos1)
        PindelInvDfFilMer$pos2 <- as.numeric(PindelInvDfFilMer$pos2)
    } else {
        PindelInvDfFilMer <- NULL
    }

    ## filtering and merging tandom duplications
    if (is.data.frame(PindelTdDf)) {
        PindelTdIrange <- GRanges(seqnames=PindelTdDf$chromosome, 
                              ranges=IRanges(start=PindelTdDf$BP_left, 
                                             end=PindelTdDf$BP_right))
        PindelTdIrangeRes <- findOverlaps(PindelTdIrange, reduce(PindelTdIrange))
        PindelTdDf$clu <- subjectHits(PindelTdIrangeRes)
        PindelTdDfFilMer <- ddply(PindelTdDf, ("clu"), PindelCluster)
        PindelTdDfFilMer <- PindelTdDfFilMer[, c(3:5, 2)]
        names(PindelTdDfFilMer)[2:4] <- c("pos1", "pos2", "size")
        PindelTdDfFilMer$size <- as.numeric(PindelTdDfFilMer$size)
        PindelTdDfFilMer$pos1 <- as.numeric(PindelTdDfFilMer$pos1)
        PindelTdDfFilMer$pos2 <- as.numeric(PindelTdDfFilMer$pos2)
    } else {
        PindelTdDfFilMer <- NULL
    }

    retuRes <- list(del=PindelDelDfFilMer, inv=PindelInvDfFilMer, 
                dup=PindelTdDfFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes);
}
