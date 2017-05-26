
## merging overlapped SVs predicted by DELLY
dellyCluster <- function(df) 
{
    maxReadPairSupp <- max(df$rp_support)
    dfFil <- df[df$rp_support>=(maxReadPairSupp/2), ]
    dfFilIrange <- IRanges(start=dfFil$start, end=dfFil$end)
    outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
    dfFil$cludelly <- subjectHits(outTmp)
    dfFilRes <- ddply(dfFil, ("cludelly"), function(x){
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
readDelly <- function(dataDir=".", regSizeLowerCutoff=100, 
                      regSizeUpperCutoff=1000000, readsSupport=3,
                      method="DELLY", pass=TRUE, minMappingQuality=20) 
{
    ## reading in SV predictions
    dellyFileList <- list.files(dataDir, full.names=T)

    dellyFileCont <- lapply(dellyFileList, function(x) {
        dellyData <- try(read.table(x, as.is=T), silent=T)
        if (is.data.frame(dellyData)) {
            dellyData <- dellyData[, c(1:2, 7:8, 10)]
            names(dellyData) <- c("chr", "start", "pass", 
                              "info", "detail")
            return(dellyData)
        } else {
            return(NULL)
        }
    })

    dellyCont <- do.call(rbind, dellyFileCont)
    
    if (pass) {
      dellyCont <- dellyCont[dellyCont$pass=="PASS", ]
    }
    
    dellyCont$end <- as.numeric(gsub(".+;END=(\\d+);.+", "\\1", dellyCont$info))
    dellyCont$type <- gsub(".+;SVTYPE=([A-Z]+);.+", "\\1", dellyCont$info)
    dellyCont$pe_support <- as.numeric(gsub(".+;PE=(\\d+);.+", "\\1", dellyCont$info))
    dellyCont$sr_support <- 0
    dellyCont$sr_support[grepl(";SR=", dellyCont$info)] <- 
              as.numeric(gsub(".+;SR=(\\d+);.+", "\\1", dellyCont$info[grepl(";SR=", dellyCont$info)]))
    
    dellyCont$rp_support <- dellyCont$pe_support + dellyCont$sr_support
    dellyCont$map_quality <- as.numeric(gsub(".+;MAPQ=(\\d+).*", "\\1", dellyCont$info))
    dellyCont$size <- as.numeric(dellyCont$end - dellyCont$start)
    
    dellyCont <- dellyCont[dellyCont$map_quality>=minMappingQuality&
                          dellyCont$rp_support>=readsSupport&
                            dellyCont$size>=regSizeLowerCutoff&
                            dellyCont$size<=regSizeUpperCutoff, ]
    
    dellyCont <- dellyCont[, c("chr", "start", "end", "size", "rp_support",
                               "pe_support", "sr_support", "map_quality", "type")]

    dellyDelDf <- dellyCont[dellyCont$type=="DEL", ]

    ## filtering and merging deletions
    if (!is.null(dellyDelDf) && nrow(dellyDelDf)>0) {
        dellyDelIrange <- GRanges(seqnames=dellyDelDf$chr, 
                              ranges=IRanges(start=dellyDelDf$start, 
                                             end=dellyDelDf$end))
        dellyDelIrangeRes <- findOverlaps(dellyDelIrange, reduce(dellyDelIrange))
        dellyDelDf$clu <- subjectHits(dellyDelIrangeRes)
        dellyDelDfFilMer <- ddply(dellyDelDf, ("clu"), dellyCluster)
        dellyDelDfFilMer <- dellyDelDfFilMer[, c("chr", "start", "end", "size", "rp_support",
                                           "pe_support", "sr_support")]
        names(dellyDelDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport",
                                           "pe_support", "sr_support")
	dellyDelDfFilMer$readsSupport <- NULL
	dellyDelDfFilMer$info <- paste0("PE=", dellyDelDfFilMer$pe_support, ";",
				       "SR=", dellyDelDfFilMer$sr_support)
	dellyDelDfFilMer$pe_support <- NULL
	dellyDelDfFilMer$sr_support <- NULL
    } else {
        dellyDelDfFilMer <- NULL
    }

    ## reading duplications
    dellyDupDf <- dellyCont[dellyCont$type=="DUP", ]

    ## filtering and merging duplications
    if (!is.null(dellyDupDf) && nrow(dellyDupDf)>0) {
        dellyDupIrange <- GRanges(seqnames=dellyDupDf$chr, 
                              ranges=IRanges(start=dellyDupDf$start, 
                                             end=dellyDupDf$end))
        dellyDupIrangeRes <- findOverlaps(dellyDupIrange, reduce(dellyDupIrange))
        dellyDupDf$clu <- subjectHits(dellyDupIrangeRes)
        dellyDupDfFilMer <- ddply(dellyDupDf, ("clu"), dellyCluster)
        dellyDupDfFilMer <- dellyDupDfFilMer[, c("chr", "start", "end", "size", "rp_support",
                                           "pe_support", "sr_support")]
        names(dellyDupDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport",
                                           "pe_support", "sr_support")
	dellyDupDfFilMer$readsSupport <- NULL
	dellyDupDfFilMer$info <- paste0("PE=", dellyDupDfFilMer$pe_support, ";",
				       "SR=", dellyDupDfFilMer$sr_support)
	dellyDupDfFilMer$pe_support <- NULL
	dellyDupDfFilMer$sr_support <- NULL
    } else {
        dellyDupDfFilMer <- NULL
    }

    ## reading inversions
    dellyInvDf <- dellyCont[dellyCont$type=="INV", ]
    
    ## filtering and merging inversions
    if (!is.null(dellyInvDf) && nrow(dellyInvDf)>0) {
        dellyInvIrange <- GRanges(seqnames=dellyInvDf$chr, 
                              ranges=IRanges(start=dellyInvDf$start, 
                                             end=dellyInvDf$end))
        dellyInvIrangeRes <- findOverlaps(dellyInvIrange, reduce(dellyInvIrange))
        dellyInvDf$clu <- subjectHits(dellyInvIrangeRes)
        dellyInvDfFilMer <- ddply(dellyInvDf, ("clu"), dellyCluster)
        dellyInvDfFilMer <- dellyInvDfFilMer[, c("chr", "start", "end", "size", "rp_support",
                                           "pe_support", "sr_support")]
        names(dellyInvDfFilMer) <- c("chromosome", "pos1", "pos2", "size", "readsSupport",
                                           "pe_support", "sr_support")
	dellyInvDfFilMer$readsSupport <- NULL
	dellyInvDfFilMer$info <- paste0("PE=", dellyInvDfFilMer$pe_support, ";",
				       "SR=", dellyInvDfFilMer$sr_support)
	dellyInvDfFilMer$pe_support <- NULL
	dellyInvDfFilMer$sr_support <- NULL
    } else {
        dellyInvDfFilMer <- NULL
    }

    retuRes <- list(del=dellyDelDfFilMer, dup=dellyDupDfFilMer, 
                inv=dellyInvDfFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes)
}
