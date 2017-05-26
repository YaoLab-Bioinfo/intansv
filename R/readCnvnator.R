

## Merging overlapped SVs predicted by CNVnator
CnvnatorCluster <- function(df)
{
    leftMin <- min(df$CNV_left)
    rightMax <- max(df$CNV_right)
    dfTmp <- NULL
    dfTmp$CNV_type <- df$CNV_type[1]
    dfTmp$length <- rightMax - leftMin
    dfTmp$chr <- df$chr[1]
    dfTmp$CNV_left <- leftMin
    dfTmp$CNV_right <- rightMax
    dfTmp$eval1 <- df$p_val1
    dfTmp$eval2 <- df$p_val2
    dfTmp$eval3 <- df$p_val3
    dfTmp$eval4 <- df$p_val4
    return(as.data.frame(dfTmp))
}

## Reading in the predicted SVs given by CNVnator
readCnvnator <- function(dataDir=".", regSizeLowerCutoff=100, 
                         regSizeUpperCutoff=1000000, method="CNVnator") 
{
    CnvnatorOutputList <- list.files(dataDir, full.names=T)
    CnvnatorPredc <- lapply(CnvnatorOutputList, function(x){
        dataTmp <- read.table(x, as.is=T)
        dataTmp$chr <- gsub(":.+", "", dataTmp$V2)
        dataTmp$CNV_left <- gsub("-\\d+", "", sub(".+:", "", dataTmp$V2))
        dataTmp$CNV_right <- gsub("\\d+-", "", sub(".+:", "", dataTmp$V2))
        dataTmp$V2 <- NULL
        names(dataTmp) <- c("CNV_type", "length", "normalized_RD", "p_val1", 
                            "p_val2", "p_val3", "p_val4", "q0", "chr", 
                            "CNV_left", "CNV_right")
        return(dataTmp)
    })

    CnvnatorRes <- do.call(rbind, CnvnatorPredc)
    CnvnatorRes$CNV_left <- as.numeric(CnvnatorRes$CNV_left)
    CnvnatorRes$CNV_right <- as.numeric(CnvnatorRes$CNV_right)
    CnvnatorRes <- CnvnatorRes[CnvnatorRes$length>=regSizeLowerCutoff & 
                               CnvnatorRes$length<=regSizeUpperCutoff, ]
    CnvnatorDel <- CnvnatorRes[CnvnatorRes$CNV_type=="deletion", ]
    CnvnatorDup <- CnvnatorRes[CnvnatorRes$CNV_type=="duplication", ]

    ## filtering and merging deletions
    CnvnatorDelIrange <- GRanges(seqnames=CnvnatorDel$chr, 
                                 ranges=IRanges(start=CnvnatorDel$CNV_left, 
                                                end=CnvnatorDel$CNV_right))
    CnvnatorDelIrangeRes <- findOverlaps(CnvnatorDelIrange, 
                                         reduce(CnvnatorDelIrange))
    CnvnatorDel$clu <- subjectHits(CnvnatorDelIrangeRes)
    CnvnatorDelFilMer <- ddply(CnvnatorDel, ("clu"), CnvnatorCluster)
    CnvnatorDelFilMer$clu <- NULL
    names(CnvnatorDelFilMer) <- c("CNV_type", "size", "chromosome", 
                                  "pos1", "pos2", "eval1", "eval2", "eval3", "eval4")
    CnvnatorDelFilMer <- CnvnatorDelFilMer[, c(3:5, 2, 6:9)]
    CnvnatorDelFilMer$info <- paste0("eval1=", CnvnatorDelFilMer$eval1, ";",
				     "eval2=", CnvnatorDelFilMer$eval2, ";",
				     "eval3=", CnvnatorDelFilMer$eval3, ";",
				     "eval4=", CnvnatorDelFilMer$eval4)
    CnvnatorDelFilMer <- CnvnatorDelFilMer[, -c(5:8)]
    CnvnatorDelFilMer$chromosome <- as.character(CnvnatorDelFilMer$chromosome)

    ## filtering and merging duplications
    CnvnatorDupIrange <- GRanges(seqnames=CnvnatorDup$chr, 
                                 ranges=IRanges(start=CnvnatorDup$CNV_left, 
                                                end=CnvnatorDup$CNV_right))
    CnvnatorDupIrangeRes <- findOverlaps(CnvnatorDupIrange, 
                                         reduce(CnvnatorDupIrange))
    CnvnatorDup$clu <- subjectHits(CnvnatorDupIrangeRes)
    CnvnatorDupFilMer <- ddply(CnvnatorDup, ("clu"), CnvnatorCluster)
    CnvnatorDupFilMer$clu <- NULL
    names(CnvnatorDupFilMer) <- c("CNV_type", "size", "chromosome", 
                                  "pos1", "pos2", "eval1", "eval2", "eval3", "eval4")
    CnvnatorDupFilMer <- CnvnatorDupFilMer[, c(3:5, 2, 6:9)]
    CnvnatorDupFilMer$info <- paste0("eval1=", CnvnatorDupFilMer$eval1, ";",
				     "eval2=", CnvnatorDupFilMer$eval2, ";",
				     "eval3=", CnvnatorDupFilMer$eval3, ";",
				     "eval4=", CnvnatorDupFilMer$eval4)
    CnvnatorDupFilMer <- CnvnatorDupFilMer[, -c(5:8)]
    CnvnatorDupFilMer$chromosome <- as.character(CnvnatorDupFilMer$chromosome)

    retuRes <- list(del=CnvnatorDelFilMer, dup=CnvnatorDupFilMer)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes);
}


