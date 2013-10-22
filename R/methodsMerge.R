
## Merging overlapped SVs predicted by different methods
methodsCluster <- function(df, methodsName, overLapPer=0.8, numMethodsSup=2)
{
    if(nrow(df)<2) {
        return(NULL)
    } else {
        dfTmp <- as.data.frame(t(as.matrix(df)), stringsAsFactors=F)
        DistMat <- sapply(dfTmp, function(x) {
            return(sapply(dfTmp, function(y) {
                OverlapLen <- min(as.numeric(as.character(x[3])), 
                    as.numeric(as.character(y[3]))) - 
                    max(as.numeric(as.character(x[2])), 
                    as.numeric(as.character(y[2])))
                OverlapPercent <- min(OverlapLen/(as.numeric(as.character(x[3]))-
                    as.numeric(as.character(x[2]))), 
                    OverlapLen/(as.numeric(as.character(y[3]))-
                    as.numeric(as.character(y[2]))))
                if (OverlapPercent<overLapPer) {
                    return(1)
                } else {
                    return(0)
                }
            }))
        })

        hc <- hclust(as.dist(DistMat))
        cl <- cutree(hc, h=0)
        df$cl <- unname(cl)
        dfRes <- ddply(df, ("cl"), function(dfTemp){
            if (length(unique(dfTemp$method))<numMethodsSup) {
                return(NULL)
            } else {
                tmpRes <- c(dfTemp$chromosome[1], round(mean(dfTemp$pos1)), 
                            round(mean(dfTemp$pos2)))
                MethodStat <- sapply(methodsName, function(x){
                                     return(ifelse(any(grepl(x, dfTemp$method)),
                                                   "Y", "N"))})
                return(c(tmpRes, MethodStat))
            }
        })
        return(dfRes)
    }
}



methodsMerge <- function(breakdancer, pindel, cnvnator, 
                         delly, svseq, others=NULL, overLapPer=0.8, numMethodsSup=2) 
{
    ## collecting all deletions predicted by different methods
    DeletionList <- list(breakdancer$del, pindel$del, 
                         cnvnator$del, delly$del, svseq$del)
    DeletionDf <- do.call(rbind, DeletionList)
    DeletionDf <- rbind(DeletionDf, others[others$type=="del", ][, 1:4])
    DeletionDf$method <- c(rep(c("breakdancer", "pindel", "cnvnator", "delly", 
        "svseq"), unlist(lapply(DeletionList, function(x){
            ifelse(is.null(nrow(x)), 0, nrow(x))}))), 
            others[others$type=="del", ]$methods)

    ## collecting all duplications predicted by different methods
    DuplicationList <- list(breakdancer$dup, pindel$dup, cnvnator$dup, 
                            delly$dup, svseq$dup)
    DuplicationDf <- do.call(rbind, DuplicationList)
    DuplicationDf <- rbind(DuplicationDf, others[others$type=="dup", ][, 1:4])
    DuplicationDf$method <- c(rep(c("breakdancer", "pindel", 
                                    "cnvnator", "delly", "svseq"), 
        unlist(lapply(DuplicationList, function(x){
            ifelse(is.null(nrow(x)), 0, nrow(x))}))), 
            others[others$type=="dup", ]$methods)

    ## collecting all inversions predicted by different methods
    InversionList <- list(breakdancer$inv, pindel$inv, cnvnator$inv, 
                          delly$inv, svseq$inv)
    InversionDf <- do.call(rbind, InversionList)
    InversionDf <- rbind(InversionDf, others[others$type=="inv", ][, 1:4]);
    InversionDf$method <- c(rep(c("breakdancer", "pindel", 
                                  "cnvnator", "delly", "svseq"), 
        unlist(lapply(InversionList, function(x){
            ifelse(is.null(nrow(x)), 0, nrow(x))}))), 
            others[others$type=="inv", ]$methods)

    MethodsName <- union(c("breakdancer", "pindel", 
                           "cnvnator", "delly", "svseq"), 
        names(table(others$method)))
    
    ## merging inversions predicted by different methods
    InversionIrange <- GRanges(seqnames=InversionDf$chromosome, 
        ranges=IRanges(start=InversionDf$pos1, end=InversionDf$pos2))
    InversionRes <- findOverlaps(InversionIrange, reduce(InversionIrange))
    InversionDf$class <- subjectHits(InversionRes)
    InversionDfMerge <- ddply(InversionDf, ("class"), 
                              methodsCluster, methodsName=MethodsName, 
                              overLapPer=overLapPer, numMethodsSup=numMethodsSup)
    InversionDfMerge$class <- NULL
    InversionDfMerge$cl <- NULL
    names(InversionDfMerge)[1:3] <- c("chromosome", "pos1", "pos2")
    InversionDfMerge$pos1 <- as.numeric(InversionDfMerge$pos1)
    InversionDfMerge$pos2 <- as.numeric(InversionDfMerge$pos2)

    ## merging deletions predicted by different methods
    DeletionIrange <- GRanges(seqnames=DeletionDf$chromosome, 
        ranges=IRanges(start=DeletionDf$pos1, end=DeletionDf$pos2))
    DeletionRes <- findOverlaps(DeletionIrange, reduce(DeletionIrange))
    DeletionDf$class <- subjectHits(DeletionRes)
    DeletionDfMerge <- ddply(DeletionDf, ("class"), 
                             methodsCluster, methodsName=MethodsName,
                             overLapPer=overLapPer, numMethodsSup=numMethodsSup)
    DeletionDfMerge$class <- NULL
    DeletionDfMerge$cl <- NULL
    names(DeletionDfMerge)[1:3] <- c("chromosome", "pos1", "pos2")
    DeletionDfMerge$pos1 <- as.numeric(DeletionDfMerge$pos1)
    DeletionDfMerge$pos2 <- as.numeric(DeletionDfMerge$pos2)


    ## merging duplications predicted by different methods
    DuplicationIrange <- GRanges(seqnames=DuplicationDf$chromosome, 
        ranges=IRanges(start=DuplicationDf$pos1, end=DuplicationDf$pos2))
    DuplicationRes <- findOverlaps(DuplicationIrange, reduce(DuplicationIrange))
    DuplicationDf$class <- subjectHits(DuplicationRes)
    DuplicationDfMerge <- ddply(DuplicationDf, ("class"), 
                                methodsCluster, methodsName=MethodsName,
                                overLapPer=overLapPer, numMethodsSup=numMethodsSup)
    DuplicationDfMerge$class <- NULL
    DuplicationDfMerge$cl <- NULL
    names(DuplicationDfMerge)[1:3] <- c("chromosome", "pos1", "pos2")
    DuplicationDfMerge$pos1 <- as.numeric(DuplicationDfMerge$pos1)
    DuplicationDfMerge$pos2 <- as.numeric(DuplicationDfMerge$pos2)

    return(list(del=DeletionDfMerge, dup=DuplicationDfMerge, 
                inv=InversionDfMerge))
}
