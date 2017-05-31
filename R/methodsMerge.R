
## Merging overlapped SVs predicted by different methods
methodsCluster <- function(df, methodsName, overLapPer=0.8, numMethodsSup=2)
{
  if ( (nrow(df)<2) & (numMethodsSup>=2) ) {
    return(NULL)
  } else if ( (nrow(df)<2) & (numMethodsSup==1) ) {
    MethodStat <- paste(df$method, sep="", collapse = ":")
    infoStat <- paste(df$info, sep="", collapse = ":")
    return(as.data.frame(t(c(df$chromosome[1], df$pos1[1], df$pos2[1],
                             MethodStat, infoStat)), stringsAsFactors=FALSE))
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
        MethodStat <- paste(dfTemp$method, sep="", collapse = ":")
        infoStat <- paste(dfTemp$info, sep=":", collapse = ":")
        return(c(tmpRes, MethodStat, infoStat))
      }
    })
    return(dfRes)
  }
}


methodsMerge <- function(..., others=NULL, overLapPerDel=0.8, 
                         overLapPerDup=0.8, overLapPerInv=0.8, 
                         numMethodsSupDel=2, numMethodsSupDup=2, numMethodsSupInv=2) 
{
    svIn <- list(...);
    
    ## collecting all deletions predicted by different methods
    DeletionList <- lapply(svIn, function(x){return(x$del)})
    DeletionDf <- do.call(rbind, DeletionList)
    svInMethod <- sapply(svIn, function(x){return(attr(x, "method"))})
    DeletionDf <- rbind(DeletionDf, others[others$type=="del", ][, 1:4])
    DeletionDf$method <- c(rep(svInMethod, unlist(lapply(DeletionList, function(x){
            ifelse(is.null(nrow(x)), 0, nrow(x))}))), 
            others[others$type=="del", ]$methods)

    ## collecting all duplications predicted by different methods
    DuplicationList <- lapply(svIn, function(x){return(x$dup)})
    DuplicationDf <- do.call(rbind, DuplicationList)
    DuplicationDf <- rbind(DuplicationDf, others[others$type=="dup", ][, 1:4])
    DuplicationDf$method <- c(rep(svInMethod, 
        unlist(lapply(DuplicationList, function(x){
            ifelse(is.null(nrow(x)), 0, nrow(x))}))), 
            others[others$type=="dup", ]$methods)

    ## collecting all inversions predicted by different methods
    InversionList <- lapply(svIn, function(x){return(x$inv)})
    InversionDf <- do.call(rbind, InversionList)
    InversionDf <- rbind(InversionDf, others[others$type=="inv", ][, 1:4]);
    InversionDf$method <- c(rep(svInMethod, 
        unlist(lapply(InversionList, function(x){
            ifelse(is.null(nrow(x)), 0, nrow(x))}))), 
            others[others$type=="inv", ]$methods)

    MethodsName <- union(svInMethod, 
        names(table(others$method)))
    
    ## merging inversions predicted by different methods
    InversionIrange <- GRanges(seqnames=InversionDf$chromosome, 
        ranges=IRanges(start=InversionDf$pos1, end=InversionDf$pos2))
    InversionRes <- findOverlaps(InversionIrange, reduce(InversionIrange))
    InversionDf$class <- subjectHits(InversionRes)
    if (!is.null(nrow(InversionDf)) && nrow(InversionDf)>0) {
      InversionDfMerge <- ddply(InversionDf, ("class"), 
                                methodsCluster, methodsName=MethodsName, 
                                overLapPer=overLapPerInv, numMethodsSup=numMethodsSupInv)
      if (nrow(InversionDfMerge)>0) {
        InversionDfMerge$class <- NULL
        InversionDfMerge$cl <- NULL
        names(InversionDfMerge)[1:5] <- c("chromosome", "pos1", "pos2", "methods", "info")
        InversionDfMerge$pos1 <- as.numeric(InversionDfMerge$pos1)
        InversionDfMerge$pos2 <- as.numeric(InversionDfMerge$pos2)
      } else {
        InversionDfMerge <- NULL
      }
    } else {
      InversionDfMerge <- NULL
    }
    

    ## merging deletions predicted by different methods
    DeletionIrange <- GRanges(seqnames=DeletionDf$chromosome, 
        ranges=IRanges(start=DeletionDf$pos1, end=DeletionDf$pos2))
    DeletionRes <- findOverlaps(DeletionIrange, reduce(DeletionIrange))
    DeletionDf$class <- subjectHits(DeletionRes)
    if (!is.null(nrow(DeletionDf)) && nrow(DeletionDf)>0) {
      DeletionDfMerge <- ddply(DeletionDf, ("class"), 
                               methodsCluster, methodsName=MethodsName,
                               overLapPer=overLapPerDel, numMethodsSup=numMethodsSupDel)
      if (nrow(DeletionDfMerge)>0) {
        DeletionDfMerge$class <- NULL
        DeletionDfMerge$cl <- NULL
        names(DeletionDfMerge)[1:5] <- c("chromosome", "pos1", "pos2", "methods", "info")
        DeletionDfMerge$pos1 <- as.numeric(DeletionDfMerge$pos1)
        DeletionDfMerge$pos2 <- as.numeric(DeletionDfMerge$pos2)
      } else {
        DeletionDfMerge <- NULL
      }
    } else {
      DeletionDfMerge <- NULL
    }

    ## merging duplications predicted by different methods
    DuplicationIrange <- GRanges(seqnames=DuplicationDf$chromosome, 
        ranges=IRanges(start=DuplicationDf$pos1, end=DuplicationDf$pos2))
    DuplicationRes <- findOverlaps(DuplicationIrange, reduce(DuplicationIrange))
    DuplicationDf$class <- subjectHits(DuplicationRes)
    if (!is.null(nrow(DuplicationDf)) && nrow(DuplicationDf)>0) {
      DuplicationDfMerge <- ddply(DuplicationDf, ("class"), 
                                  methodsCluster, methodsName=MethodsName,
                                  overLapPer=overLapPerDup, numMethodsSup=numMethodsSupDup)
      if (nrow(DuplicationDfMerge)>0) {
        DuplicationDfMerge$class <- NULL
        DuplicationDfMerge$cl <- NULL
        names(DuplicationDfMerge)[1:5] <- c("chromosome", "pos1", "pos2", "methods", "info")
        DuplicationDfMerge$pos1 <- as.numeric(DuplicationDfMerge$pos1)
        DuplicationDfMerge$pos2 <- as.numeric(DuplicationDfMerge$pos2)
      } else {
        DuplicationDfMerge <- NULL
      }
    } else {
      DuplicationDfMerge <- NULL
    }

    return(list(del=DeletionDfMerge, dup=DuplicationDfMerge, 
                inv=InversionDfMerge))
}
