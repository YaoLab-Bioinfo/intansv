
svAnnotation <- function(structuralVariation, genomeAnnotation) 
{
  if (is.null(structuralVariation)) {
    return(NULL)
  } else if (!is.data.frame(structuralVariation)) {
    stop("structuralVariation should be a data frame!\n")
  } else if (!identical(names(structuralVariation)[1:3], 
                        c("chromosome", "pos1", "pos2"))) {
    stop("structuralVariation should have appropriate names!\n")
  } else if (!is.data.frame(genomeAnnotation)) {
    stop("genomeAnnotation should be a data frame!\n")
  } else if (nrow(genomeAnnotation) == 0) {
    stop("The content of genomeAnnotation is NULL!\n")
  }
  
  if (is.null(structuralVariation) || nrow(structuralVariation)<=0) {
    return(NULL)
  } else {
    
    sv.gr <- GRanges(seqnames=structuralVariation$chromosome, 
                     ranges=IRanges(start=as.numeric(structuralVariation$pos1), 
                                    end=as.numeric(structuralVariation$pos2)))
    anno.gr <- GRanges(genomeAnnotation$chr, IRanges(genomeAnnotation$start,
                                                     genomeAnnotation$end))
    sv.anno.op <- findOverlaps(sv.gr, anno.gr)
    sv.anno.op.df <- data.frame(sv.anno.op)
    names(sv.anno.op.df) <- c("que", "sub")
    
    structuralVariation$que <- 1:nrow(structuralVariation)
    genomeAnnotation$sub <- 1:nrow(genomeAnnotation)
    sv.anno <- merge(genomeAnnotation, sv.anno.op.df, by="sub")
    sv.anno <- merge(structuralVariation, sv.anno, by="que", all.x=TRUE)
    
    sv.anno$que <- NULL
    sv.anno$sub <- NULL
    sv.anno$chr <- NULL
    
    return(sv.anno)
  }
}

