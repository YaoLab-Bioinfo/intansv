
## Visualization of SVs in a specific genomic region
plotRegion <- function(structuralVariation, genomeAnnotation, 
                       regionChromosome, regionStart, regionEnd)
{
  RegionNewIrange <- GRanges(seqnames=regionChromosome, 
                             IRanges(start=1, end=regionEnd-regionStart+1))
  seqlengths(RegionNewIrange) <- regionEnd - regionStart + 1
  p <- ggplot() + layout_circle(RegionNewIrange, geom="ideo", 
                                fill="gray70", radius=25, trackWidth=2)
  
  RegionIrange <- GRanges(seqnames=regionChromosome, 
                          IRanges(start=regionStart, end=regionEnd))
  ##  genes in specific region
  RegionTarget <- genomeAnnotation[subjectHits(findOverlaps(RegionIrange, 
                                                            genomeAnnotation))]
  RegionTargetGene <- RegionTarget[RegionTarget$type=="gene"]
  RegionTargetGeneDf <- NULL
  RegionTargetGeneDf$chr <- rep(seqnames(RegionTargetGene)@values, 
                                seqnames(RegionTargetGene)@lengths)
  RegionTargetGeneDf$start <- start(ranges(RegionTargetGene))
  RegionTargetGeneDf$end <- end(ranges(RegionTargetGene))
  RegionTargetGeneDf$loc <- RegionTargetGene$Name
  RegionTargetGeneDf <- as.data.frame(RegionTargetGeneDf)
  RegionTargetGeneDf <- RegionTargetGeneDf[RegionTargetGeneDf$start>=regionStart&
                                             RegionTargetGeneDf$end<=regionEnd, ]
  RegionTargetGeneDf$start <- RegionTargetGeneDf$start - regionStart + 1
  RegionTargetGeneDf$end <- RegionTargetGeneDf$end - regionStart + 1
  
  RegionTargetGeneDfIrange <- 
    GRanges(seqnames=as.character(RegionTargetGeneDf$chr), 
            IRanges(RegionTargetGeneDf$start, 
                    RegionTargetGeneDf$end), 
            name=RegionTargetGeneDf$loc)
  
  if (length(RegionTargetGeneDfIrange)>0) {
    seqlengths(RegionTargetGeneDfIrange) <- regionEnd - regionStart + 1
    p <- p + layout_circle(RegionTargetGeneDfIrange, geom="rect", 
                           color="blue", fill="blue", radius=29, trackWidth=2)
    ##  genes names
    p <- p + layout_circle(RegionTargetGeneDfIrange, geom="text", 
                           aes_string(label='name', angle=90), 
                           radius=33, trackWidth=6, vjust=0, size=2)
  }
  
  ##  deletions in specific region
  if (any(names(structuralVariation)=="del")) {
    RegionDel <- 
      structuralVariation$del[
        (structuralVariation$del)$chromosome==regionChromosome &
          (structuralVariation$del)$pos1>=regionStart &
          (structuralVariation$del)$pos2<=regionEnd, ]
    RegionDel$pos1 <- RegionDel$pos1 - regionStart + 1
    RegionDel$pos2 <- RegionDel$pos2 - regionStart + 1
    
    if (nrow(RegionDel)>0) {
      RegionDelIrange <- GRanges(seqnames=RegionDel$chromosome, 
                                 IRanges(RegionDel$pos1, RegionDel$pos2))
      seqlengths(RegionDelIrange) <- regionEnd - regionStart + 1
      p <- p + layout_circle(RegionDelIrange, geom="rect", 
                             color="red", fill="red", radius=22, 
                             trackWidth=2)
    }
  }
  
  ##  inversions in specific region
  if (any(names(structuralVariation)=="inv")) {
    RegionInv <- 
      structuralVariation$inv[(structuralVariation$inv)$chromosome==regionChromosome & 
                                (structuralVariation$inv)$pos1>=regionStart & 
                                (structuralVariation$inv)$pos2<=regionEnd, ]
    RegionInv$pos1 <- RegionInv$pos1 - regionStart + 1
    RegionInv$pos2 <- RegionInv$pos2 - regionStart + 1
    if (nrow(RegionInv)>0) {
      RegionInvIrange <- GRanges(seqnames=RegionInv$chromosome, 
                                 IRanges(RegionInv$pos1, RegionInv$pos2))
      seqlengths(RegionInvIrange) <- regionEnd - regionStart + 1
      p <- p + layout_circle(RegionInvIrange, geom="rect", 
                             color="purple", fill="purple", 
                             radius=16, trackWidth=2)
    }
  }
  
  ##  duplications in specific region
  if (any(names(structuralVariation)=="dup")) {
    RegionDup <- 
      structuralVariation$dup[(structuralVariation$dup)$chromosome==regionChromosome & 
                                (structuralVariation$dup)$pos1>=regionStart & 
                                (structuralVariation$dup)$pos2<=regionEnd, ]
    RegionDup$pos1 <- RegionDup$pos1 - regionStart + 1
    RegionDup$pos2 <- RegionDup$pos2 - regionStart + 1
    if (nrow(RegionDup)>0) {
      RegionDupIrange <- GRanges(seqnames=RegionDup$chromosome, 
                                 IRanges(RegionDup$pos1, RegionDup$pos2))
      seqlengths(RegionDupIrange) <- regionEnd - regionStart + 1
      p <- p + layout_circle(RegionDupIrange, geom="rect", 
                             color="green4", fill="green4", 
                             radius=19, trackWidth=2)
    }
  }
  
  return(p)
}


