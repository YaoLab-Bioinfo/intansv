
## Visualization of SVs in a specific genomic region
plotRegion <- function(structuralVariation, genomeAnnotation, 
                       regionChromosome, regionStart, regionEnd)
{
  RegionNewIrange <- GRanges(seqnames=regionChromosome, 
                             IRanges(start=1, end=regionEnd-regionStart+1))
  seqlengths(RegionNewIrange) <- regionEnd - regionStart + 1
  p <- ggbio()
  
  RegionIrange <- GRanges(seqnames=regionChromosome, 
                          IRanges(start=regionStart, end=regionEnd))
  
  anno.gr <- GRanges(genomeAnnotation$chr[genomeAnnotation$tag=="gene"], 
                     IRanges(genomeAnnotation$start[genomeAnnotation$tag=="gene"],
                             genomeAnnotation$end[genomeAnnotation$tag=="gene"]),
                     Name=genomeAnnotation$ID[genomeAnnotation$tag=="gene"])
  
  ##  genes in specific region
  RegionTarget <- anno.gr[subjectHits(findOverlaps(RegionIrange, anno.gr))]
  RegionTargetGene <- RegionTarget
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
      p <- p + circle(RegionDupIrange, geom="rect", 
                      color="green4", fill="green4")
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
      p <- p + circle(RegionInvIrange, geom="rect", 
                      color="purple", fill="purple")
    }
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
      p <- p + circle(RegionDelIrange, geom="rect", 
                      color="red", fill="red")
    }
  }
  
  p <- p + circle(RegionNewIrange, geom="ideo", 
                  fill="gray70")
  
  if (length(RegionTargetGeneDfIrange)>0) {
    seqlengths(RegionTargetGeneDfIrange) <- regionEnd - regionStart + 1
    p <- p + circle(RegionTargetGeneDfIrange, geom="rect", 
                    color="blue", fill="blue")
    
    ##  genes names
    p <- p + circle(RegionTargetGeneDfIrange, geom="text", 
                    aes(label=name, label.text.angle=90), 
                    vjust=0, size=2)
  }
  
  return(p)
}


