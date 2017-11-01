# A function to calculate the overlaping percentage between query and subject in GRanges or IRange format
tellOLPercantage <- function(query, subject, ...) {
  if(!((class(query)=="GRanges" & class(subject)=="GRanges") | (class(query)=="IRanges" & class(subject)=="IRanges"))) {
    stop("Query and subject need to be of same class, either GRanges or IRanges!")
  }
  hits <- findOverlaps(query, subject)
  overlaps <- pintersect(query[queryHits(hits)],subject[subjectHits(hits)])
  olp_percentage <- width(overlaps)/width(subject[subjectHits(hits)])
  return(olp_percentage)
}


# MergeOLCNVs: A function to merge overlapping CNVs similar in size
mergeOLCNVs <- function(df, percentage = 0.6, software="", ...){
  if(!class(df)=="data.frame"){
    stop("a dataframe is required as the input")
  }

  GRSubj <- makeGRangesFromDataFrame(df,keep.extra.columns = T)
  self_overlap <- countOverlaps(GRSubj, GRSubj , type="any")
  
  UniqCNVs    <- GRSubj[self_overlap == 1,] #to be returned
  NonUniqCNVs <- GRSubj[self_overlap >  1,] 
  
  #self_overlap_cycle_2 (so2)
  so2<- tellOLPercantage(NonUniqCNVs, reduce(NonUniqCNVs))
  
  #PassedCNVs_in_cycle_2 (Pass_2)
  #FailedCNVs_in_cycle_2 (Fail_2)
  
  Pass2nd <- NonUniqCNVs[so2 <= percentage, ] #to be returned
  Fail2nd <- NonUniqCNVs[so2 >  percentage, ] 
  
  Fail2nd$clu <- subjectHits(findOverlaps(Fail2nd, reduce(Fail2nd)))

  #Passed_CNVs_in_cycle_3 (pc3)
  Pass3rd <- GRanges()
  for (i in unique(Fail2nd$clu)){
    CNV_in_clu <- Fail2nd[Fail2nd$clu==i,]
    if(length(CNV_in_clu) == 1){
      CNV_in_clu$clu <- NULL
      Pass3rd<- append(Pass3rd, CNV_in_clu)
    }else{
      best_supported <- CNV_in_clu[which.max(CNV_in_clu$rp_support), ]
      best_supported$clu <- NULL
      Pass3rd <- append(Pass3rd, best_supported)
    }
  }
  
  FinalSet <- GRanges()
  FinalSet <- append(FinalSet,UniqCNVs)
  FinalSet <- append(FinalSet,Pass2nd)
  FinalSet <- append(FinalSet,Pass3rd)
  if(software == "DELLY"){
    info_str <- paste0("PE=", FinalSet$pe_support, ";", "SR=", FinalSet$sr_support)
    }else if(software == "BreakDancer"){
    info_str <- paste0("score=", FinalSet$score, ";", "PE=", FinalSet$rp_support)
    }

  reformt_df <- data.frame(chromosome = seqnames(FinalSet), 
                           pos1 = start(FinalSet), 
                           pos2 = end(FinalSet),
                           size = FinalSet$size,
                           info = info_str,
                           stringsAsFactors=FALSE
                           )
  reformt_df$chromosome = as.character(reformt_df$chromosome)
  return(reformt_df)
}
