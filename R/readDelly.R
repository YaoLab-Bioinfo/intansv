## Reading in the predicted SVs given by DELLY
readDelly <- function(dataDir=".", regSizeLowerCutoff=100, 
                      regSizeUpperCutoff=10000000, readsSupport=3,
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
    
    
    dellyDelDf <- mergeOLCNVs(dellyCont[dellyCont$type =="DEL", ], software= method)
    dellyDupDf <- mergeOLCNVs(dellyCont[dellyCont$type =="DUP", ], software= method)
    dellyINVDf <- mergeOLCNVs(dellyCont[dellyCont$type =="INV", ], software= method)


    retuRes <- list(del=dellyDelDf, dup= dellyDupDf, inv=dellyINVDf)
    attributes(retuRes) <- c(attributes(retuRes), list(method=method))
    
    return(retuRes)
}
