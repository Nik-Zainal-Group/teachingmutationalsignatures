
#' @export
fetchData <- function(datasetname){
  if(is.null(dataObj[[datasetname]])){
    message("[error fetchData] datatset name ",datasetname," unknown. Available datasets are: ",
            paste(names(dataObj),collapse = ","),". Nothing returned.")
    return(NULL)
  }else{
    return(dataObj[[datasetname]]$catalogues)
  }
  
}

#' @export
checkPerformanceSignatures <- function(estimated_signatures,
                                       datasetname,
                                       outfile = NULL){
  if(is.null(dataObj[[datasetname]])){
    message("[error fetchData] datatset name ",datasetname," unknown. Available datasets are: ",
            paste(names(dataObj),collapse = ","),". Nothing done.")
    return(NULL)
  }else{
    resPerf <- signature.tools.lib::evaluatePerformanceSignatureSimilarity(true_signatures = dataObj[[datasetname]]$signatures,
                                                                           estimated_signatures = estimated_signatures,
                                                                           true_exposures = t(dataObj[[datasetname]]$exposures),
                                                                           outfile = outfile)
    return(resPerf)
  }
}

#' @export
checkPerformanceSignaturesList <- function(estimated_signatures_list,
                                           datasetname,
                                           outfile = NULL){
  if(is.null(dataObj[[datasetname]])){
    message("[error fetchData] datatset name ",datasetname," unknown. Available datasets are: ",
            paste(names(dataObj),collapse = ","),". Nothing done.")
    return(NULL)
  }else{
    resPerf <- signature.tools.lib::evaluatePerformanceSignatureSimilarityList(true_signatures = dataObj[[datasetname]]$signatures,
                                                                               estimated_signatures_list = estimated_signatures_list,
                                                                               true_exposures = t(dataObj[[datasetname]]$exposures),
                                                                               outfile = outfile)
    return(resPerf)
  }
}


#' @export
checkPerformanceExposures <- function(estimated_exposures,
                                      datasetname,
                                      outfile = NULL){
  if(is.null(dataObj[[datasetname]])){
    message("[error fetchData] datatset name ",datasetname," unknown. Available datasets are: ",
            paste(names(dataObj),collapse = ","),". Nothing done.")
    return(NULL)
  }else{
    # check if the signatures and samples match, and possibly transpose
    transposed <- FALSE
    okmatch <- FALSE
    quitloop <- FALSE
    while(!quitloop){
      checkRowNames <- length(intersect(rownames(dataObj[[datasetname]]$exposures),rownames(estimated_exposures)))>0
      checkColNames <- length(intersect(colnames(dataObj[[datasetname]]$exposures),colnames(estimated_exposures)))>0
      okmatch <- checkRowNames & checkColNames
      if(okmatch){
        quitloop <- TRUE
      }else if(!okmatch & !transposed){
        estimated_exposures <- t(estimated_exposures)
        transposed <- TRUE
      }else if(!okmatch & transposed){
        quitloop <- TRUE
      }
    }
    
    # let's remove the unassigned column if present
    if ("unassigned" %in% rownames(estimated_exposures)){
      pos <- which("unassigned" == rownames(estimated_exposures))
      estimated_exposures <- estimated_exposures[-pos,,drop=F]
    }
    
    # let's check if there are common/rare signatures and we need to compute
    # separate performance metrics
    commonSigsNames <- colnames(dataObj[[datasetname]]$signatures)[grepl(colnames(dataObj[[datasetname]]$signatures),pattern = "common")]
    commonSigsNames <- union(commonSigsNames,rownames(estimated_exposures)[grepl(rownames(estimated_exposures),pattern = "common")])
    rareSigsNames <- colnames(dataObj[[datasetname]]$signatures)[!grepl(colnames(dataObj[[datasetname]]$signatures),pattern = "common")]
    rareSigsNames <- union(rareSigsNames,rownames(estimated_exposures)[!grepl(rownames(estimated_exposures),pattern = "common")])
    
    if(length(commonSigsNames)==0) commonSigsNames <- NULL
    if(length(rareSigsNames)==0) rareSigsNames <- NULL
    
    resPerf <- signature.tools.lib::evaluatePerformanceExposures(true_exposures = t(dataObj[[datasetname]]$exposures),
                                                                 estimated_exposures = t(estimated_exposures),
                                                                 commonNames = commonSigsNames,
                                                                 rareNames = rareSigsNames,
                                                                 outfile = outfile)
    return(resPerf)
  }
}

#' @export
checkPerformanceExposuresList <- function(estimated_exposures_list,
                                          datasetname,
                                          outfile = NULL){
  if(is.null(dataObj[[datasetname]])){
    message("[error fetchData] datatset name ",datasetname," unknown. Available datasets are: ",
            paste(names(dataObj),collapse = ","),". Nothing done.")
    return(NULL)
  }else{
    
    # check and update the estimated exposures list
    for (g in names(estimated_exposures_list)){
      estimated_exposures <- estimated_exposures_list[[g]]
      # check if the signatures and samples match, and possibly transpose
      transposed <- FALSE
      okmatch <- FALSE
      quitloop <- FALSE
      while(!quitloop){
        checkRowNames <- length(intersect(rownames(dataObj[[datasetname]]$exposures),rownames(estimated_exposures)))>0
        checkColNames <- length(intersect(colnames(dataObj[[datasetname]]$exposures),colnames(estimated_exposures)))>0
        okmatch <- checkRowNames & checkColNames
        if(okmatch){
          quitloop <- TRUE
        }else if(!okmatch & !transposed){
          estimated_exposures <- t(estimated_exposures)
          transposed <- TRUE
        }else if(!okmatch & transposed){
          quitloop <- TRUE
        }
      }
      
      # let's remove the unassigned column if present
      if ("unassigned" %in% rownames(estimated_exposures)){
        pos <- which("unassigned" == rownames(estimated_exposures))
        estimated_exposures <- estimated_exposures[-pos,,drop=F]
      }
      
      # update
      estimated_exposures_list[[g]] <- t(estimated_exposures)
    }
    
    # let's check if there are common/rare signatures and we need to compute
    # separate performance metrics
    commonSigsNames <- colnames(dataObj[[datasetname]]$signatures)[grepl(colnames(dataObj[[datasetname]]$signatures),pattern = "common")]
    for(g in names(estimated_exposures_list)) commonSigsNames <- union(commonSigsNames,colnames(estimated_exposures_list[[g]])[grepl(colnames(estimated_exposures_list[[g]]),pattern = "common")])
    rareSigsNames <- colnames(dataObj[[datasetname]]$signatures)[!grepl(colnames(dataObj[[datasetname]]$signatures),pattern = "common")]
    for(g in names(estimated_exposures_list)) rareSigsNames <- union(rareSigsNames,colnames(estimated_exposures_list[[g]])[!grepl(colnames(estimated_exposures_list[[g]]),pattern = "common")])
    
    if(length(commonSigsNames)==0) commonSigsNames <- NULL
    if(length(rareSigsNames)==0) rareSigsNames <- NULL
    
    resPerf <- signature.tools.lib::evaluatePerformanceExposuresList(true_exposures = t(dataObj[[datasetname]]$exposures),
                                                                     estimated_exposures_list = estimated_exposures_list,
                                                                     commonNames = commonSigsNames,
                                                                     rareNames = rareSigsNames,
                                                                     outfile = outfile)
    return(resPerf)
  }
}



