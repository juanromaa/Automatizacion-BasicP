###########################################################
normalitza <- function (my.data, method)
{
  switch(method,
         RMA =  affy::rma(my.data),
         GCRMA = gcrma(my.data),
         MAS5 = mas5 (my.data))
}
###########################################################
#' normalization
#'
#' Function that normalizes raw data.
#'
#' @param my.data Raw data to be normalized.
#' @param method Method used for the normalization.
#' @param targetsinfo Target information about raw data.
#' @param inputDir Path of the normalized data file.
#' @param loadFile If TRUE the normalized data is loaded.
#' @param normalizedFName Name of the eset normalized file. By default "normalized.Rda"
#' @param outputDir Path to store the raw data normalized.
#' @param exonSt Default value FALSE. If TRUE the raw data is exon array.
#' @importFrom affy rma
#' @return A data set called my.norm with the raw data normalized.
#' @examples 
#' \dontrun{
#' load("rawData.Rda")
#' rawData <- my.raw
#' normMethod <- "RMA"
#' my.targets <- read.AnnotatedDataFrame("./celfiles/targets.txt", header = TRUE, row.names = 1)
#' celFilesDir <-"./celfiles"
#' loadFile <- FALSE 
#' normalized.eset.FileName <-  "normalizedData.Rda"   
#' outputDir <- "./ResultsDir"
#' exonStudy <- FALSE

#' eset_norm <- normalization(my.data = rawData, method = normMethod, 
#' targetsinfo = my.targets, inputDir = celFilesDir, loadFile = loadFile , 
#' normalizedFName = normalized.eset.FileName, outputDir = outputDir, 
#' exonSt = exonStudy)}
#' @export

normalization <- function(my.data = NULL,
                          method = NULL,
                          targetsinfo,
                          inputDir,
                          loadFile = FALSE,
                          normalizedFName = "normalized.Rda",
                          outputDir,
                          exonSt=FALSE)
{
  if (!loadFile)
  {
    if(exonSt)
    {
      my.norm <- affy::rma(my.data, target = ifelse(is.null(method), "core", method))
    }else{
        my.norm <- normalitza(my.data, method)
    }
    save(my.norm, file = file.path(outputDir, normalizedFName))
  }else{
    load(file = file.path(outputDir, normalizedFName))
  }
  return(my.norm)
}

