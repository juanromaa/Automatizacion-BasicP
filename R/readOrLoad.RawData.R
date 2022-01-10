#'readOrLoad.RawData
#'
#'Function that reads the raw data or if the eset is already created it is loaded.
#'
#' @param readCELS If TRUE it reads and stores de novo the raw data.
#' @param phenoDat Name of the file that contains phenodata information.
#' @param fileNames List of the file names to import.
#' @param dataFName Name of the file that stores raw data.
#' @param outputDir Path where the my.raw will be stored.
#' @param exonSt Default value FALSE. If TRUE the raw data is exon array.
#' @param cdf name of "cel definition file" package
#' @return my.raw expresion set object with the raw data.
#' @importFrom affy ReadAffy
#' @importFrom oligo read.celfiles
#' @examples
#' \dontrun{
#' require(BasicP)
#' require(affy)
#' require(oligo) #ExonStudy
#' readCELS <- TRUE
#' my.targets <- "./celfiles/targets.txt"
#' targets<- read.table("./celfiles/targets.txt", head=TRUE, sep="\t", row.names = 1)
#' my.fileNames <-paste("./celfiles/",rownames(targets),sep="")
#' rawDataFileName <- "rawData.Rda"
#' my.outputDir <- "."
#' isExonStudy <- FALSE
#' orgPackage <- "org.Hs.eg"
#'
#' rawData <- readOrLoad.RawData(readCELS = readCELS, phenoDat = my.targets,
#' fileNames = my.fileNames, dataFName =rawDataFileName,outputDir = my.outputDir,
#' exonSt = isExonStudy)}
#' @export

readOrLoad.RawData <- function (readCELS, phenoDat, fileNames, dataFName, outputDir, exonSt=FALSE, cdf=NULL)
{
  if (readCELS)
  {
    if(exonSt)
    {
      my.raw <- read.celfiles(filenames = fileNames, verbose = TRUE)
    }else{
      my.raw <- ReadAffy(filenames = fileNames, phenoData = phenoDat,
                         verbose = TRUE, cdfname=cdf)
    }
     save(my.raw, file = file.path(outputDir, dataFName))
  }else{
    load(file = file.path(outputDir, dataFName))
  }
  return(my.raw)
}

