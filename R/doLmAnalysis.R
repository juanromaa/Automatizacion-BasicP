##############################################################
loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}
##############################################################
write2csv <- function(my.data, fileName, csv = c("csv2", "csv", "txt", "xls"), outputDir)
{
  fileName<- file.path(outputDir, paste(fileName, substr(csv[1], 1, 3) , sep = "."))
  switch (csv[1],
          "csv" = write.csv(my.data, file = fileName, quote = F),
          "csv2" = write.csv2(my.data, file = fileName, quote = F),
          "txt" = write.table(my.data, file = fileName, quote = F))
}

##############################################################
#' add2parsList
#' @param oneList List in which the new parameter will be introduced.
#' @param object Object to add in the list.
#' @return Returns the list with the new parameter.
#' @export



add2parsList <- function(oneList,object)
{
  pos <- length(oneList) + 1
  oneList[[pos]] <- object
  names(oneList)[pos] <- oneList[[pos]]$comparisonName
  return(oneList)
}
##############################################################

#' doLmAnalysis
#' 
#' Function to make the Linear model analysis. The unique parameter is lmPar list.
#' @param lmPar list object that contains the parameters needed to carry out the analysis.
#' @importFrom links2File addToLinksFile
#' @exemples
#' \dontrun{
#' lmParsList <- list()
#' Estudi <- list(dades = NULL,
#'                expresFileName = "exprs.filtered.Rda",
#'                targets = targets,
#'                designMat = design,
#'                contMat = cont.matrix,
#'                whichContrasts = 1:ncol(cont.matrix),
#'                anotPack = NULL,
#'                outputDir = outputDir,
#'                ExpressionsAndTop = TRUE,
#'                showLmParams = FALSE, 
#'                use.dupCorr = FALSE,
#'                block = NULL,
#'                nDups = 1,
#'                comparisonName = comparison,  
#'                ENTREZIDs = "entrezTable",
#'                SYMBOLIDs = "symbolsTable",
#'                fileOfLinks = linksFile,
#'                fitFileName = fitFileName,
#'                csvType=csvType,
#'                rows2HTML = NULL,
#'                anotFilename = anotFileName
#' )
#' 
#' lmParsList <- add2parsList(lmParsList, Estudi)
#' 
#' for(ix in 1:length(lmParsList))
#' {
#'   fit.Main   <- BasicP::doLmAnalysis(lmParsList[ix])
#' }
#' }
#' 
#' @export

doLmAnalysis <- function (lmPar)
{
  p <- lmPar[[1]]
  if(!is.null(p$ENTREZIDs)){
    EntrezIDs <-  p$ENTREZIDs
  }

  if(!is.null(p$SYMBOLIDs)){
    SymbolIDs <-  p$SYMBOLIDs
  }

  if (!is.null(p$expresFileName)){
    expres <- loadFromFile (file.path(p$outputDir, p$expresFileName))
  }else{
    if (!is.null(p$dades))
    {
      expres <- eval(parse(text = p$dades)) # Posar-hi un tryCatch per poder sortir si d error!!!
    }else{
      stop("Error, Cal definir o les dades o el nom de l'arxiu")
    }
  }

  if (is.null(p$whichContrasts))
  {
    contrasts2test <- 1:ncol(p$contMat)
  }else{
    contrasts2test <- p$whichContrasts
  }

  fitMain <- lmAnalysis(exprs.filtered = expres,
                        design = p$designMat,
                        cont.matrix = p$contMat,
                        contrasts2test = contrasts2test,
                        anotPackage = p$anotPack,
                        outputDir = p$outputDir,
                        comparison = p$comparisonName,
                        Expressions_And_Top = p$ExpressionsAndTop,
                        showParams = p$showLmParams,
                        use.dupCorr = p$use.dupCorr,
                        block = p$block,
                        nDups = p$nDups ,
                        ENTREZIDs = EntrezIDs,
                        SYMBOLIDs = SymbolIDs,
                        linksFile = p$fileOfLinks,
                        fitFileName = p$fitFileName,
                        csvType=p$csvType,
                        rows2HTML = p$rows2HTML,
                        anotFileName = p$anotFilename
  )


  designMatrixName = paste("designMatrix",p$comparisonName, sep=".")
  contrastMatrixName = paste("contrastMatrix",p$comparisonName, sep=".")

  write2csv(p$designMat, fileName = designMatrixName, csv = p$csvType, outputDir = p$outputDir)
  write2csv(p$contMat, fileName = contrastMatrixName, csv = p$csvType, outputDir = p$outputDir)

  csvType <- ifelse(is.null(p$csvType), "csv2", p$csvType)
  addToLinksFile(p$fileOfLinks, paste(designMatrixName, substr(csvType, 1, 3),  sep="."), categ = "ANALYSIS",
                 desc = paste("Design Matrix for comparison ", p$comparisonName, sep=""))
  addToLinksFile(p$fileOfLinks, paste(contrastMatrixName, substr(csvType, 1, 3),  sep="."), categ = "ANALYSIS",
                 desc = paste("Contrast Matrix for comparison ", p$comparisonName, sep=""))

  return (fitMain)
}
