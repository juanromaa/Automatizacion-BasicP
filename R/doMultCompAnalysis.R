###################################################
loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, envir = tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}
###################################################
#' doMultCompAnalysis
#'
#' Function for make the comparative analysis.
#' @param mcPar List object that contains the parameters.
#' @examples 
#' \dontrun{
#' mcParsList <- list()
#' compName <- c("Group1", "Group2")
#' wCont <- 1:ncol(cont.matrix)
#' pValCutOff <- c(0.01, 0.01) 
#' adjMethod <- c("none", "none")
#' minLogFoldChange <- c(1, 1) 
#' for(i in 1:length(compName))
#' {  
#'   mci <- list(fitMain = fitMain,
#'               fitFileName = fitFileName,
#'               whichContrasts = wCont[[i]],
#'               comparisonName = compName[i],
#'               titleText = paste("for",ifelse(adjMethod[i]=="none",
#'                                              "p-values","adj. p-values"), "<", pValCutOff[i], "and |logFC| >",
#'                                 minLogFoldChange[i], sep = " "),
#'               anotPackage = "org.Hs.eg",
#'               my.symbols = symbolsTable,
#'               outputDir = outputDir,
#'               fileOfLinks = linksFile,
#'               multCompMethod = "separate",
#'               adjustMethod = adjMethod[i], 
#'               selectionType = "any",
#'               P.Value.cutoff = pValCutOff[i],
#'               plotVenn = TRUE,
#'               colsVenn = NULL,
#'               vennColors= c("red","yellow","green","blue","pink"),  
#'               cexVenn = 1,
#'               geneListFName = paste("geneList",compName[i],
#'                                     ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
#'                                     "LT",pValCutOff[i],"Rda",sep = "."),
#'               minLogFC = minLogFoldChange[i],
#'               csvType = csvType)
#'   
#'   mcParsList <- add2parsList(mcParsList, mci)
#' }
#' 
#' for(ix in 1:length(mcParsList))
#' {
#'   geneList.MCi <- BasicP::doMultCompAnalysis(mcParsList[ix])   
#' }
#' }
#' @export


doMultCompAnalysis <- function(mcPar)
{

  p <- mcPar[[1]]
  if (!is.null(p$fitFileName)){
    fitMain <- loadFromFile(file.path(p$outputDir,p$fitFileName))
  }else{
    if (!is.null(p$fitMain)) {
      fitMain <- eval(parse(text = p$fitMain))
      # Posar-hi un tryCatch per poder sortir si dona error!!!
    }else{
      stop("Error, Cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

  geneList <-  multipleComp(fitMain = fitMain,
                            whichContrasts = p$whichContrast,
                            comparisonName = p$comparisonName,
                            titleText = p$titleText,
                            outputDir = p$outputDir,
                            anotPackage = p$anotPackage,
                            my.symbols = p$my.symbols,
                            linksFile = p$fileOfLinks,
                            multCompMethod = p$multCompMethod,
                            adjustMethod = p$adjustMethod,
                            selectionType = p$selectionType,
                            P.Value.cutoff = p$P.Value.cutoff,
                            plotVenn = p$plotVenn,
                            colsVenn = p$colsVenn,
                            vennColors = p$vennColors,
                            cexVenn = p$cexVenn,
                            geneListFName=p$geneListFName,
                            csvType=p$csvType,
                            minLFC=p$minLogFC)

  return (geneList)
}
