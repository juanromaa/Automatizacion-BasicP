############################################################

loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}
############################################################
#'doKEGGAnalysis
#'
#' Function to interate with KEGGAnalysis
#' @param KEGGPar List object that contains the parameters needed to carry out the analysis.
#' @return KEGGResult
#' @examples  
#' \dontrun{
#' KEGGParsList <- list()
#' KEGGPar <- list(fitFileName = "fit.Rda",
#'                whichContrasts = 1:3,
#'                 comparisonName = "Estudi",
#'                 outputDir = "./ResultsDir",
#'                 anotPackage = "org.Hs.eg", 
#'                 organisme = "mmu",    
#'                 my.IDs = "entrezTable",
#'                 addGeneNames = TRUE,
#'                 fileOfLinks = linksFileName,
#'                 fitMain = NULL,
#'                 cutoffMethod = "unadjusted",
#'                 P.Value.cutoff = rep(0.01, length(wCont)), 
#'                 pvalKEGGterms = 0.05, 
#'                 minNumGens = 0) 
#' 
#' KEGGParsList <- add2parsList (KEGGParsList, KEGGPar)
#' 
#' 
#' for(i in 1:length(KEGGParsList))
#' {
#'   KEGGList <- BasicP::doKEGGAnalysis(KEGGParsList[i])
#' }
#' }
#' @export
#'
doKEGGAnalysis <- function(KEGGPar)
{

  p <- KEGGPar[[1]]
  if(!is.null(p$my.IDs)){
    EntrezIDs <-   p$my.IDs #  Seran el EntrezTable
  }

  if (!is.null(p$fitFileName))
  {
    fitMain <- loadFromFile(file.path(p$outputDir, p$fitFileName))
  }else{
    if (!is.null(p$fitMain))
    {
      fitMain <- eval(parse(text = p$fitMain)) # Posar-hi un tryCatch per poder sortir si dona error!!!
    }else{
      stop("Error, cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

  KEGGResult <- KEGGAnalysis(fitMain = fitMain,
                             whichContrasts = p$whichContrasts,
                             comparison.Name = p$comparisonName,
                             outputDir = p$outputDir,
                             anotPackage = p$anotPackage,
                             organisme = p$organisme,
                             my.IDs = EntrezIDs, # era p$my.IDs
                             addGeneNames = p$addGeneNames,
                             fileOfLinks = p$fileOfLinks,
                             cutoffMethod = p$cutoffMethod,
                             P.Value.cutoff = p$P.Value.cutoff,
                             pval = p$pvalKEGGterms,
                             thrLogFC = p$minLogFC,
                             minNumGens = p$minNumGens)

  return(KEGGResult)
}
