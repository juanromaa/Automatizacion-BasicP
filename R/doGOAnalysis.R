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
extractSinonims <- function(my.strings)
{
  my.sinonims <- list()
  if (runMulticore ==1 || runMulticore ==3) {
    my.sinonims <- mclapply(my.strings, function(x) unlist(strsplit(x, " /// ")))

  } else {
    my.sinonims <- lapply(my.strings, function(x) unlist(strsplit(x, " /// ")))
  }

  return(my.sinonims)
}
###################################################
midSinonims <- function(my.IDs)
{
  my.indexes <- grep(" /// ", my.IDs)
  my.sinonimIDs <- my.IDs[my.indexes]

  my.IDs[my.indexes] <- extractSinonims(my.sinonimIDs)
  return(my.IDs)
}

###################################################

#'doGOAnalysis
#'
#'Function to interate with GOAnalysis
#'
#'@param GOParsList List object that contains the parameters needed to carry out the analysis.
#'@return GOResult
#'@exemples
#' \dontrun{
#' GOResult<-BasicP::GOAnalysis(fitMain = fitMain, whichContrasts = wCont, 
#' comparison.Name = "Estudi", outputDir = outputDir, anotPackage = "org.Hs.eg", 
#' my.IDs = entrezTable, addGeneNames = TRUE, fileOfLinks = linksFile, thrLogFC = 1, 
#' cutoffMethod = "adjusted", P.Value.cutoff = rep(0.05, length(wCont)), pval = 0.01,
#' min.count = 3, ontologias = c("MF", "BP", "CC"), testDirections = c("over", "under"), 
#' minNumGens = 0)
#'}
#'@export


doGOAnalysis <- function(GOPar)
{

  p <- GOPar[[1]]

  if(!is.null(p$my.IDs)){
    EntrezIDs <-  p$my.IDs #  Si se hace el eval solo guarda el ultimo id
  }

  if (!is.null(p$fitFileName))
  {
    fitMain <- loadFromFile(file.path(p$outputDir, p$fitFileName))
  }else{
    if (!is.null(p$fitMain))
    {
      fitMain <- p$fitMain # Posar-hi un tryCatch per poder sortir si d??na error!!!
    }else{
      stop("Error, Cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

  if (!is.null(p$my.IDs))
  {
    my.IDs <- midSinonims(p$my.IDs)
  }

  GOResult <- GOAnalysis(fitMain = fitMain,
                         whichContrasts = p$whichContrasts,
                         comparison.Name = p$comparisonName,
                         outputDir = p$outputDir,
                         anotPackage = orgPackage,
                         my.IDs = EntrezIDs,
                         addGeneNames = p$addGeneNames,
                         fileOfLinks = p$fileOfLinks,
                         thrLogFC = p$minLogFC,
                         cutoffMethod = p$cutoffMethod,
                         P.Value.cutoff = p$P.Value.cutoff,
                         pval = p$pvalGOterms,
                         min.count = p$min.count,
                         ontologias = p$ontologias,
                         testDirections = p$testDirections,
                         minNumGens = p$minNumGens)

  return(GOResult)
}
