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
#'doClusterAnalysis
#'
#'Function to perform cluster analysis
#'@param clustPar List object that contains the parameters needed to carry out the analysis.
#'@examples
#'\dontrun{
#' clustParsList <- list()
#' for(i in 1:length(compName))
#' {
#'   clustPar <- list(expres = NULL,
#'                    expresFileName = "exprs.filtered.Rda",
#'                    geneListFName = paste("geneList",  compName[i],     ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"), "LT", pValCutOff[i], "Rda", sep = "."),
#'                    genes2cluster = NULL, 
#'                    samples2cluster = s2clust, 
#'                    sampleNames = as.character(targets$ShortName)[s2clust], 
#'                    comparisonName = compName[i], 
#'                    anotPackage = "org.Hs.eg",
#'                    my.symbols = symbolsTable,
#'                    outputDir = outputDir,
#'                    fileOfLinks = linksFile,
#'                    numClusters = 2,
#'                    rowDistance = NULL,
#'                    colDistance = NULL,
#'                    RowVals = TRUE,
#'                    ColVals = FALSE,
#'                    escala = "row",
#'                    colorsSet = pal,
#'                    densityInfo = "density",
#'                    colsForGroups = c("pink","pink","pink","pink","pink","blue","blue","blue","blue","blue"),
#'                    cexForColumns = 0.8,
#'                    cexForRows = 0.8,
#'                    Title = paste(compName[i],
#'                                  "with",
#'                                  ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
#'                                  "<",
#'                                  pValCutOff[i], ifelse(minLogFoldChange[i]==0, "", paste("\n and |logFC|>=", minLogFoldChange[i], sep=""))),
#'                    paste("Comparison:", compName[i], sep=" "),
#'                    csvType = csvType)
#'   
#'   clustParsList <- add2parsList(clustParsList, clustPar)
#' }
#' 
#' for(ix in 1:length(clustParsList))
#' {
#'  hm.Estudi <- BasicP::doClusterAnalysis(clustParsList[ix])
#' }
#' }
#'@export


doClusterAnalysis <- function(clustPar)
{

  p <- clustPar[[1]]
  if (!is.null(p$expresFileName)){
    expres <- loadFromFile (file.path(p$outputDir, p$expresFileName))
  }else{
    if (!is.null(p$dades)) {
      expres <- eval(parse(text = p$dades)) # Posar-hi un tryCatch per poder sortir si d??na error!!!
    }else{
      stop("Error, Cal definir o les dades o el nom de l'arxiu")
    }
  }

  if (!is.null(p$geneListFName)){
    genes2cluster <- loadFromFile (file.path(p$outputDir, p$geneListFName))
  }else{
    if (is.null(p$genes2cluster)) {
      stop("Error, Cal definir l'arxiu que conte la llista de gens o passar una variable que la contingui")
    }else{
      genes2cluster <- p$genes2cluster
    }
  }

  clust <- clusterAnalysis(expres = expres,
                           genes = genes2cluster,
                           samples = p$samples2cluster,
                           sampleNames = p$sampleNames,
                           comparisonName = p$comparisonName,
                           anotPackage = p$anotPackage,
                           my.symbols = p$my.symbols,
                           outputDir = p$outputDir,
                           fileOfLinks = p$fileOfLinks,
                           numClusters = p$numClusters,
                           rowDistance = p$rowDistance,
                           colDistance = p$colDistance,
                           RowVals = p$RowVals,
                           ColVals = p$ColVals,
                           escala = p$escala,
                           colorsSet = p$colorsSet,
                           densityInfo = p$densityInfo,
                           colsForGroups = p$colsForGroups,
                           cexForColumns = p$cexForColumns,
                           cexForRows = p$cexForRows,
                           Title = p$Title,
                           csvType = p$csvType)

  return(clust)
}
