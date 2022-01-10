sets_doclusterAnalysis = 
  function(compName, adjMethod, pValCutOff, minLogFoldChange){
    for(i in 1:length(compName))
    {
      if (i==1){
        clustParsList = list()
      }
      clustPar <- list(expres = NULL,
                       expresFileName = "expres.filtered.Rda",
                       geneListFName = paste("geneList",  compName[i], ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"), "LT", pValCutOff[i], "Rda", sep = "."),
                       genes2cluster = NULL, 
                       samples2cluster = tar_read(s2clust), 
                       sampleNames = as.character(tar_read(targets)$ShortName)[tar_read(s2clust)],
                       comparisonName = compName[i], 
                       anotPackage = "org.Hs.eg",
                       my.symbols = get(load("./ResultsDir/Symbols.Rda")),
                       outputDir = "./ResultsDir",
                       fileOfLinks = "Links.txt",
                       numClusters = 2,
                       rowDistance = NULL,
                       colDistance = NULL,
                       RowVals = TRUE,
                       ColVals = FALSE,
                       escala = "row",
                       colorsSet = colorpanel(n = 32, low = "green", mid = "white", high = "magenta"),
                       densityInfo = "density",
                       colsForGroups = c("pink","pink","pink","pink","pink","blue","blue","blue","blue","blue"),
                       cexForColumns = 0.8,
                       cexForRows = 0.8,
                       Title = paste(compName[i],
                                     "with",
                                     ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
                                     "<",
                                     pValCutOff[i], ifelse(minLogFoldChange[i]==0, "", paste("\n and |logFC|>=", minLogFoldChange[i], sep=""))),
                       paste("Comparison:", compName[i], sep=" "),
                       csvType = "csv2")
      
      clustParsList <- add2parsList(clustParsList, clustPar)
    }
    clustParsList <<- clustParsList
    return(clustParsList)
}