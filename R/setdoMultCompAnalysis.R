sets_doMultCompAnalysis = 
  function(wCont, compName, adjMethod, pValCutOff, minLogFoldChange){
    for (i in 1:length(compName)){
      if(i==1){
        mcParsList = list()
      }
      mci = list(fitMain = get(load("./ResultsDir/fit.Rda")),
                 fitFileName = "fit.Rda",
                 whichContrast = wCont[[i]],
                 comparisonName = compName[i],
                 titleText = paste("for",ifelse(adjMethod[i]=="none",
                                                "p-values","adj. p-values"), "<", pValCutOff[i], "and |logFC| >",
                                   minLogFoldChange[i], sep = " "),
                 anotPackage = "org.Hs.eg",
                 my.symbols = get(load("./ResultsDir/Symbols.Rda")),
                 outputDir = "./ResultsDir",
                 fileOfLinks = "Links.txt",
                 multCompMethod = "separate",
                 adjustMethod = adjMethod[i],
                 selectionType = "any",
                 P.Value.cutoff = pValCutOff[i],
                 plotVenn = TRUE,
                 colsVenn = NULL,
                 vennColors = c("red","yellow","green","blue","pink"),
                 cexVenn = 1,
                 geneListFName = paste("geneList",compName[i],
                                       ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
                                       "LT",pValCutOff[i],"Rda",sep = "."),
                 minLogFC = minLogFoldChange[i],
                 csvType = "csv")
      mcParsList <- add2parsList(oneList=mcParsList, object=mci)
    }
    mcParsList<<-mcParsList
    return(mcParsList)
  }