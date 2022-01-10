#Carga de targets y establecer directorio principal
library(targets)
setwd("")

#Carga de archivos
source("R/readOrLoad.rawData.R")
source("R/createOrLoadAnnotations.R")
source("R/normalization.R")
source("R/normplots2File.R")
source("R/filterData.R")
source("R/saveData.R")
source("R/lmAnalysis.R")
source("R/doLmAnalysis.R")
source("R/GeneAnnotation.R")
source("R/doGeneAnnotation.R")
source("R/multipleComp.R")
source("R/doMultCompAnalysis.R")
source("R/setdoMultCompAnalysis.R")
source("R/clusterAnalysis.R")
source("R/doClusterAnalysis.R")
source("R/setdoClusterAnalysis.R")
source("R/GoAnalysis.R")
source("R/doGoAnalysis.R")
source("R/KEGGAnalysis.R")
source("R/doKEGGAnalysis.R")

#Valores predeterminados
runMulticore = 3
toTIFF = FALSE

#Paquetes a utilizar
tar_option_set(
  packages = c(
    "affy",
    "oligo",
    "limma",
    "annotate",
    "BiocManager",
    "Biobase",
    "BiocGenerics",
    "parallel",
    "GOstats",
    "SortableHTMLTables",
    "org.Hs.eg.db",
    "VennDiagram",
    "grDevices",
    "grid",
    "stats",
    "gplots",
    "DBI"
  ))

#listado de targets
list(
  #ReadOrLoad
  tar_target(readOrLoad.rawData.pars, 
             list(readCELS = TRUE, 
                  phenoDat = "./celfiles/targets.txt",
                  fileNames = paste("./celfiles/", 
                                    rownames(read.table("./celfiles/targets.txt", 
                                                        head=TRUE, sep="\t",  
                                                        row.names = 1)), sep=""),
                  dataFName = "rawData.Rda",
                  outputDir = ".",
                  exonSt = FALSE),
             priority = 1),
  
  tar_target(readOrLoad.rawData.func, 
             do.call(readOrLoad.RawData, readOrLoad.rawData.pars),
             priority = 0.99),
  
  #outputDir 
  tar_target(outputDir, "./ResultsDir", priority = 0.98),
  
  
  #CreateOrLoad
  tar_target(createOrLoadAnnotations.pars, list(loadAnnotations = FALSE,
                                                chipPackAvailable = TRUE,
                                                platformDesignPackAvailable = NULL,
                                                chipPackage = "hgu133a2",
                                                platformDesignPackage = NULL,
                                                outputDir = outputDir,
                                                annotationsFileName = "Annotations",
                                                entrezTableFileName = "Entrezs.Rda",
                                                symbolsTableFileName = "Symbols.Rda",
                                                controlsTableFileName = "controls.Rda"),
             priority = 0.97),
  
  tar_target(createOrLoadAnnotations.func, 
             do.call(createOrLoadAnnotations, createOrLoadAnnotations.pars),
             priority = 0.96),
  
  
  
  #Normalization
  tar_target(normalization.pars, 
             list(my.data = get(load("./rawData.Rda")),
                  method = "RMA",
                  targetsinfo = read.AnnotatedDataFrame("./celfiles/targets.txt",
                                                        header = TRUE, 
                                                        row.names = 1),
                  inputDir = "./celfiles",
                  loadFile = FALSE,
                  normalizedFName = "normalizedData.Rda",
                  outputDir = outputDir,
                  exonSt = FALSE),
             priority = 0.95),
  
  tar_target(normalization.func, 
             do.call(normalization, normalization.pars),
             priority = 0.94),
  
  #Normalized
  tar_target(normalized, 
             get(load("./ResultsDir/normalizedData.Rda")), 
             priority = 0.93),
  
  #Normplots2File
  tar_target(normplots2File.pars, 
             list(my.data = normalized,
                  sampleNames = pData(normalized)$ShortName,
                  my.colors = rainbow(length(sampleNames(normalized))),
                  my.groups = pData(normalized)$Group,
                  my.method = "average",
                  my.cex = 0.8,
                  posText = 2,
                  dim3 = FALSE,
                  fileName = "NormalizedPlots.pdf",
                  outputDir = outputDir,
                  PCAPlots=TRUE,
                  csv = "csv2",
                  lFile = NULL),
             priority = 0.92),
  
  tar_target(normplots2File.func, 
             do.call(normplots2File, normplots2File.pars),
             priority = 0.91),

  #Entrez, controls y targets
  tar_target(entrez, get(load("./ResultsDir/Entrezs.Rda")), priority = 0.90),
  
  tar_target(controls, get(load("./ResultsDir/controls.Rda")), priority = 0.89),
  
  tar_target(targets,read.AnnotatedDataFrame("./celfiles/targets.txt", 
                                             header = TRUE,row.names = 1), 
             priority = 0.88),
  
  #FilterData
  tar_target(filterData.pars, 
             list(expres = exprs(normalized)[!duplicated(exprs(normalized), MARGIN=1),],
                  controls = names(controls),
                  removeNAs = TRUE,
                  entrezs = entrez,
                  bySignal = TRUE,
                  signalThr = 50,
                  grups = pData(normalized)$Group, 
                  sigFun.Name = "filt.by.Signal", 
                  sigThr.as.perc = TRUE,
                  byVar = TRUE,
                  variabilityThr = 50,
                  varFun.Name = "sdf",
                  varThr.as.perc = TRUE,
                  pairingFun.Name = NULL,
                  targets = targets,
                  doReport = TRUE,          
                  outputDir = outputDir,                
                  filteringReportFName = "FilteringReport.txt"),
             priority = 0.87),
  
  tar_target(filterData.func, 
             do.call(filterData, filterData.pars),
             priority = 0.86),

  #Symbols y linkfile
  tar_target(symbols, get(load("./ResultsDir/Symbols.Rda")), priority = 0.85),
  tar_target(linkfile, "Links.txt", priority = 0.845),
  
  #SaveData de todos los genes
  tar_target(saveData.pars, 
             list(expres = exprs(normalized)[!duplicated(exprs(normalized),MARGIN=1),], 
                  expresNames = colnames(normalized),
                  expres.csv.FileName = "normalized.all", 
                  csvType = "csv2", 
                  description = "Normalized Values for all genes",
                  anotPackage = NULL, 
                  SYMBOL = "SYMBOL",
                  symbolsVector = symbols,
                  expres.bin.FileName = "expres.Rda",
                  linksFile = linkfile, 
                  outputDir = outputDir),
             priority = 0.84),
  
  tar_target(saveData.func, 
             do.call(saveData, saveData.pars),
             priority = 0.83),
  
  #saveData de genes filtrados
  tar_target(saveData.filt.pars, 
             list(expres = filterData.func, 
                  expresNames = colnames(filterData.func),
                  expres.csv.FileName = "normalized.filtered", 
                  csvType = "csv2", 
                  description="Normalized Values for filtered genes",
                  anotPackage = NULL, 
                  SYMBOL = "SYMBOL",
                  symbolsVector = symbols,
                  expres.bin.FileName = "expres.filtered.Rda",
                  linksFile = linkfile, 
                  outputDir = outputDir),
             priority = 0.82),
  
  tar_target(saveData.filt.func, 
             do.call(saveData, saveData.filt.pars),
             priority = 0.81),

  #designmatrix
  tar_target(designmatrix, function(datos, columnas, filas){
    design = model.matrix(~ 0 + datos)
    colnames(design) = columnas
    rownames(design) = filas
    return(design)
  },
  priority = 0.80),
  
  #expres.filtered y targets.designmatrix
  tar_target(expres.filtered, get(load("./ResultsDir/expres.filtered.Rda")), 
             priority = 0.79),
  tar_target(targets.designmatrix, read.table("./celfiles/targets.txt",
                                              head=TRUE, sep="\t"), priority = 0.78),
  
  #lmAnalysis
  tar_target(lmAnalysis.pars, 
             list(exprs.filtered = expres.filtered,
                  
                  design = designmatrix(datos = targets.designmatrix[,5], 
                                        columnas = unique(targets.designmatrix[,5]), 
                                        filas = targets.designmatrix$ShortName),
                  
                  cont.matrix = makeContrasts(AvsB= (A-B), AvsL= (A-L), BvsL= (B-L),
                                              levels = designmatrix(
                                                datos = targets.designmatrix[,5],
                                                columnas = unique(targets.designmatrix[,5]),
                                                filas = targets.designmatrix$ShortName)),
                  
                  contrasts2test = 1:ncol(makeContrasts(
                    AvsB= (A-B), AvsL= (A-L), BvsL= (B-L),
                    levels = unique(read.table("./celfiles/targets.txt", 
                                               head=TRUE, sep="\t",row.names = 1)[,4]))),
                  anotPackage = NULL,
                  outputDir = outputDir,
                  comparison = "Estudi",
                  Expressions_And_Top = TRUE,
                  showParams = FALSE,
                  use.dupCorr = FALSE,
                  block = NULL, 
                  nDups = 1,
                  ENTREZIDs = entrez,
                  SYMBOLIDs = symbols,
                  linksFile = linkfile,
                  fitFileName = "fit.Rda",
                  csvType = "csv",
                  rows2HTML = NULL,
                  anotFileName = "Annotations"),
             priority = 0.77),
  
  tar_target(lmAnalysis.func, 
             do.call(lmAnalysis, lmAnalysis.pars),
             priority = 0.76),

  #doLmAnalysis
  tar_target(doLmAnalysis.pars, 
             add2parsList(lmParslist<-list(),
                          Estudi <- 
                            list(dades = NULL,
                                 expresFileName = "expres.filtered.Rda",
                                 targets = targets.designmatrix,
                                 designMat = lmAnalysis.pars$design,
                                 contMat = lmAnalysis.pars$cont.matrix,
                                 whichContrasts = 1:ncol(lmAnalysis.pars$cont.matrix),
                                 anotPack = NULL,
                                 outputDir = outputDir,
                                 ExpressionsAndTop = TRUE,
                                 showLmParams = FALSE, 
                                 use.dupCorr = FALSE,
                                 block = NULL,
                                 nDups = 1,
                                 comparisonName = "Estudi",
                                 ENTREZIDs = entrez,
                                 SYMBOLIDs = symbols, 
                                 fileOfLinks = linkfile,
                                 fitFileName = "fit.Rda",
                                 csvType= "csv",
                                 rows2HTML = NULL,
                                 anotFilename = "Annotations")),
             priority = 0.75),
  
  tar_target(doLmAnalysis.func, 
             for(ix in 1:length(doLmAnalysis.pars)){
               fit.Main <- doLmAnalysis(doLmAnalysis.pars[ix])
             },
             priority = 0.74),
 
  #fitmain
  tar_target(fitmain, get(load("./ResultsDir/fit.Rda")), priority = 0.73),
  
  #GeneAnnotation
  tar_target(GeneAnnotation.pars, list(egIDs = entrez[unique(rownames(fitmain$p.value))],
                                       anotPackage = "org.Hs.eg",
                                       toHTML = TRUE,
                                       outputDir = outputDir,
                                       filename = "Annotations",
                                       myTitle = "Annotations for filtered genes",
                                       specie = "Homo_sapiens",
                                       info2show = c("Affymetrix", "EntrezGene",
                                                     "GeneSymbol", "GeneName", 
                                                     "KEGG", "GO"),
                                       linksFile = linkfile,
                                       maxGenes = NULL),
             priority = 0.72),
  
  tar_target(GeneAnnotation.func, 
             do.call(GeneAnnotation,GeneAnnotation.pars),
             priority = 0.71),
  
  #doGeneAnnotation
  tar_target(doGeneAnnotation.pars, list(list(fitMain = NULL,
                                              fitFileName = "fit.Rda",
                                              my.IDs = entrez,
                                              anotPackage = "org.Hs.eg",
                                              toHTML = TRUE,
                                              outputDir = outputDir,
                                              anotFilename = "Annotations",
                                              titleAnotations = "Annotations for filtered genes",
                                              specie = "homo sapiens",
                                              info2show = c("Affymetrix", "EntrezGene", 
                                                            "GeneSymbol", "GeneName", 
                                                            "KEGG", "GO"),
                                              linksFile = linkfile,
                                              numGenesPerPage = NULL,
                                              organisme = "hsa")),
             priority = 0.70),
  
  tar_target(doGeneAnnotations.func, 
             do.call(doGeneAnnotation,doGeneAnnotation.pars),
             priority = 0.69),

  #multipleComp.var
  tar_target(multipleComp.var, list(adjMethod = c("none"),
                                    pValcutoff = 0.01,
                                    compName = c("Group1"),
                                    minLogFoldChange = 1), 
             priority = 0.68),
  
  #multiplecomp
  tar_target(multipleComp.pars, 
             list(fitMain = fitmain,
                  whichContrasts = 1:3,
                  comparisonName = multipleComp.var$compName,
                  titleText = paste("for", 
                                    ifelse(multipleComp.var$adjMethod==
                                             "none","p-values","adj. p-values"),
                                    "<", 
                                    multipleComp.var$pValcutoff, 
                                    "and |logFC| >", 
                                    multipleComp.var$minLogFoldChange, sep = " "),
                  outputDir = outputDir,
                  anotPackage = "org.Hs.eg",
                  my.symbols = symbols,
                  linksFile = linkfile,
                  multCompMethod = "separate",
                  adjustMethod = c("none"),
                  selectionType = "any",
                  P.Value.cutoff = 0.01,
                  plotVenn = TRUE,
                  colsVenn = NULL,
                  vennColors = c("red","yellow","green","blue","pink"),
                  cexVenn = 1,
                  geneListFName = paste("geneList",
                                        multipleComp.var$compName,
                                        ifelse(multipleComp.var$adjMethod==
                                                 "none","pvalues","adj-pvalues"),
                                        "LT",
                                        multipleComp.var$pValcutoff,
                                        "Rda",sep = "."),
                  csvType = "csv",
                  minLFC = 1),
             priority = 0.67),
  
  tar_target(multipleComp.func, 
             do.call(multipleComp, multipleComp.pars),
             priority = 0.66),

  #doMultiCompAnalysis.var
  tar_target(doMultCompAnalysis.var, list(wCont = (1:ncol(lmAnalysis.pars$cont.matrix)),
                                          compName = c("Group1", "Group2"),
                                          adjMethod = c("none", "none"),
                                          pValCutOff = c(0.01, 0.01),
                                          minLogFoldChange = c(1, 1)),
             priority = 0.65),
  
  #doMultCompAnalysis
  tar_target(doMultCompAnalysis.pars, do.call(sets_doMultCompAnalysis,doMultCompAnalysis.var),
             priority = 0.64),
  
  tar_target(doMultCompAnalysis.func,
             for (ix in 1:length(doMultCompAnalysis.pars))
             {
               geneList.MCi = doMultCompAnalysis(doMultCompAnalysis.pars[ix])
             }, 
             priority = 0.63),

  #s2clust
  tar_target(s2clust, 
             which(as.logical(apply(lmAnalysis.pars$design[
               ,
               as.logical(apply(abs(as.matrix(lmAnalysis.pars$cont.matrix[
                 ,
                 (1:ncol(lmAnalysis.pars$cont.matrix))[[1]]])),1,sum))],1,sum))),
             priority = 0.626),
  
  #clusterAnalysis
  tar_target(clusterAnalysis.pars, 
             list(expres = expres.filtered,
                  genes = get(load("./ResultsDir/geneList.Group1.pvalues.LT.0.01.Rda")),
                  samples = s2clust,
                  sampleNames = as.character(targets$ShortName)[s2clust],
                  comparisonName = "Compar 1",
                  anotPackage = "org.Hs.eg",
                  my.symbols = symbols,
                  outputDir = outputDir,
                  fileOfLinks = linkfile,
                  numCluster = 2,
                  rowDistance = NULL,
                  colDistance = NULL,
                  RowVals = TRUE,
                  ColVals = TRUE,
                  escala = "row",
                  colorsSet = colorpanel(n = 32, 
                                         low = "green", 
                                         mid = "white", 
                                         high = "magenta"),
                  densityInfo = "denisty",
                  colsForGroups = c("pink","pink","pink","pink","pink",
                                    "blue","blue","blue","blue","blue"),
                  cexForColumns = 0.8,
                  cexForRows = 0.8,
                  Title = paste("Compar 1 with",
                                ifelse(c("none") == "none", "pvalues", "adj-pvalues"),
                                "<",
                                c(0.01), ifelse(1==0, "", 
                                                paste("\n and |logFC|>=", 
                                                      1, sep=""))),
                  csvType = "csv2"),
             priority = 0.625),
  
  tar_target(clusterAnalysis.func, do.call(clusterAnalysis, clusterAnalysis.pars),
             priority = 0.623),

  #doclusterAnalysis.var
  tar_target(doClusterAnalysis.var, list(compName = c("Group1", "Group2"),
                                         adjMethod = c("none", "none"),
                                         pValCutOff = c(0.01, 0.01),
                                         minLogFoldChange = c(1, 1)),
             priority = 0.625),
  
  #doClusterAnalysis
  tar_target(doClusterAnalysis.pars, do.call(sets_doclusterAnalysis, doClusterAnalysis.var),
             priority = 0.62),
  
  tar_target(doclusterAnalysis.func,
             for (ix in 1:length(doClusterAnalysis.pars))
             {
               hm.Estudi <- doClusterAnalysis(doClusterAnalysis.pars[ix])
             },
             priority = 0.61),

  #GoAnalysis  
  tar_target(GoAnalysis.pars, list(fitMain = fitmain,
                                   whichContrasts = 1:ncol(lmAnalysis.pars$cont.matrix),
                                   comparison.Name = "Estudi",
                                   outputDir = outputDir,
                                   anotPackage = "org.Hs.eg",
                                   my.IDs = entrez,
                                   addGeneNames = TRUE,
                                   fileOfLinks = linkfile,
                                   thrLogFC = 1,
                                   cutoffMethod = "adjusted",
                                   P.Value.cutoff = rep(0.05, 
                                                        length(1:ncol(
                                                          lmAnalysis.pars$cont.matrix))),
                                   pval = 0.01,
                                   min.count = 3, 
                                   ontologias = c("MF", "BP", "CC"),
                                   testDirections = c("over", "under"),
                                   minNumGens = 0),
             priority = 0.60),
  
  tar_target(GOAnalysis.func,  do.call(GOAnalysis, GoAnalysis.pars),
             priority = 0.59),
  
  #doGoAnalysis
  tar_target(doGOAnalysis.pars,
             add2parsList(GOParsList <- list(),
                          GOPar <- list(fitMain = fitmain,
                                        whichContrasts = 1:ncol(lmAnalysis.pars$cont.matrix),
                                        comparison.Name = "Estudi",
                                        outputDir = outputDir,
                                        anotPackage = "org.Hs.eg",
                                        my.IDs = entrez,
                                        addGeneNames = TRUE,
                                        fileOfLinks = linkfile,
                                        thrLogFC = 1,
                                        cutoffMethod = "adjusted",
                                        P.Values.cutoff = rep(0.05, 
                                                              length(1:ncol(
                                                                lmAnalysis.pars$cont.matrix))),
                                        pval = 0.01,
                                        min.count = 3,
                                        ontologias = c("MF", "BP", "CC"),
                                        testDirections = c("over", "under"),
                                        minNumGens = 0)),
             priority = 0.58),
  tar_target(doGOAnalysis.func, for (i in 1:length(doGOAnalysis.pars))
  {
    GOList = doGOAnalysis(doGOAnalysis.pars[i])
  },
  priority = 0.57),

  #KEGGAnalysis
  tar_target(KEGGAnalysis.pars, list(fitMain = fitmain,
                                     whichContrasts = 1:3,
                                     comparison.Name = "Estudi",
                                     outputDir = outputDir,
                                     anotPackage = "org.Hs.eg",
                                     organisme = "hsa",
                                     my.IDs = entrez,
                                     addGeneNames = TRUE,
                                     fileOfLinks = linkfile,
                                     cutoffMethod = c("unadjusted"), 
                                     P.Value.cutoff = rep(0.05, length(1:3)),
                                     pval = 0.05,
                                     thrLogFC = NULL,
                                     minNumGens = 0), 
             priority = 0.56),
  tar_target(KEGGAnalysis.func, do.call(KEGGAnalysis, KEGGAnalysis.pars),
             priority = 0.55),
  
  #doKEGGAnalysis
  tar_target(doKEGGAnalysis.pars, 
             add2parsList(KEGGParsList <- list(),
                          KEGGPar <- list(fitFileName = "fit.Rda",
                                          whichContrasts = 1:3,
                                          comparisonName = "Estudi",
                                          outputDir = outputDir,
                                          anotPackage = "org.Hs.eg", 
                                          organisme = "mmu", 
                                          my.IDs = entrez, 
                                          addGeneNames = TRUE,
                                          fileOfLinks = linkfile,
                                          fitMain = NULL,
                                          cutoffMethod = "unadjusted",
                                          P.Value.cutoff = rep(0.01, length(1:3)),
                                          pvalKEGGterms = 0.05,
                                          minNumGens = 0)),
             priority = 0.54),
  
  tar_target(doKEGGAnalysis.func, for (i in 1:length(doKEGGAnalysis.pars))
  {
    KEGGList = doKEGGAnalysis(doKEGGAnalysis.pars[i])
  },
  priority = 0.53))

  