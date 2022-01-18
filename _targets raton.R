#Cargar targets
library(targets)

#Establecer directorio principal
setwd("C:/Users/Juan/Desktop/NOVO/mouse")

#Parámetros generales
runMulticore = 3
toTIFF = FALSE

specie_used = "Mus_musculus" #especie para anotacion de genes con repositorios KEGG y ENSEMBLE
# actualmente para KEGG solo estan disponibles los siguientes:
#Homo_sapiens, Mus_musculus, Rattus_norvegicus, Bos_taurus, Danio_rerioï, Sus_scrofa

orgPackage_used = "org.Mm.eg" #orgPackage que se utilizara
annotation_data_affymetrix = "mouse4302.db" #libreria del organismo correspondiente
#da problemas entre ponerle o no .db al final

#Librerías necesarias
libraries = c(
  "BiocManager",
  "parallel",
  "GOstats",
  "SortableHTMLTables",
  "VennDiagram",
  "grDevices",
  "grid",
  "stats",
  "gplots",
  "DBI",
  "affy",
  "annotate",
  "Biobase",
  "BiocGenerics",
  "limma",
  "oligo",
  "pd.clariom.s.mouse.ht", #importante para el raton
  "org.Mm.eg.db" #importanto para el raton
  )

##Intento de automatizar la comprobación e insalación de los paquetes necesarios
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("mouse4302.db")
#BiocManager::install("pd.clariom.s.mouse.ht") para normalization
#maybe? no ha dado warning ni nada pero no se 100%

#for (i in normal_libraries){
#  p_load(i)
#}
#p_load(normal_libraries, install = TRUE)

#for (i in libraries){
#  if(!require(i)){
#    BiocManager::install(i)}
#}

#Carga de todos los scripts en la carpeta R
sapply(paste("R/", list.files(paste("C:/Users/dirección carpeta R", sep="")), sep=""), source)

#librerías a utilizar
tar_option_set(
  packages = libraries
  )

#Listado de targets
list(
  #ReadOrLoad.RawData
  tar_target(readOrLoad.rawData.pars, 
             list(readCELS = TRUE,
                  phenoDat = "./celfiles/targets.txt",
                  fileNames = paste("./celfiles/", 
                                    rownames(read.table("./celfiles/targets.txt", 
                                                        head=TRUE, sep="\t",  
                                                        row.names = 1)),sep=""),
                  dataFName = "rawData.Rda",
                  outputDir = ".",
                  exonSt = TRUE), #cambiado para el raton
             priority = 1),

  tar_target(readOrLoad.rawData.func, 
             do.call(readOrLoad.RawData, readOrLoad.rawData.pars),
             priority = 0.99),
  
  ####outputDir####
  tar_target(outputDir, "./ResultsDir", priority = 0.98),
  
  #CreateOrLoadAnnotations
  tar_target(createOrLoadAnnotations.pars, list(loadAnnotations = TRUE, 
                                                chipPackAvailable = NULL, #este como null
                                                platformDesignPackAvailable = TRUE, #y este pasa a ser true
                                                chipPackage = NULL, #este pasaa ser null
                                                platformDesignPackage = "pd.clariom.s.mouse.ht", #este es el paquete que se usa
                                                outputDir = outputDir,
                                                annotationsFileName = "Annotations",
                                                entrezTableFileName = "Entrezs.Rda",
                                                symbolsTableFileName = "Symbols.Rda",
                                                controlsTableFileName = "controls.Rda"),
             priority = 0.97),
  
  tar_target(createOrLoadAnnotations.func, 
             do.call(createOrLoadAnnotations,createOrLoadAnnotations.pars),
             priority = 0.96),
  
  #Normalization
  tar_target(normalization.pars, 
             list(my.data = get(load("./rawData.Rda")),
                  method = "RMA",

                  inputDir = "./celfiles",
                  loadFile = FALSE,
                  normalizedFName = "normalizedData.Rda",
                  outputDir = outputDir,
                  exonSt = FALSE),
             priority = 0.95),
  
  tar_target(normalization.func, 
             do.call(normalization, normalization.pars),
             priority = 0.94),
  
  ####Normalized####
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
  
  ####Entrez, Controls y targets####
  tar_target(entrez, get(load("./ResultsDir/Entrezs.Rda")), priority = 0.90),
  
  tar_target(controls, get(load("./ResultsDir/controls.Rda")), priority = 0.89),
  
  tar_target(targets, read.AnnotatedDataFrame("./celfiles/targets.txt", 
                                              header = TRUE,row.names = 1), 
             priority = 0.88),
  
  #FilterData
  tar_target(filterData.pars, #notas 
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

  ####Symbols y linkfile####
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

  ####designmatrix####
  tar_target(designmatrix, function(datos, columnas, filas){
    design = model.matrix(~ 0 + datos)
    colnames(design) = columnas
    rownames(design) = filas
    return(design)
  },
  priority = 0.80),
  
  ####expres.filtered y targets.designmatrix####
  tar_target(expres.filtered, get(load("./ResultsDir/expres.filtered.Rda")), priority = 0.79),
  tar_target(targets.designmatrix, read.table("./celfiles/targets.txt",
                                        head=TRUE,sep="\t"), priority = 0.78),
  
  tar_target(gr, unique(targets.designmatrix$Group), priority = 0.775),
  
  #Esto facilitaría la etapa lmAnalysis pero solo para tres grupos distintos, si son más o menos daría problemas
  tar_target(myargs, list(paste(gr[1],"-",gr[2], sep=""), 
                          paste(gr[1],"-",gr[3], sep=""),
                          paste(gr[2],"-",gr[3], sep=""),
                          levels=designmatrix(
                            datos = targets.designmatrix$Group,
                            columnas = unique(targets.designmatrix$Group),
                            filas = targets.designmatrix$ShortName)), priority = 0.774),
  
  #lmAnalysis
  tar_target(lmAnalysis.pars, 
             list(exprs.filtered = expres.filtered, 
                  
                  design = myargs$levels, #con $ lo lee como matriz mientras que con [4] lo lee como lista y da error
                  
                  cont.matrix = do.call(makeContrasts, myargs),
                  
                  contrasts2test = 1:ncol(makeContrasts(
                    paste(tar_read(gr)[1],"-",tar_read(gr)[2], sep=""), #añadido el tar_read
                    paste(tar_read(gr)[1],"-",tar_read(gr)[3], sep=""), 
                    paste(tar_read(gr)[2],"-",tar_read(gr)[3], sep=""), 
                    levels=unique(targets.designmatrix$Group))),
                  
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
                                 ENTREZIDs = entrez, #notas
                                 SYMBOLIDs = symbols, #notas 
                                 fileOfLinks = linkfile,
                                 fitFileName = "fit.Rda",
                                 csvType= "csv",
                                 rows2HTML = NULL,
                                 anotFilename = "Annotations")),
             priority = 0.75),

  tar_target(doLmAnalysis.func, #en script lmAnalysis.R está también la funcion dolmanalysis, en esta se usa eval(parse(p$symbolids)) que da error, quito el eval(parse) y funciona
             for(ix in 1:length(doLmAnalysis.pars)){
               fit.Main <- doLmAnalysis(doLmAnalysis.pars[ix])
               },
             priority = 0.74),
  
  ####fitmain####
  tar_target(fitmain, get(load("./ResultsDir/fit.Rda")), priority = 0.73),

  #GeneAnnotation
  tar_target(GeneAnnotation.pars, list(egIDs = entrez[unique(rownames(fitmain$p.value))],
                                       anotPackage = orgPackage_used,
                                       toHTML = TRUE,
                                       outputDir = outputDir,
                                       filename = "Annotations",
                                       myTitle = "Annotations for all genes analyzed",
                                       specie = specie_used,
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
                                         anotPackage = orgPackage_used,
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
                                         organisme = "hsa")),#notas
             priority = 0.70),

  tar_target(doGeneAnnotations.func, 
             do.call(doGeneAnnotation,doGeneAnnotation.pars),
             priority = 0.69),

  ####multipleComp.var####
  tar_target(multipleComp.var, list(adjMethod = c("none"),
                                    pValcutoff = 0.01,
                                    compName = c("Group1"),
                                    minLogFoldChange = 1), 
             priority = 0.68),
  
  #multipleComp
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
                  anotPackage = orgPackage_used, #variables establecida la principio
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
                  geneListFName =  paste("geneList",
                                         multipleComp.var$compName,
                                         ifelse(multipleComp.var$adjMethod==
                                                  "none","pvalues","adj-pvalues"),
                                         "LT",
                                         multipleComp.var$pValcutoff,
                                         "Rda", sep = "."),
                  csvType = "csv",
                  minLFC = 1),
             priority = 0.67),

  tar_target(multipleComp.func, 
             do.call(multipleComp, multipleComp.pars),
             priority = 0.66),

  ####doMultCompanalysis.var####
  tar_target(doMultCompAnalysis.var, list(wCont = (1:ncol(lmAnalysis.pars$cont.matrix)),
                                          compName = c("Group1", "Group2"),
                                          adjMethod = c("none", "none"),
                                          pValCutOff = c(0.01, 0.01),
                                          minLogFoldChange = c(1, 1)),
             priority = 0.65),

  #doMultCompAnalysis 
  tar_target(doMultCompAnalysis.pars, 
             do.call(sets_doMultCompAnalysis, doMultCompAnalysis.var),
             priority = 0.64),
  
  tar_target(doMultCompAnalysis.func, 
             for (ix in 1:length(doMultCompAnalysis.pars))
               {
               geneList.MCi = doMultCompAnalysis(doMultCompAnalysis.pars[ix])
               }, 
             priority = 0.63),
  
  ####s2clust####
  tar_target(s2clust, 
             which(as.logical(apply(lmAnalysis.pars$design[
               ,
               as.logical(apply(abs(as.matrix(lmAnalysis.pars$cont.matrix[
                 ,
                 (1:ncol(lmAnalysis.pars$cont.matrix))[[1]]])),1,sum))],1,sum))),
             priority = 0.62),
  
  tar_target(clusterAnalysis.pars, 
             list(expres = expres.filtered,
                  genes = get(load("./ResultsDir/geneList.Group1.pvalues.LT.0.01.Rda")),
                  samples = s2clust,
                  sampleNames = as.character(targets$ShortName)[s2clust],
                  comparisonName = "Compar 1",
                  anotPackage = orgPackage_used,
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
                  #s2clust es un int de 8 valores, por lo que la longitud del vector colsforgroups debe coincidir con el de s2clust, 
                  #en vez de indicarlos directamente con
                  #strings, debería crear una función que extraiga la longitud de s2clust y decida los colores en base a ese número.
                  colsForGroups = c("pink","pink","pink","pink",
                                    "blue","blue","blue","blue"),
                  cexForColumns = 0.8,
                  cexForRows = 0.8,
                  Title = paste("Compar 1 with",
                                ifelse(c("none") == "none", "pvalues", "adj-pvalues"),
                                "<",
                                c(0.01), ifelse(1==0, "", 
                                                paste("\n and |logFC|>=", 
                                                      1, sep=""))),
                  csvType = "csv2"),
             priority = 0.61),

  tar_target(clusterAnalysis.func, do.call(clusterAnalysis, clusterAnalysis.pars),
             priority = 0.60),


  ####doClusterAnalysis.var####
  tar_target(doClusterAnalysis.var, list(compName = c("Group1", "Group2"),
                                         adjMethod = c("none", "none"),
                                         pValCutOff = c(0.01, 0.01),
                                         minLogFoldChange = c(1, 1)),
             priority = 0.59),
  
  #doClusterAnalysis
  tar_target(doClusterAnalysis.pars, do.call(sets_doclusterAnalysis,doClusterAnalysis.var),
             priority = 0.58),

  tar_target(doclusterAnalysis.func,
             for (ix in 1:length(doClusterAnalysis.pars))
             {
               hm.Estudi <- doClusterAnalysis(doClusterAnalysis.pars[ix])
             },
             priority = 0.57),

  #GoAnalysis  
  tar_target(GoAnalysis.pars, list(fitMain = fitmain,
                                   whichContrasts = (1:ncol(lmAnalysis.pars$cont.matrix)),
                                   comparison.Name = "Estudi",
                                   outputDir = outputDir,
                                   anotPackage = orgPackage_used,
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
             priority = 0.56),

  tar_target(GOAnalysis.func,  do.call(GOAnalysis, GoAnalysis.pars),
             priority = 0.55),

  #doGoAnalysis
  tar_target(doGOAnalysis.pars,
             add2parsList(GOParsList <- list(),
                          GOPar <- list(fitMain = fitmain,
                                        whichContrasts = 1:ncol(lmAnalysis.pars$cont.matrix),
                                        comparison.Name = "Estudi",
                                        outputDir = outputDir,
                                        anotPackage = orgPackage_used,
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
             priority = 0.54),
  
  tar_target(doGOAnalysis.func, for (i in 1:length(doGOAnalysis.pars))
    {
    GOList = doGOAnalysis(doGOAnalysis.pars[i])
  }, 
             priority = 0.53),

  #KEGGAnalysis
  tar_target(KEGGAnalysis.pars, list(fitMain = fitmain,
                                     whichContrasts = 1:3,
                                     comparison.Name = "Estudi",
                                     outputDir = outputDir,
                                     anotPackage = orgPackage_used,
                                     organisme = "mmu",
                                     my.IDs = entrez,
                                     addGeneNames = TRUE,
                                     fileOfLinks = linkfile,
                                     cutoffMethod = c("unadjusted"), 
                                     P.Value.cutoff = rep(0.05, length(1:3)),
                                     pval = 0.05,
                                     thrLogFC = NULL,
                                     minNumGens = 0), 
             priority = 0.52),
  
  tar_target(KEGGAnalysis.func, do.call(KEGGAnalysis, KEGGAnalysis.pars),
             priority = 0.51),


  #doKEGGAnalysis 
  tar_target(doKEGGAnalysis.pars, 
             add2parsList(KEGGParsList <- list(),
                          KEGGPar <- list(fitFileName = "fit.Rda",
                                          whichContrasts = 1:3,
                                          comparisonName = "Estudi",
                                          outputDir = outputDir,
                                          anotPackage = orgPackage_used, 
                                          organisme = "mmu", 
                                          my.IDs = entrez, 
                                          addGeneNames = TRUE,
                                          fileOfLinks = linkfile,
                                          fitMain = NULL,
                                          cutoffMethod = "unadjusted",
                                          P.Value.cutoff = rep(0.01, length(1:3)),
                                          pvalKEGGterms = 0.05,
                                          minNumGens = 0)),
             priority = 0.50),
  
  tar_target(doKEGGAnalysis.func, for (i in 1:length(doKEGGAnalysis.pars))
    {
    KEGGList = doKEGGAnalysis(doKEGGAnalysis.pars[i])
  },
             priority = 0.49))
