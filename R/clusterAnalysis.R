plotClust <- function(mydata, ind.cut, i, scal = FALSE, ...)
{
  my.col <- ncol(mydata)
  matData <- matrix(mydata[ind.cut==i, ], ncol = my.col)

  if(scal)
  {
    matData <- (t(scale(t(matData), scale = F)) + mean(matData))
  }

  matplot(1:my.col, t(matData), type = "l", xlab = "", ylab = "log-Expression", las = 2, axes = FALSE, main = paste("Cluster ", i, sep = ""), ...)
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], border = "black")
  axis(2, las = 2, cex.axis = 0.7)
  axis(1, at = 1:my.col, lab = colnames(mydata), las = 2, cex.axis = 0.5)
}
########################################################
write2csv <- function(my.data, fileName, csv = c("csv2", "csv", "txt", "xls"), outputDir)
{
  fileName<- file.path(outputDir, paste(fileName, substr(csv[1], 1, 3) , sep = "."))
  switch (csv[1],
          "csv" = write.csv(my.data, file = fileName, quote = F),
          "csv2" = write.csv2(my.data, file = fileName, quote = F),
          "txt" = write.table(my.data, file = fileName, quote = F))
}
########################################################
#' clusterAnalysis
#'
#' Function to make the cluster analysis.
#' @param expres Expresion oused to perform the clusters.
#' @param genes ????
#' @param samples ????
#' @param sampleNames Name of the samples of the study.
#' @param comparisonName Name of the comparison. To identify the output
#' @param anotPackage Annotation package.
#' @param my.symbols ????
#' @param outputDir Path of the file created.
#' @param fileOfLinks Name of the links file.
#' @param numClusters ????
#' @param rowDistance Row distance in the clusters.
#' @param colDistance Column distance in the clusters.
#' @param RowVals ????
#' @param ColVals ?????
#' @param colorsSet ???
#' @param colsForGroups Colors of each group for the plots.
#' @param escala ???
#' @param densityInfo ???
#' @param cexForColumns Cex for the columns.
#' @param cexForRows Cex for the rows.
#' @param plotProfiles ???
#' @param Title Title of the plot.
#' @param csvType Csv type of the output.
#' @importFrom stats hclust
#' @importFrom stats as.dendrogram
#' @importFrom gplots heatmap.2
#' @importFrom links2File addToLinksFile
#' @exemples
#'  \dontrun{
#' s2clust <-  which(as.logical(apply(design[,as.logical(apply(abs(as.matrix(cont.matrix[,wCont[[i]]])),1,sum))],1,sum)))
#' geneListFName = paste("geneList",  compName[i],     ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"), "LT", pValCutOff[i], "Rda", sep = ".")
#' pal <- colorpanel(n = 32, low = "green", mid = "white", high = "magenta")
#' load("./ResultsDir/geneList.Group1.pvalues.LT.0.01.Rda")
#' clust <- BasicP::clusterAnalysis(expres = exprs.filtered,
#'                                  genes = geneList, 
#'                                  samples = s2clust, 
#'                                  sampleNames =       as.character(targets$ShortName)[s2clust],
#'                                  comparisonName = "Compar 1",
#'                                  anotPackage = "org.Hs.eg",
#'                                  my.symbols = symbolsTable,
#'                                 outputDir = outputDir, 
#'                                  fileOfLinks = linksFile,                              
#'                                  numClusters = 2,
#'                                  rowDistance = NULL,
#'                                  colDistance = NULL,
#'                                  RowVals = TRUE,
#'                                  ColVals = FALSE,
#'                                  escala = "row",
#'                                  colorsSet = pal,
#'                                 densityInfo = "density",
#'                                  colsForGroups =    c("pink","pink","pink","pink","pink","blue","blue","blue","blue","blue"),
#'                                  cexForColumns = 0.8,
#'                                  cexForRows = 0.8,
#'                                  Title = paste("Compar 1 with",
#'                                                ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
#'                                                "<",
#'                                                pValCutOff[i], ifelse(minLogFoldChange[i]==0, "", paste("\n and |logFC|>=", minLogFoldChange[i], sep=""))),
#'                                  csvType = "csv2")
#'}
#' @export


clusterAnalysis <- function(expres,
                            genes,
                            samples,
                            sampleNames,
                            comparisonName,
                            anotPackage,
                            my.symbols=NULL,
                            outputDir,
                            fileOfLinks,
                            numClusters,
                            rowDistance,
                            colDistance,
                            RowVals = TRUE,
                            ColVals = TRUE,
                            colorsSet,
                            colsForGroups,
                            escala,
                            densityInfo = "none",
                            cexForColumns,
                            cexForRows,
                            plotProfiles = FALSE,
                            Title = "",
                            csvType = NULL)
{
  cat ("  - Comparison", Title, "\n")
  categLabel <- "CLUSTER"
  exprs2Cluster <- expres[genes, samples]
  colnames(exprs2Cluster) <- sampleNames

  dendro <- 'both'

  if (RowVals && (!is.numeric(RowVals)))
  {
    clustRow <- hclust(as.dist(1 - cor(t(exprs2Cluster))),  method = "average")
    dendroRow <- as.dendrogram(clustRow)
  }else{
    dendroRow <- RowVals
  }

  if (ColVals && (!is.numeric(ColVals)))
  {
    # Caldria mirar si 'rowDistance' es nula, i si no ho es fer-la servir aqui
    clustCol <- hclust(dist(t(exprs2Cluster)), method = "average")
    dendroCol <-  as.dendrogram(clustCol)
  }else{
    if(is.numeric(ColVals))
    {
      exprs2Cluster <- exprs2Cluster[, ColVals]
      colsForGroups <- colsForGroups[ColVals]
      ColVals <- FALSE
      dendroCol <- FALSE
    }else{
      dendroCol <- ColVals
    }
  }

  if (ColVals && RowVals)
  {
    dendro <- 'both'
  }else{
    if (ColVals && (!RowVals))
    {
      dendro <- 'column'
    }else{
      if ((!ColVals) && RowVals)
      {
        dendro <- 'row'
      }else{
        dendro <- 'none'
      }
    }
  }
  if (toTIFF == TRUE){
    heatMapFName <- paste("HeatMap", comparisonName, "tiff", sep=".")
  }else{
    heatMapFName <- paste("HeatMap", comparisonName, "pdf", sep=".")
  }
  mainTitle <- ifelse (Title == "", comparisonName, Title )

  if(!is.null(outputDir))
  {
    if (toTIFF == TRUE){
      tiff(file = file.path(outputDir, heatMapFName), width = 3200, height = 3200, units = "px", res = 800)
      par(cex.main = 1)
    }else{
      pdf(file.path(outputDir, heatMapFName))
    }
  }

  foundSymbols <- my.symbols[which(names(my.symbols) %in% rownames(exprs2Cluster))]
  newNames <- unlist(lapply(strsplit(paste(names(foundSymbols), foundSymbols ,sep="."), "//"), function(l) l[[1]][1]))
  rownames(exprs2Cluster)[rownames(exprs2Cluster) %in% names(my.symbols)] <- newNames
  #######################

  hm <- heatmap.2(exprs2Cluster,
                  col = colorsSet,
                  ColSideColors = as.character(colsForGroups),
                  scale = escala,
                  Rowv = dendroRow,
                  Colv = dendroCol,
                  dendrogram = dendro,
                  key = TRUE,
                  keysize=1.5,
                  symkey = FALSE,
                  density = "none",
                  trace = "none",                                       # Nomes a heamap.2
                  cexCol = cexForColumns,
                  #labRow = NA,						# Si no volen el nom dels gens
                  cexRow = cexForRows,
                  main = mainTitle)

  if(!is.null(outputDir))
  {
    dev.off()
  }

  addToLinksFile(fileOfLinks, heatMapFName, categ = categLabel, desc = "Heatmap made from genes selected from multiple comparisons")
  ###############################
  if  (RowVals)
  {
    cutN <- cutree(clustRow, numClusters)

    if (ColVals)
    {
      names.ord <- (clustCol$labels[clustCol$order])
    }else{
      names.ord <- 1:ncol(exprs2Cluster)
    }

    if(plotProfiles)
    {
      plotClustFName <- paste("ProfilePlots", comparisonName, "pdf", sep = ".")

      pdf(file.path(outputDir, plotClustFName))
      for(i in 1:numClusters)
      {
        plotClust(exprs2Cluster[, names.ord], cutN, i, scal = T)
      }
      dev.off()

      addToLinksFile(fileOfLinks, plotClustFName, categ = categLabel, desc = "Profiles Plots of clusterized genes selected from multiple comparisons")
    }

    ###########################################################

    raw.gNames <- rev(rownames(exprs2Cluster[hm$rowInd, ]))
    gNames <- unlist(lapply(strsplit(raw.gNames, "[.]"), function(l) l[[1]][1]))

    if(!is.null(my.symbols))
    {
      mySymbols <- my.symbols[gNames]
      my.genes <- data.frame(symbol = mySymbols, ID = gNames, GrupID = cutN[gNames], exprs2Cluster[raw.gNames, names.ord])
    }else{
      if (!is.null(anotPackage)){
        my_SYMBOL_env <- eval(parse(text = paste(anotPackage, "SYMBOL",sep = "")))
        mySymbols <- unlist(mget(gNames, my_SYMBOL_env, ifnotfound=NA))
        my.genes <- data.frame(symbol = mySymbols, ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
      }else{
        my.genes <- data.frame(ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
      }
    }
    genesInClustersFName <- paste("genesInClusters", comparisonName, sep = ".")
    csvType<-ifelse(is.null(csvType), "csv2", csvType)
    write2csv(my.genes, fileName = genesInClustersFName, csv = csvType, outputDir = outputDir)
    addToLinksFile(fileOfLinks, paste(genesInClustersFName,  substr(csvType[1], 1, 3), sep="."),
                   categ = categLabel, desc = "Assignment of each gene to a cluster of made from genes selected from multiple comparisons")
  }

  return(hm)
}
