######################################################
resSelected <- function(res)
{
  sum.res.rows <- apply(abs(res), 1, sum)
  res.selected <- res[sum.res.rows != 0, ]
  if(!is.null(res.selected))
  {
    if (is.matrix(res.selected))
    {
      resSelected <- res.selected
    }else{
      if(ncol(res)>1)
      {
        resSelected <- matrix(res.selected, nrow=1)
      }else{
        resSelected <- matrix(res.selected, ncol=1)
        rownames(resSelected)<-names(res.selected)
      }

    }
  }else{
    resSelected <- NULL
  }

  return(resSelected)
}
######################################################

write2csv <- function(my.data, fileName, csv = c("csv2", "csv", "txt", "xls"), outputDir)
{
  fileName<- file.path(outputDir, paste(fileName, substr(csv[1], 1, 3) , sep = "."))
  switch (csv[1],
          "csv" = write.csv(my.data, file = fileName, quote = F),
          "csv2" = write.csv2(my.data, file = fileName, quote = F),
          "txt" = write.table(my.data, file = fileName, quote = F))
}
######################################################
plotVennDiagram <- function(res.selected,
                            colsVenn,
                            whichContrasts,
                            vennColors = c("red", "yellow", "green", "blue", "pink"),
                            comparisonName,
                            titleText,
                            outputDir,
                            linksFile,
                            categLabel,
                            my.cex)
{
  if (is.null(colsVenn)) colsVenn <-1:length(whichContrasts)
  if (length(colsVenn)> 5) stop("venn.diagram cannot plot more than 5 elements at the same time")

  vennFName  <- paste("venn", comparisonName, sep = ".")
  ncomp <- as.character(length(colsVenn))

  switch(ncomp,
         "2" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0)),
                                          category.names = colnames(res.selected)[1:2],
                                          fill = c(vennColors[1], vennColors[2]),
                                          alpha = 0.30,
                                          main = paste(comparisonName, titleText, sep=" "),
                                          resolution = 600,
                                          cat.cex = my.cex,
                                          filename = NULL)},
         "3" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0),
                                               C = which(res.selected[,3]!=0)),
                                          category.names = colnames(res.selected)[1:3],
                                          fill = c(vennColors[1], vennColors[2], vennColors[3]),
                                          alpha = 0.30,
                                          resolution = 600,
                                          cat.cex = my.cex,
                                          main = paste(comparisonName, titleText, sep=" "),
                                          filename = NULL)},
         "4" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0),
                                               C = which(res.selected[,3]!=0), D = which(res.selected[,4]!=0)),
                                          category.names = colnames(res.selected)[1:4],
                                          fill = c(vennColors[1], vennColors[2], vennColors[3], vennColors[4]),
                                          alpha = 0.30,
                                          resolution = 600,
                                          cat.cex = my.cex,
                                          main = paste(comparisonName, titleText, sep=" "),
                                          filename = NULL)},
         "5" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0),
                                               C = which(res.selected[,3]!=0), D = which(res.selected[,4]!=0), W = which(res.selected[,5]!=0)),
                                          category.names = colnames(res.selected)[1:5],
                                          fill = c(vennColors[1], vennColors[2], vennColors[3], vennColors[4], vennColors[5]),
                                          alpha = 0.30,
                                          cat.cex = my.cex,
                                          resolution = 600,
                                          cat.pos = c(0,0,-140,-140,0),
                                          cat.default.pos = "outer",
                                          main = paste(comparisonName, titleText, sep=" "),
                                          filename = NULL)})

  pdf(file.path(outputDir,paste(vennFName, "pdf", sep = ".")))
  grid.draw(venn.plot)
  dev.off()
  addToLinksFile(linksFile, paste(vennFName, "pdf", sep = "."), categ = categLabel, desc = "Venn Diagram for genes selected from multiple comparisons")

  if(toTIFF){
    tiff(file.path(outputDir,paste(vennFName,"tiff", sep = ".")), width = 3200, height = 3200, units = "px", res = 800)
    grid.draw(venn.plot)
    dev.off()

    addToLinksFile(linksFile, paste(vennFName, "tiff", sep = "."), categ = categLabel, desc = "Venn Diagram for genes selected from multiple comparisons")
  }



}
######################################################
#' multipleComp
#'
#' Function to make multiple comparations.
#'
#' @param fitMain
#' @param whichContrasts Vector that determie the contrast to do.
#' @param comparisonName Name of the comparisons.
#' @param titleText
#' @param outputDir Path of the file created.
#' @param anotPackage Annotation package.
#' @param my.symbols
#' @param linksFile Name of the LinksFile.
#' @param multCompMethod Method to make the comparation.
#' @param adjustMethod Method to apply in the adjust.
#' @param selectionType
#' @param P.Value.cutoff Cutoff value for p.value.
#' @param plotVenn
#' @param colsVenn
#' @param vennColors
#' @param cexVenn
#' @param geneListFName
#' @param csvType
#' @param minLFC
#' @importFrom limma decideTests
#' @importFrom grDevices pdf
#' @importFrom grid grid.draw
#' @importFrom grDevices tiff
#' @importFrom VennDiagram venn.diagram
#' @importFrom links2File addToLinksFile
#' @examples 
#' \dontrun{
#' compName <- c("Group1")
#' wCont <- 1:3
#' pValCutOff <- c(0.01)
#' adjMethod <- c("none")
#' minLogFoldChange <- c(1)
#' load("./ResultsDir/Symbols.Rda")
#' titleText = paste("for", ifelse(adjMethod=="none","p-values","adj. p-values"), "<", pValCutOff, "and |logFC| >", minLogFoldChange, sep = " ")
#' geneListFName = paste("geneList",compName,ifelse(adjMethod=="none","pvalues","adj-pvalues"),"LT",pValCutOff,"Rda",sep = ".")
#' 
#' geneList <-  BasicP::multipleComp(fitMain = fitMain, whichContrasts = wCont, comparisonName = compName, titleText = titleText,
#' outputDir = outputDir,  anotPackage = "org.Hs.eg", my.symbols = symbolsTable, linksFile = linksFile, multCompMethod = "separate",
#' adjustMethod = adjMethod, selectionType = "any", P.Value.cutoff = pValCutOff, plotVenn = TRUE, colsVenn = NULL, 
#' vennColors = c("red","yellow","green","blue","pink") , cexVenn = 1, geneListFName=geneListFName, csvType=csvType, minLFC=minLogFoldChange)
#' }
#' @export

multipleComp <- function(fitMain,
                         whichContrasts,
                         comparisonName,
                         titleText,
                         outputDir,
                         anotPackage,
                         my.symbols=NULL,
                         linksFile,
                         multCompMethod = "separate",
                         adjustMethod = "none",
                         selectionType = c("any", "all", "anyDown", "anyUp", "allDown", "allUp"),
                         P.Value.cutoff = 0.05,
                         plotVenn,
                         colsVenn = NULL,
                         vennColors,
                         cexVenn = 1,
                         geneListFName=paste("geneList", comparisonName,"Rda", sep="."),
                         csvType = NULL,
                         minLFC = 0)
{
  geneList <- NULL
  categLabel <- "MULTCOMP"
  if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[,whichContrasts])) && (is.matrix(fitMain$contrasts[,whichContrasts])))
  {
    res <- decideTests(fitMain[, whichContrasts], method = multCompMethod, adjust.method = adjustMethod, p.value = P.Value.cutoff, lfc = minLFC)
    res.selected <- resSelected(res)

    if (!is.null(res.selected))
    {

      if (!is.null(my.symbols))
      {
        gNames <- rownames(res.selected)

        my.indNA<-which(is.na(my.symbols[gNames]))
        my.symbols[gNames[my.indNA]]<-my.symbols[gNames[my.indNA]]

        symbols.selected <- my.symbols[gNames]
        res.selected2 <- cbind(SYMBOLS=symbols.selected, res.selected)
      }else{
        if (!is.null(anotPackage))
        {
          stopifnot(require(old2db(anotPackage), character.only = T))
          myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
          symbols.selected <- unlist(mget(rownames(res.selected), env = myenvirSYMBOL, ifnotfound = NA))
          res.selected2 <- cbind(SYMBOLS=symbols.selected, res.selected)
        }else{
          res.selected2 <- res.selected
        }
      }
      selectedFName <- paste("multComp", comparisonName, sep = ".")

      csvType<-ifelse(is.null(csvType), "csv2", csvType)
      write2csv(res.selected2, fileName = selectedFName, csv = csvType, outputDir = outputDir)
      addToLinksFile(linksFile, paste(selectedFName, substr(csvType[1], 1, 3), sep="."), categ = categLabel, desc="Genes selected from multiple comparisons")

      if (plotVenn)
      {
        plotVennDiagram(res.selected, colsVenn, vennColors, whichContrasts = whichContrasts, comparisonName = comparisonName, titleText, outputDir = outputDir, linksFile, categLabel, my.cex = cexVenn)
      }
    }

    geneList <- rownames(res.selected)
    save(geneList, file = file.path(outputDir, geneListFName))
  }else{
    if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[,whichContrasts])) && (!is.matrix(fitMain$contrasts[,whichContrasts])))
    {
      res <- decideTests(fitMain[, whichContrasts], method = multCompMethod, adjust.method = adjustMethod, p.value = P.Value.cutoff, lfc = minLFC)
      res.selected <- resSelected(res)
      geneList <- rownames(res.selected)
      save(geneList, file = file.path(outputDir, geneListFName))
    }
  }

  return(geneList)
}

