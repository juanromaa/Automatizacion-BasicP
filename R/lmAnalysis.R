##################################################################
is.oneSampleContrast <-function (row.cont.matrix){
  (length(row.cont.matrix[row.cont.matrix!=0])==1) && (sum(row.cont.matrix)==1)
}
###################################################################
is.twoSampleContrast <-function (row.cont.matrix){
  (length(row.cont.matrix[row.cont.matrix!=0])==2) && (sum(row.cont.matrix)==0)
}
###################################################################
genesSelectable <- function (topTab, adj0, adj1, adj2, P1, P2)
{
  upBelowB <- sum(topTab$B > 0  & topTab$t > 0)
  downBelowB <- sum(topTab$B >0 & topTab$t < 0)

  upBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t > 0)
  downBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t < 0)

  upBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t > 0)
  downBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t < 0)

  upBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t > 0)
  downBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t < 0)

  upBelowP1 <- sum(topTab$P.Value < P1 & topTab$t > 0)
  downBelowP1 <- sum (topTab$P.Value < P1 & topTab$t < 0)

  upBelowP2 <- sum(topTab$P.Value < P2 & topTab$t > 0)
  downBelowP2 <- sum(topTab$P.Value < P2 & topTab$t < 0)

  return(c(upBelowB = upBelowB,downBelowB = downBelowB,
           upBelowAdj0 = upBelowAdj0, downBelowAdj0 = downBelowAdj0,
           upBelowAdj1 = upBelowAdj1, downBelowAdj1 = downBelowAdj1,
           upBelowAdj2 = upBelowAdj2, downBelowAdj2 = downBelowAdj2,
           upBelowP1 = upBelowP1, downBelowP1 = downBelowP1,
           upBelowP2 = upBelowP2, downBelowP2 = downBelowP2))
}

###################################################################
genesSelected <- function (topTab = NULL, adj0=0.01, adj1 = 0.05, adj2 = 0.25, P1 = 0.01, P2 = 0.05)
{
  if (is.null(topTab))
  {
    seleccio <- data.frame(Col1 = integer(length = 12),
                           row.names = c("upReg-B>0", "downReg-B>0",
                                         paste("upReg-Adjusted-p-val", adj0, sep = " < "),
                                         paste("downReg-Adjusted-p-val", adj0, sep = " < "),
                                         paste("upReg-Adjusted-p-val", adj1, sep = " < "),
                                         paste("downReg-Adjusted-p-val", adj1, sep = " < "),
                                         paste("upReg-Adjusted-p-val", adj2, sep = " < "),
                                         paste("downReg-Adjusted-p-val", adj2, sep = " < "),
                                         paste("upReg-P value", P1, sep = " < "),
                                         paste("downReg-P value", P1, sep = " < "),
                                         paste("upReg-P value", P2, sep = " < "),
                                         paste("downReg-P value", P2, sep = " < ")))[, -1]
  }else{
    genesSelectable(topTab, adj0, adj1, adj2, P1, P2)
  }

  return(seleccio)
}

###################################################################
knowledge2gNames <- function(my.genes, knowledgeFile)
{
  my.genes <- as.character(my.genes)

  my.annotation <- character()

  if (!is.null(knowledgeFile))
  {
    for (i in 1:length(my.genes))
    {
      my.annotation[i] <- paste("<A HREF=\"", knowledgeFile,"#", my.genes[i] ,"\" TARGET=\"_blank\">", my.genes[i],"</A>", sep = "")
    }
  }
  return(my.annotation)
}
###################################################################
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
###################################################################
midSinonims <- function(my.IDs)
{
  my.indexes <- grep(" /// ", my.IDs)
  my.sinonimIDs <- my.IDs[my.indexes]

  my.IDs[my.indexes] <- extractSinonims(my.sinonimIDs)
  return(my.IDs)
}

###################################################################

annotateTopTable2 <- function (topTab,
                               fName,
                               Title = "Genes selected",
                               anotPackage,
                               EntrezIDs = NULL,
                               SymbolIDs = NULL,
                               #comparison = "",
                               anotFilename = "annotations"
)
{
  require("annotate",  character.only = TRUE)
  if (!is.null(topTab$ID)){
    gNames <- as.character(as.integer(topTab$ID))
  }else{
    gNames <- as.character(rownames(topTab))
  }

  linkedGeneNanes <- knowledge2gNames(gNames, knowledgeFile =  paste(anotFilename, "html", sep="."))

  if (is.null(EntrezIDs))
  {
    myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
    EntrezIDs <- unlist(mget(gNames, env = myenvirENTREZID, ifnotfound = NA))
  }else{
    EntrezIDs <- midSinonims(EntrezIDs[gNames])
  }

  if (is.null(SymbolIDs))
  {
    myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
    SymbolIDs <- unlist(mget(gNames, env = myenvirSYMBOL, ifnotfound = NA))
  }else{
    SymbolIDs <- midSinonims(SymbolIDs[gNames])
  }

  aux.SymbolIDs <- as.character(SymbolIDs)
  names(aux.SymbolIDs) <- names(SymbolIDs)
  if (runMulticore ==1 || runMulticore ==3) {
    lS <- mclapply(aux.SymbolIDs, function(l){unlist(strsplit(l, "//"))})
  } else {
    lS <- lapply(aux.SymbolIDs, function(l){unlist(strsplit(l, "//"))})
  }

  linkedList <- list(en = EntrezIDs)
  if (!is.null(topTab$ID)){ topTab <- topTab[, which(names(topTab)=="ID")]}

  # Si existeixen sinonims "otherNames" ha de ser una llista i no un data.frame com estava
  otherNames = list(affyIDs = linkedGeneNanes, GeneSymbols = lS, topTab)
  htmlpage(linkedList, filename = fName, title = Title, othernames = otherNames,
           table.head = c("EntrezID", "affyIDs", "GeneSymbols", names(topTab)),
           table.center = TRUE, repository = list("en"), digits=4)
}

###################################################################
write2csv <- function(my.data, fileName, csv = c("csv2", "csv", "txt", "xls"), outputDir)
{
  fileName<- file.path(outputDir, paste(fileName, substr(csv[1], 1, 3) , sep = "."))
  switch (csv[1],
          "csv" = write.csv(my.data, file = fileName, quote = F),
          "csv2" = write.csv2(my.data, file = fileName, quote = F),
          "txt" = write.table(my.data, file = fileName, quote = F))
  # }
}
###################################################################

escriuTop_i_Express <- function(expres,
                                topTab,
                                grup1 = NULL,
                                grup2 = NULL,
                                nom1,
                                nom2,
                                fName = NULL,
                                fit = NULL,
                                fitCoef = NULL,
                                anotPackage,
                                my.symbols = NULL,
                                my.entrezs = NULL,
                                csvType,
                                outputDir)
{
  if (!is.null(grup2)) {
    expresA <- expres[, c(grup1, grup2)]

    means2groups <- function(x)
    {
      mean1 <- mean(x[1:length(grup1)])
      mean2 <- mean(x[(length(grup1) + 1):(length(grup1) + length(grup2))])
      foldC <- mean1 - mean2

      return(c(foldC, mean1, mean2))
    }

    meansA <- t(apply(expresA, 1, means2groups))
    colnames(meansA) <- c("logFC validation", nom1, nom2)

  }else if (!is.null(grup1)) {   # en aquest cas sols hi ha grup1
    expresA <- expres[, grup1]

    means1groups <- function(x) {
      foldC <- mean(x[1:length(grup1)])
      return(c(foldC, foldC))
    }

    meansA <- t(apply(expresA, 1, means1groups))
    colnames(meansA) <- c("logFC validation", "logFC validation")

  }else {
    meansA <- NULL
    expresA <- expres
  }

  if (!is.null(topTab$ID)){
    topA.order <- as.character(as.integer(topTab$ID))
  }else{
    topA.order <- rownames(topTab)
  }

  if (!is.null(meansA)) { meansA.ord <- meansA[topA.order, ] }
  expresA.ord <- expresA[topA.order, ]
  fitA.ord <- fit[topA.order, ]

  if (! is.null(anotPackage))
  {
    SymbolsA <- getSYMBOL(topA.order, old2db (anotPackage))
  }else{
    if (!is.null(my.symbols))
    {
      gNames <- topA.order
      SymbolsA <- my.symbols[gNames]
      EntrezsA <- my.entrezs[gNames]
    }
  }

  # calculem els parametres del limma
  if(!is.null(fit) && !is.null(fitCoef))
  {
    s.post <- sqrt(fitA.ord$s2.post)
    s.unscaled <- fitA.ord$stdev.unscaled[, fitCoef]
    df.res <- fitA.ord$df.residual

    s.res <- fitA.ord$sigma
    s.prior <- rep(sqrt(fitA.ord$s2.prior), length(s.res))
    df.prior <- rep(fitA.ord$df.prior, length(s.res))

    fit.coef <- fitA.ord$coef[, fitCoef]
    t.ord <-  fit.coef / (s.unscaled * s.res)
    p.ord <- 2 * pt (abs(t.ord), df = df.res, lower.tail = FALSE)

    limmaPars <- data.frame(s.post, s.unscaled, t.ord, p.ord, df.res, s.res, s.prior, df.prior)
  }

  if (!is.null(fit) && !is.null(fitCoef)) {  # mirem si tenim limmaPars
    if (! is.null(my.symbols)) {  # mirem si tenim symbols
      if (! is.null(my.entrezs)) {  # mirem si tenim entrezs
        if (! is.null(meansA)) {  # mirem si hi ha mitjanes (o sigui que els contrasts son simples, e.g: A - B )
          combinedA <- cbind(SymbolsA, EntrezsA, topTab, meansA.ord, expresA.ord, limmaPars)
        }else {  # si no hi ha mitjanes (o sigui que els contrasts son simples, e.g: (A-C) - (B-C) )
          combinedA <- cbind(SymbolsA, EntrezsA, topTab, expresA.ord, limmaPars)
        }
      }else { # si no tenim entrezs
        if (! is.null(meansA)) {  # si hi ha mitjanes (o sigui que els contrasts son simples, e.g: A - B )
          combinedA <- cbind(SymbolsA, topTab, meansA.ord, expresA.ord, limmaPars)
        }else {  # si no hi ha mitjanes (o sigui que els contrasts son simples, e.g: (A-C) - (B-C) )
          combinedA <- cbind(SymbolsA, topTab, expresA.ord, limmaPars)
        }
      }
    }else {  # si no tenim symbols
      if (! is.null(my.entrezs)) {
        if (!is.null(meansA)) {
          combinedA <- cbind(EntrezsA, topTab, meansA.ord, expresA.ord, limmaPars)
        }else {
          combinedA <- cbind(EntrezsA, topTab, expresA.ord, limmaPars)
        }
      }else {
        if (! is.null(meansA)) {  # si hi ha mitjanes (o sigui que els contrasts son simples, e.g: A - B )
          combinedA <- cbind(topTab, meansA.ord, expresA.ord, limmaPars)
        }else {  # si no hi ha mitjanes (o sigui que els contrasts son simples, e.g: (A-C) - (B-C) )
          combinedA <- cbind(topTab, expresA.ord, limmaPars)
        }
      }
    }
  }else {  # si no tenim limmaPars
    if (!is.null(my.symbols)) {
      if (! is.null(my.entrezs)) {
        if (! is.null(meansA)) {
          combinedA <- cbind(SymbolsA, EntrezsA, topTab, meansA.ord, expresA.ord)
        }else {
          combinedA <- cbind(SymbolsA, EntrezsA, topTab, expresA.ord)
        }
      }else {
        if (! is.null(meansA)) {
          combinedA <- cbind(SymbolsA, topTab, meansA.ord, expresA.ord)
        }else {
          combinedA <- cbind(SymbolsA, topTab, expresA.ord)
        }
      }
    }else {
      if (! is.null(my.entrezs)) {
        if (!is.null(meansA)) {
          combinedA <- cbind(EntrezsA, topTab, meansA.ord, expresA.ord)
        }else {
          combinedA <- cbind(EntrezsA, topTab, expresA.ord)
        }
      }else {
        if (!is.null(meansA)) {
          combinedA <- cbind(topTab, meansA.ord, expresA.ord)
        }else {
          combinedA <- cbind(topTab, expresA.ord)
        }
      }
    }
  }

  if (!is.null(fName)) write2csv(combinedA, fileName = fName, csv = csvType, outputDir = outputDir)

  return(combinedA)
}
###################################################################
### lmAnalysis
###################################################################

#'lmAnalysis
#'
#'Function to make the linear model analysis.
#'
#' @param exprs.filtered Dataset of filtered and normalized data.
#' @param design Design matrix.
#' @param cont.matrix Contrast matrix.
#' @param contrasts2test Numeric vector that contains the number of the columns for the contrasts to analyze.
#' @param anotPackage Annotation package.
#' @param Expressions_And_Top If is TRUE it will do "Expresions" and the toptable.
#' @param showParams If FALSe it won't show the lmFit parameters in the ExpressionsAndToptable.
#' @param use.dupCorr  ???
#' @param block ???
#' @param nDups ???
#' @param comparison Name of the comparations in all contrasts.
#' @param outputDir Path of the file created.
#' @param ENTREZIDs Name of the file for Entrez genes.
#' @param SYMBOLIDs Name of the file for gene Symbols .
#' @param linksFile Name of the file linksFile.
#' @param fitFileName ???
#' @param csvType type of csv to store the data.
#' @param rows2HTML ???
#' @param anotFileName ???
#' @importFrom limma duplicateCorrelation
#' @importFrom limma lmFit
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom limma volcanoplot
#' @importFrom annotate htmlpage
#' @importFrom annotate getSYMBOL
#' @importFrom links2File addToLinksFile
#' @exmaples
#' \dontrun{
#' load("./ResultsDir/exprs.filtered.Rda")
#' contrasts2test <- 1:ncol(cont.matrix)
#' anotPackage = NULL
#' comparison =  "Estudi"
#' outputDir = "./ResultsDir"
#' ENTREZIDs = "entrezTable"
#' SYMBOLIDs = "symbolsTable"
#' linksFile = "Links.txt"
#' fitFileName = "fit.Rda"
#' csvType= "csv"
#' rows2HTML= NULL
#' anotFileName <- "Annotations"
#' runMulticore = 0 
#' toTIFF= FALSE
#' fitMain <- BasicP::lmAnalysis(exprs.filtered = exprs.filtered, design = design, 
#' cont.matrix = cont.matrix, contrasts2test = contrasts2test, anotPackage = anotPackage,
#' outputDir = outputDir, comparison = comparison, Expressions_And_Top = TRUE ,
#' showParams = FALSE , use.dupCorr = FALSE, block = NULL, nDups = 1 , ENTREZIDs = ENTREZIDs, 
#' SYMBOLIDs = SYMBOLIDs, linksFile = linksFile,fitFileName = fitFileName , csvType=csvType, 
#' rows2HTML = NULL, anotFileName = anotFileName)
#' }
#' @export

lmAnalysis <- function(exprs.filtered,
                       design, cont.matrix,
                       contrasts2test = 1:ncol(cont.matrix),
                       anotPackage,
                       Expressions_And_Top = TRUE,
                       showParams = FALSE,
                       use.dupCorr = FALSE,
                       block = NULL, nDups = 1,
                       comparison = "",
                       outputDir = ".",
                       ENTREZIDs = NULL,
                       SYMBOLIDs = NULL,
                       linksFile,
                       fitFileName,
                       csvType,
                       rows2HTML = NULL,
                       anotFileName)
{
  categLabel ='ANALYSIS'

  ### 1. Ajust del model lineal

  if (use.dupCorr && (!is.null(block)))
  {
    stopifnot(require(statmod))

    corfit <- duplicateCorrelation(exprs.filtered, ndups = nDups, block = blocs, design = design)
    fit < -lmFit(exprs.filtered, design = design, block = blocs, cor = corfit$consensus)
  }else{
    fit <- lmFit(exprs.filtered, design)
  }

  fit.main <- contrasts.fit(fit, cont.matrix)
  fit.main <- eBayes(fit.main)

  ### Taula resum dels resultats : Pendent: numGenesChanged

  numGenesChanged <- genesSelected(NULL)
  if ( (is.null(ENTREZIDs)) & (!is.null(anotPackage))){
    stopifnot(require(old2db (anotPackage), character.only = T))
    envirName <- paste(anotPackage, "ENTREZID", sep = "")
    myenvirENTREZID <- eval(parse(text = envirName))
    ENTREZIDs <- unlist(AnnotationDbi::mget(rownames(exprs.filtered), myenvirENTREZID, ifnotfound = NA))
  }

  if ( (is.null(SYMBOLIDs))  & (!is.null(anotPackage)))
  {
    stopifnot(require(old2db (anotPackage), character.only = T))
    myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
    SYMBOLIDs <- unlist(AnnotationDbi::mget(rownames(exprs.filtered), env = myenvirSYMBOL, ifnotfound = NA))
  }

  for (i in contrasts2test){
    top.Diff <- topTable (fit.main, coef=i, n=nrow(fit.main$t), adjust="fdr")
    fitCoefName<- colnames(cont.matrix)[i]
    contrastTitle <- paste(comparison,colnames(cont.matrix)[i],sep=".")
    atitle<-paste("Genes analyzed in comparison",contrastTitle,sep=" ")
    aFName0<-paste(contrastTitle,"html",sep=".")
    aFName<-paste(file.path(outputDir,contrastTitle), "html",sep=".")

    if((!is.null(anotPackage)) |((!is.null(SYMBOLIDs)) & (!is.null(ENTREZIDs))) )
    {
      if(is.null(rows2HTML))
      {
        top.Diff2HTML <- top.Diff
      }else{
        len.rows2HTML <- length(rows2HTML)
        if (len.rows2HTML==1)
        {
          cat(paste("Only", rows2HTML, "rows selected in html gene lists... \n", sep = " "))
          top.Diff2HTML <- top.Diff[1:rows2HTML, ]
        }else{
          cat(paste("Only",length(rows2HTML), "rows selected in html gene lists... \n", sep = " "))
          top.Diff2HTML <- top.Diff[rows2HTML, ]
        }
      }

      outNUL <- annotateTopTable2(top.Diff2HTML,
                                  aFName,
                                  atitle,
                                  anotPackage=anotPackage,
                                  EntrezIDs=ENTREZIDs,
                                  SymbolIDs=SYMBOLIDs,
                                  anotFilename=anotFileName
      )
    }

    addToLinksFile (linksFile, aFName0, categ=categLabel,
                    desc = paste("Test parameters and ranked list of genes for the comparison", contrastTitle, sep=" "))

    numGenesChanged <- cbind(numGenesChanged, genesSelectable(top.Diff, 0.01, 0.05, 0.25, 0.01, 0.05))


    vFName0<-paste("volcano", contrastTitle, "pdf",sep=".")
    if (toTIFF == TRUE){
      vFName<-file.path(outputDir,paste("volcano", contrastTitle, "tiff",sep="."))
      tiff(file = vFName, width = 3200, height = 3200, units = "px", res = 800)
    }else{
      vFName<-file.path(outputDir,paste("volcano", contrastTitle, "pdf",sep="."))
      pdf(vFName)
    }
    opt <- par(cex.lab = 0.7)
    if(!is.null(SYMBOLIDs)) {
      volcanoNames<-SYMBOLIDs[rownames(exprs.filtered)]
    }else {
      volcanoNames<-rownames(exprs.filtered)
    }

    volcanoplot(fit.main, coef=i, highlight=10, names=volcanoNames,
                main=paste("Differentially expressed genes",contrastTitle, sep="\n"))
    abline(v=c(-1,1))
    par(opt)
    dev.off()
    addToLinksFile (linksFile, vFName0, categ = "ANALYSIS",
                    desc = paste("Volcano Plot for the comparison", contrastTitle, sep=" "))

    if (Expressions_And_Top)
    {
      if (is.twoSampleContrast(cont.matrix[,i])) { # Si la comparacio inclou 2 mostres s'ha de definir grup1 i grup2
        nom1 <- rownames(cont.matrix)[cont.matrix[,i]==1]
        nom2 <- rownames(cont.matrix)[cont.matrix[,i]==-1]
        grup1 <- which(design[, colnames(design) == nom1] == 1)
        grup2 <- which(design[, colnames(design) == nom2] == 1)
      } else if (is.oneSampleContrast(cont.matrix[,i])) { # si inclou sols 1 mostra s'ha de definir sols grup1
        nom1 <- rownames(cont.matrix)[cont.matrix[,i]==1]
        nom2 <- NULL
        grup1 <- which(design[, colnames(design) == nom1] == 1)
        grup2 <- NULL
      } else { # si la comparacio es mes complexa s'han de posar ambdos grups i els seus noms a NULL
        grup1 <- NULL
        grup2 <- NULL
        nom1 <- NULL
        nom2 <- NULL
      }

      topFileName <- paste("ExpressAndTop", contrastTitle, sep=".")
      csvType <- ifelse(is.null(csvType), "csv2", csvType)
      if (showParams) {
        combined <- escriuTop_i_Express(exprs.filtered, top.Diff,
                                        grup1 = grup1, grup2 = grup2, nom1 = nom1, nom2 = nom2,
                                        #fName=file.path(outputDir, topFileName),
                                        fName = topFileName,
                                        fit = fit.main,
                                        fitCoef = fitCoefName,
                                        anotPackage = anotPackage,
                                        my.symbols = SYMBOLIDs,
                                        my.entrezs = ENTREZIDs,
                                        csvType = csvType,
                                        outputDir = outputDir)
      }else {
        combined <- escriuTop_i_Express(exprs.filtered, top.Diff,
                                        grup1 = grup1, grup2 = grup2, nom1 = nom1, nom2 = nom2,
                                        #fName=file.path(outputDir, topFileName),
                                        fName = topFileName,
                                        fit = NULL,
                                        fitCoef = NULL,
                                        anotPackage = anotPackage,
                                        my.symbols = SYMBOLIDs,
                                        my.entrezs = ENTREZIDs,
                                        csvType = csvType,
                                        outputDir = outputDir)
      }
      addToLinksFile (linksFile, paste(topFileName, substr(csvType, 1, 3), sep="."), categ=categLabel,
                      desc = paste("TopTable and normalized expression values for the comparison", fitCoefName, sep=" "))
    }
  }

  colnames(numGenesChanged) <- colnames(cont.matrix)[contrasts2test]
  chgdFName0 <- paste("numGenesChanged", comparison, substr(csvType, 1, 3), sep = ".")
  chgdFName  <- paste("numGenesChanged", comparison,  sep = ".")
  write2csv(numGenesChanged, fileName = chgdFName, csv = csvType, outputDir = outputDir)

  addToLinksFile (linksFile, chgdFName0, categ = categLabel, desc = paste("Number of selectable genes in comparisons", comparison, sep = " "))
  save(fit.main, file = file.path(outputDir, fitFileName))

  return(fit.main)
}
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
### doLmAnalysis
############################################################

#' doLmAnalysis
#'
#'
#' @param lmPar list object that contains the parameters
#' @export

doLmAnalysis <- function (lmPar)
{
  p <- lmPar[[1]]
  if(!is.null(p$ENTREZIDs)){
    EntrezIDs <-  eval(parse(text = p$ENTREZIDs))
  }

  if(!is.null(p$SYMBOLIDs)){
    SymbolIDs <-  eval(parse(text = p$SYMBOLIDs))
  }

  if (!is.null(p$expresFileName)){
    expres <- loadFromFile (file.path(p$outputDir, p$expresFileName))
  }else{
    if (!is.null(p$dades))
    {
      expres <- eval(parse(text = p$dades)) # Posar-hi un tryCatch per poder sortir si d error!!!
    }else{
      stop("Error, Cal definir o les dades o el nom de l'arxiu")
    }
  }

  if (is.null(p$whichContrasts))
  {
    contrasts2test <- 1:ncol(p$contMat)
  }else{
    contrasts2test <- p$whichContrasts
  }

  fitMain <- lmAnalysis(exprs.filtered = expres,
                        design = p$designMat,
                        cont.matrix = p$contMat,
                        contrasts2test = contrasts2test,
                        anotPackage = p$anotPack,
                        outputDir = p$outputDir,
                        comparison = p$comparisonName,
                        Expressions_And_Top = p$ExpressionsAndTop,
                        showParams = p$showLmParams,
                        use.dupCorr = p$use.dupCorr,
                        block = p$block,
                        nDups = p$nDups ,
                        ENTREZIDs = EntrezIDs,
                        SYMBOLIDs = SymbolIDs,
                        linksFile = p$fileOfLinks,
                        fitFileName = p$fitFileName,
                        csvType=p$csvType,
                        rows2HTML = p$rows2HTML,
                        anotFileName = p$anotFilename
  )


  designMatrixName = paste("designMatrix",p$comparisonName, sep=".")
  contrastMatrixName = paste("contrastMatrix",p$comparisonName, sep=".")

  write2csv(p$designMat, fileName = designMatrixName, csv = p$csvType, outputDir = p$outputDir)
  write2csv(p$contMat, fileName = contrastMatrixName, csv = p$csvType, outputDir = p$outputDir)

  csvType <- ifelse(is.null(p$csvType), "csv2", p$csvType)
  addToLinksFile(p$fileOfLinks, paste(designMatrixName, substr(csvType, 1, 3),  sep="."), categ = "ANALYSIS",
                 desc = paste("Design Matrix for comparison ", p$comparisonName, sep=""))
  addToLinksFile(p$fileOfLinks, paste(contrastMatrixName, substr(csvType, 1, 3),  sep="."), categ = "ANALYSIS",
                 desc = paste("Contrast Matrix for comparison ", p$comparisonName, sep=""))

  return (fitMain)
}
