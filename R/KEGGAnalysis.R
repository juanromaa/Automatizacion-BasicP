####################################################
old2db <- function(anot){paste(anot, "db", sep = ".")}
#####################################################
write.htmltable <- function (x, filename, title = "", sortby = NULL, decreasing = TRUE,
                             open = "wt", formatNumeric = function(x) paste(signif(x,
                                                                                   3)))
{
  if (!is.null(sortby)) {
    if (!sortby %in% colnames(x))
      stop(paste("Invalid argument \"sortby\": could not find a column in data frame x with name",
                 sortby))
    soby = x[, sortby]
    if (!is.numeric(soby))
      stop("Invalid argument \"sortby\": column is not numeric")
    x = x[order(soby, decreasing = decreasing), ]
  }
  outfile <- file(paste(filename, ".html", sep = ""), open = open)
  cat("<html>", "<STYLE>", "<!--TD { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 14px;}-->",
      "<!--H1 { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 22px;}-->",
      "</STYLE>", "<head>", paste("<TITLE>", title, "</TITLE>",
                                  sep = ""), "</head>", "<body bgcolor=#ffffff>", file = outfile,
      sep = "\n")
  if (title != "")
    cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n",
        file = outfile, sep = "\n")
  cat("<CENTER> \n", file = outfile)
  cat("<TABLE BORDER=0>", file = outfile, sep = "\n")
  cat("<TR>", file = outfile)
  for (j in 1:ncol(x)) cat("<TD BGCOLOR=\"", c("#e0e0ff", "#d0d0f0")[j%%2 +
                                                                       1], "\"><B>", colnames(x)[j], "</B></TD>\n", sep = "",
                           file = outfile)
  cat("</TR>", file = outfile)
  for (i in 1:nrow(x)) {
    cat("<TR>", file = outfile)
    for (j in 1:ncol(x)) {
      txt = switch(class(x[[j]]), numeric = formatNumeric(x[i,
                                                            j]), as.character(x[i, j]))
      if (length(grep("^http:", txt)) > 0) {
        txt <- sub(";$", "", txt)
        s <- unlist(strsplit(txt, "[/?=]"))
        txt <- paste("<A HREF=\"", txt, "\" TARGET=\"z\">",
                     s[length(s)], "</A>", sep = "")
      }
      cat("<TD BGCOLOR=\"", c("#e0e0ff", "#d0d0f0", "#f0f0ff",
                              "#e0e0f0")[i%%2 * 2 + j%%2 + 1], "\">", txt,
          "</TD>\n", sep = "", file = outfile)
    }
    cat("</TR>", file = outfile)
  }
  cat("</TABLE></CENTER>", "</body>", "</html>", sep = "\n",
      file = outfile)
  close(outfile)
}
#####################################################

GOTerms2Genes.sql <- function(hgResult, anotPackage)
{
  selectedGOTerms <- intersect(names(geneIdUniverse(hgResult)), summary(hgResult)[, 1])
  selectedGO<- geneIdUniverse(hgResult)[selectedGOTerms]
  if (runMulticore ==1 || runMulticore ==3) {
    selectedGenes <- mclapply(selectedGO, function(x) {intersect(geneIds(hgResult),x)})
  } else {
    selectedGenes <- lapply(selectedGO, function(x) {intersect(geneIds(hgResult),x)})
  }

  sql.ENTREZSYMBOL <- "SELECT gene_id, symbol
                       FROM genes, gene_info
                       WHERE genes._id=gene_info._id"

  if (regexpr(".db", anotPackage) < 0)
  {
    genesAnnot <- dbGetQuery(eval(parse(text = paste(anotPackage, "_dbconn()", sep = ""))), sql.ENTREZSYMBOL)
  }else{
    genesAnnot <- dbGetQuery(eval(parse(text = paste(substr(anotPackage, 1, nchar(anotPackage) - 3), "_dbconn()", sep = ""))), sql.ENTREZSYMBOL)
  }

  if (runMulticore ==1 || runMulticore ==3) {
    selectedSymb <- mclapply(selectedGenes, function(x) {genesAnnot[which(unlist(genesAnnot) %in% x), ]})
    selectedSymb <- mclapply(selectedSymb, function(x) {x <- x[, -1]})
  } else {
    selectedSymb <- lapply(selectedGenes, function(x) {genesAnnot[which(unlist(genesAnnot) %in% x), ]})
    selectedSymb <- lapply(selectedSymb, function(x) {x <- x[, -1]})
  }

  return(selectedSymb)
}

#####################################################
KEGGEnrAn <- function(EntrezIDs,
                      anotFName,
                      anotPackage,
                      organisme = "hsa",
                      outputDir,
                      universe,
                      pval,
                      min.count,
                      addGeneNames,
                      contrast.Name)
{
  universe <- unique(universe)
  EntrezIDs <- unique(EntrezIDs)

  param <- new("KEGGHyperGParams",
               geneIds = EntrezIDs,
               universeGeneIds = universe,
               annotation = anotPackage,
               pvalueCutoff = pval)

  anotPack <- annotation(param)

  getSymbol <- function (x)
  {
    if (length(x)>0)
    {
      simbols <- getSYMBOL(x, anotPack)
    }else{
      simbols <- NULL
    }

    return(simbols)
  }

  hyperRes <- hyperGTest(param)
  sumari <- summary(hyperRes, p=param@pvalueCutoff)
  fName <- paste("KEGG Enrichment Analysis", contrast.Name, sep = " ") # Informe en HTML

  if (addGeneNames)
  {
    EnrichedKEGGTerms <- as.character(sumari[, 1])

    if (length(EnrichedKEGGTerms) > 0)
    {
      selectedSymbols <- GOTerms2Genes.sql(hyperRes, anotPackage)
      genesInHg <- sapply(selectedSymbols, function(x) paste(x, collapse = ", "))
      Report <- cbind(sumari, GeneNames = genesInHg)
    }else{
      Report <- sumari
    }
  }else{
    Report <- sumari
  }

  Report[, 1] <- paste("<a href=\"http://www.genome.jp/dbget-bin/www_bget?path:", organisme, Report[, 1], "\">", Report[, 1], "</a>", sep="")

  ReportSig <- Report[1:nrow(sumari),]

  write.htmltable(x = Report[1:nrow(sumari),],
                  file = file.path(outputDir, anotFName),
                  title =  paste("KEGG Enrichment Analysis", contrast.Name, sep = " "),)
  sortable.html.table(df = Report[1:nrow(sumari),],
                      output.file = paste0(anotFName, "-sortable.html"),
                      output.directory = outputDir,
                      page.title = "KEGG Enrichment Analysis")

}
#####################################################
#'KEGGAnalysis
#'
#' Function to make KEGG analysis.
#' @param fitMain ????
#' @param whichContrasts ????
#' @param comparison.Name ????
#' @param outputDir ????
#' @param anotPackage ????
#' @param organisme ????
#' @param my.IDs ????
#' @param addGeneNames ????
#' @param fileOfLinks ????
#' @param cutoffMethod ????
#' @param P.Value.cutoff ????
#' @param pval ????
#' @param thrLogFC ????
#' @param minNumGens ????
#' @importFrom limma topTable
#' @importFrom annotate getSYMBOL
#' @importFrom GOstats hyperGTest
#' @importFrom SortableHTMLTables sortable.html.table
#' @importFrom DBI dbGetQuery
#' @examples  
#' \dontrun{
#' load("./ResultsDir/fit.Rda")
#' wCont <- 1:3
#' comparison.Name <-"Estudi"
#' outputDir <- "./ResultsDir"
#' orgPackage <-"org.Hs.eg"
#' organisme <-"hsa"
#' load("./ResultsDir/Entrezs.Rda")
#' addGeneNames <- TRUE
#' linksFileName <- "Links.txt"
#' cutoffMethod <-"unadjusted"
#' pval <- 0.05
#' minNumGens <-0
#' runMulticore <- 0
#' 
#' KEGGResult <- BasicP::KEGGAnalysis(fitMain = fit.main, whichContrasts = wCont, 
#' comparison.Name = comparison.Name, outputDir = outputDir, anotPackage = orgPackage, 
#' organisme = organisme, my.IDs = entrezTable, addGeneNames = addGeneNames, 
#' fileOfLinks = linksFileName, cutoffMethod = cutoffMethod, 
#' P.Value.cutoff = rep(0.05, length(wCont)), pval = pval,thrLogFC = 1, 
#' minNumGens = minNumGens)
#' }
#' @export

KEGGAnalysis <- function(fitMain,
                         whichContrasts,
                         comparison.Name,
                         outputDir,
                         anotPackage,
                         organisme,
                         my.IDs,
                         addGeneNames = TRUE,
                         fileOfLinks,
                         cutoffMethod = c("adjusted", "unadjusted"),
                         P.Value.cutoff = rep(0.05, length(whichContrasts)),
                         pval = 0.05,
                         thrLogFC = NULL,
                         minNumGens = 0)
{
  stopifnot(require(old2db(anotPackage), character.only = TRUE))

  categLabel <- "KEGG"

  if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[, whichContrasts])))
  {
    if (is.null(my.IDs))
    {
      myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
      geneUniverse <- unlist(mget(unique( fitMain$genes[,1]), env = myenvirENTREZID, ifnotfound = NA))
    }else{
      geneUniverse <- unlist(my.IDs)
      geneUniverse <- unique(geneUniverse[geneUniverse != "---"])
    }

    for (i in whichContrasts)
    {
      if(is.null(thrLogFC))
      {
        thrLogFC <- 0
      }else{
        thrLogFC <- abs(thrLogFC)
      }

      top.Diff <- topTable(fitMain, coef = i, n = nrow(fitMain$t), adjust = "fdr", lfc = thrLogFC)

      if (cutoffMethod=="adjusted")
      {
        top.Diff.selected <- top.Diff[top.Diff$adj.P.Val < P.Value.cutoff[i], ]
      }else{
        top.Diff.selected <- top.Diff[top.Diff$P.Value < P.Value.cutoff[i], ]
      }

      if (!is.null(top.Diff.selected$ID)){
        gNames <- top.Diff.selected.up$ID
      }else{
        gNames <- rownames(top.Diff.selected)
      }


      if (!(is.null(my.IDs)))
      {
        selectedEntrezIds <- unlist(my.IDs[gNames])
        selectedEntrezIds <- unique(selectedEntrezIds[selectedEntrezIds != "---"])
      }else{
        myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
        selectedEntrezIds <- unlist(mget(unique(gNames), env = myenvirENTREZID, ifnotfound = NA))
      }

      contrast.Name <- colnames(fitMain$contrasts)[i]

      if (is.null(selectedEntrezIds)|(length(selectedEntrezIds) <= minNumGens))
      {
        warning(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " in comparison", contrast.Name, sep = " "))
        cat(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " in comparison", contrast.Name, "\n", sep = " "))
      }else{
        eaFName <- paste("SignificantKEGG", comparison.Name, contrast.Name, sep = ".")
        eA <- KEGGEnrAn(EntrezIDs = selectedEntrezIds,
                        anotFName = eaFName,
                        anotPackage = anotPackage,
                        organisme = organisme,
                        outputDir = outputDir,
                        universe = geneUniverse,
                        pval = pval,
                        min.count = minNumGens,
                        addGeneNames = addGeneNames,
                        contrast.Name = paste("for regulated genes in comparison", contrast.Name, sep = ": "))

        addToLinksFile(fileOfLinks,
                       paste(eaFName,"html", sep="."),
                       categ = categLabel,
                       desc = paste("KEGG Enrichment Analysis regulated genes in comparison", contrast.Name, sep = ": "))
      }
    }
  }
  return(NULL)
}
