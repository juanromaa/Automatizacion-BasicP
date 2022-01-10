###########################################################
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
###########################################################
goLinksTest <- function(my.GOIDs)
{
  if (runMulticore ==1 || runMulticore ==3) {
    GOIDslinked <- unlist(mclapply(my.GOIDs, function(x) {paste("<a href=\"http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&query=", x, "\">",
                                                                x,
                                                                "</a>",
                                                                sep = "")}))
  } else {
    GOIDslinked <- unlist(lapply(my.GOIDs, function(x) {paste("<a href=\"http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&query=", x, "\">",
                                                              x,
                                                              "</a>",
                                                              sep = "")}))
  }

  return(as.character(GOIDslinked))
}
###########################################################
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
###########################################################

enrichment_Analysis <- function(EntrezIDs,
                                anotFName,
                                anotPackage,
                                outputDir,
                                universe = NULL,
                                pval = 0.05,
                                min.count = 3,
                                addGeneNames = TRUE,
                                ontologias = c("MF", "BP", "CC"),
                                testDirections = c("over", "under"),
                                contrast.Name)
{

  universe <- unique(universe)
  EntrezIDs <- unique(EntrezIDs)

  params <- new("GOHyperGParams",
                geneIds = EntrezIDs,
                universeGeneIds = universe,
                annotation = anotPackage,
                ontology = "MF",
                pvalueCutoff = pval,
                testDirection = "over")

  resCum <- NULL; resum <- NULL
  ontoTitle <- c(MF="Molecular Function", BP="Biological Process", CC="Cellular Component")

  for (onto in ontologias)
  {
    params@ontology <- onto

    dir.i <- 1

    for (direction in testDirections)
    {
      params@testDirection <- direction

      hgResult <- hyperGTest(params)

      if (addGeneNames)
      {
        EnrichedGOTerms <- as.character(summary(hgResult)[,1])
        if (length(EnrichedGOTerms) > 0)
        {
          selectedSymbols <- GOTerms2Genes.sql(hgResult, anotPackage)
          genesInHg <- sapply(selectedSymbols, function(x) paste(x, collapse=", "))
          reshgTest <- cbind(summary(hgResult), GeneNames = genesInHg)
        }else{
          reshgTest <- summary(hgResult)
        }
      }else{
        reshgTest <- summary(hgResult)
      }

      colorTest <- ifelse(direction=="over", "red;\">", "green;\">")
      reshgTest <- cbind(reshgTest,
                         OverUnder = rep(paste("<center><span style=\"color:", colorTest, direction, "</span></center>", sep = ""), dim(reshgTest)[1]))
      reshgTest[, 1] <- goLinksTest(reshgTest[, 1])
      colnames(reshgTest)[1] <- "GOID"

      if (dir.i==1)
      {
        reshgTest2table <- reshgTest
      }else{
        reshgTest2table <- rbind(reshgTest2table, reshgTest)
      }

      dir.i <- dir.i + 1
    }

    if (dim(reshgTest2table)[1]!=0)
    {
      reshgTest2table <- cbind(Ontology = paste("<center>", onto, "</center>", sep=""),
                               reshgTest2table)

      resum <- rbind(resum, reshgTest2table)
    }
  }

  if (!is.null(resum))
  {
    if (addGeneNames)
    {
      my.table <- resum[, c(1, 2, 8, 9, 7, 6, 5, 4, 3, 10)]
    }else{
      my.table <- resum[, c(1, 2, 8, 7, 6, 5, 4, 3, 9)]
    }

    write.htmltable(x = my.table,
                    file = file.path(outputDir, anotFName),
                    title =  paste("GO Enrichment Analysis", contrast.Name, sep = " "),
                    open = "wt")
    sortable.html.table(df = my.table,
                        output.file = paste0(anotFName, "-sortable.html"),
                        output.directory = outputDir,
                        page.title = "GO Enrichment Analysis" )
  }

}
###########################################################
#' GOAnalysis
#'
#' Function that do GOAnalysis
#' @param fitMain
#' @param whichContrasts
#' @param comparison.Name Name of the comparison. To identify the output.
#' @param outputDir Path of the file created.
#' @param anotPackage Annotation package.
#' @param my.IDs
#' @param addGeneNames Vector of the gene names.
#' @param fileOfLinks Name of the links file.
#' @param thrLogFC Threshold for the logarithm Fold Change
#' @param cutoffMethod
#' @param P.Value.cutoff Maximum value for the p.value.
#' @param pval
#' @param min.count
#' @param ontologias
#' @param testDirections
#' @param minNumGens Mimimum number of genes.
#' @importfrom topTable limma
#' @importFrom GOstats hyperGTest
#' @importFrom SortableHTMLTables sortable.html.table
#' @importFrom DBI dbGetQuery
#' @importFrom links2File addToLinksFile
#' @import GO.db
#' @import GOstats
#' @examples  
#' \dontrun{
#' library(GOstats)
#' GOResult<-BasicP::GOAnalysis(fitMain = fitMain, whichContrasts = wCont, 
#' comparison.Name = "Estudi", outputDir = outputDir, anotPackage = "org.Hs.eg", 
#' my.IDs = entrezTable, addGeneNames = TRUE, fileOfLinks = linksFile, thrLogFC = 1, 
#' cutoffMethod = "adjusted", P.Value.cutoff = rep(0.05, length(wCont)), pval = 0.01, 
#' min.count = 3, ontologias = c("MF", "BP", "CC"), testDirections = c("over", "under"), 
#' minNumGens = 0)
#' }
#' @export


GOAnalysis <- function(fitMain,
                       whichContrasts,
                       comparison.Name,
                       outputDir,
                       anotPackage,
                       my.IDs,
                       addGeneNames = TRUE,
                       fileOfLinks,
                       thrLogFC=NULL,
                       cutoffMethod = c("adjusted", "unadjusted"),
                       P.Value.cutoff = rep(0.05, length(whichContrasts)),
                       pval = 0.05,
                       min.count = 3,
                       ontologias = c("MF", "BP", "CC"),
                       testDirections = c("over", "under"),
                       minNumGens = 0)
{

  categLabel <- "GO"

  if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[, whichContrasts])))
  {
    if (is.null(my.IDs))
    {
      stopifnot(require(old2db (anotPackage), character.only=T))
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

      top.Diff <- topTable(fitMain, coef = i, n = nrow(fitMain$t), adjust = "fdr", lfc=thrLogFC)  # Seleccionem per FDR y mes


      if (cutoffMethod=="adjusted")
      {
        top.Diff.selected.up  <- top.Diff[(top.Diff$adj.P.Val < P.Value.cutoff[i]) & (top.Diff$t >0), ]
        top.Diff.selected.down <- top.Diff[(top.Diff$adj.P.Val < P.Value.cutoff[i]) & (top.Diff$t < 0), ]
      }else{
        top.Diff.selected.up  <- top.Diff[(top.Diff$P.Value < P.Value.cutoff[i]) & (top.Diff$t >0), ]
        top.Diff.selected.down <- top.Diff[(top.Diff$P.Value < P.Value.cutoff[i]) & (top.Diff$t < 0), ]
      }

      if (!is.null(top.Diff$ID)){

        gNames.up <- top.Diff.selected.up$ID
        gNames.down <- top.Diff.selected.down$ID
      }else{
        gNames.up <- rownames(top.Diff.selected.up)
        gNames.down <- rownames(top.Diff.selected.down)
      }
      if (!(is.null(my.IDs)))
      {
        selectedEntrezIds.up <- unlist(my.IDs[gNames.up])
        selectedEntrezIds.up <- unique(selectedEntrezIds.up[selectedEntrezIds.up != "---"])

        selectedEntrezIds.down <- unlist(my.IDs[gNames.down])
        selectedEntrezIds.down <- unique(selectedEntrezIds.down[selectedEntrezIds.down != "---"])

      }else{

        stopifnot(require(old2db (anotPackage), character.only=T))

        myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))

        selectedEntrezIds.up <- unlist(mget(unique(gNames.up), env = myenvirENTREZID, ifnotfound = NA))
        selectedEntrezIds.down <- unlist(mget(unique(gNames.down), env = myenvirENTREZID, ifnotfound = NA))
      }

      contrast.Name <- colnames(fitMain$contrasts)[i]

      if (is.null(selectedEntrezIds.up)|(length(selectedEntrezIds.up)<=minNumGens))
      {
        warning(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t > 0 in comparison", contrast.Name, sep = " "))
        cat(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t > 0 in comparison", contrast.Name, "\n", sep = " "))
      }else{
        eaUpFName <- paste("SignificantGO", comparison.Name, contrast.Name, "Up", sep = ".")
        eA.up <- enrichment_Analysis(EntrezIDs = selectedEntrezIds.up,
                                     anotFName = eaUpFName,
                                     anotPackage = anotPackage,
                                     outputDir = outputDir,
                                     universe = geneUniverse,
                                     pval = pval,
                                     min.count = min.count,
                                     addGeneNames = addGeneNames,
                                     ontologias = ontologias,
                                     testDirections = testDirections,
                                     contrast.Name = paste("for up-regulated genes in comparison", contrast.Name, sep = ": "))


        addToLinksFile(fileOfLinks,
                       paste(eaUpFName,"html", sep="."),
                       categ = categLabel,
                       desc = paste("Enrichment Analysis for up-regulated genes in comparison", contrast.Name, sep = ": "))


      }

      #      if (is.null(selectedEntrezIds.down))
      if (is.null(selectedEntrezIds.down)|(length(selectedEntrezIds.down)<=minNumGens))
      {
        warning(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t < 0 in comparison", contrast.Name, sep = " "))
        cat(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t < 0 in comparison", contrast.Name, "\n", sep = " "))
      }else{
        eaDownFName <- paste("SignificantGO", comparison.Name, contrast.Name, "Down", sep = ".")
        eA.down <- enrichment_Analysis(EntrezIDs = selectedEntrezIds.down,
                                       anotFName = eaDownFName,
                                       anotPackage = anotPackage,
                                       outputDir = outputDir,
                                       universe = geneUniverse,
                                       pval = pval,
                                       min.count = min.count,
                                       addGeneNames = addGeneNames,
                                       ontologias = ontologias,
                                       testDirections = testDirections,
                                       contrast.Name = paste("for down-regulated genes in comparison", contrast.Name, sep = ": "))

        addToLinksFile(fileOfLinks, paste(eaDownFName, "html", sep="."),
                       categ = categLabel,
                       desc = paste("Enrichment Analysis for down-regulated genes in comparison", contrast.Name, sep = ": "))
      }
    }
  }

  return(NULL)
}
