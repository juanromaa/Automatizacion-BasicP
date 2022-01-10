###########################################################
trim <- function(x)
{
  gsub("^[[:space:]]+|[[:space:]]+$", "", x)
}
###########################################################
repofun <- function(x, y, z = NULL, ...)
{
  NAids <- which(is.na(x))
  blnksIDs <- which(x=="&nbsp;")
  x <- trim(x)

  specie <- z
  if(y %in% c("SYMBOL", "GENENAME", "MAP", "ENZYME"))
  {
    out <- x
  }else{
    if(y=="PMID")
    {
      out <- paste("<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/",
                   paste(x, collapse = ","), "?tool=bioconductor\" target=\"_blank\">", length(x), "</a>", sep = "")
    }else{
      out <- switch(y,
                    "ACCNUM" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/protein/", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    #  "ACCNUM" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=nucleotide&term=",
                    #                   x, "[ACCN]&doptcmdl=GenBank\" target=\"_blank\">", x, "</a>", sep = ""),
                    "ENTREZ" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=",
                                     x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "REFSEQ" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/nuccore/", x, "?\" target=\"_blank\">", x, "</a>", sep = ""),
                    "UNIGENE" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=",
                                      substr(x, 4, nchar(x)), "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "ENSEMBL" = paste("<a href=\"http://www.ensembl.org/", z, "/Gene/Summary?g=", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "UNIPROT" = paste("<a href=\"http://www.uniprot.org/uniprot/", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "GO" =  paste("<a href=\"http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "PATH" = paste("<a href=\"http://www.genome.jp/dbget-bin/www_bget?", specie, x, "\" target=\"_blank\">", x, "</a>", sep = ""))
    }
  }

  out[NAids] <- "&nbsp;"
  out[blnksIDs] <- "&nbsp;"

  return(out)
}
###########################################################
translateIDs <- function(geneIDs,
                         anotPackage,
                         typeID = "REFSEQ",
                         toHTML = TRUE,
                         specie = NULL,
                         ...)
{
  annot <- eval(parse(text = paste(anotPackage, typeID, sep = "")))
  newIDs <- AnnotationDbi::mget(geneIDs, annot, ifnotfound = NA)

  if(typeID=="GO")
  {
    if (runMulticore ==1 || runMulticore ==3) {
      goids <- mclapply(1:length(names(newIDs)), function(x, y){names(newIDs[[x]])}, y = newIDs)
    } else {
      goids <- lapply(1:length(names(newIDs)), function(x, y){names(newIDs[[x]])}, y = newIDs)
    }
    names(goids) <- names(newIDs)
    if (runMulticore ==1 || runMulticore ==3) {
      idxs <- which(unlist(mclapply(goids, is.null))==TRUE)
    } else {
      idxs <- which(unlist(lapply(goids, is.null))==TRUE)
    }
    goids[idxs] <- NA
    newIDs <- goids
  }

  if(toHTML)
  {
    if (runMulticore ==1 || runMulticore ==3) {
      newIDs <- mclapply(newIDs, function(x, y, z) {repofun(x, y, z)}, y = typeID, z = specie)
    } else {
      newIDs <- lapply(newIDs, function(x, y, z) {repofun(x, y, z)}, y = typeID, z = specie)
    }

    mysep <- ", "
  }else{
    mysep <- " /// "
  }

  if (runMulticore ==1 || runMulticore ==3) {
    out <- mclapply(newIDs, FUN = function(x, y){paste(x, collapse = mysep)}, y = mysep)
  } else {
    out <- lapply(newIDs, FUN = function(x, y){paste(x, collapse = mysep)}, y = mysep)
  }

  return(out)
}
###########################################################
#'GeneAnnotation
#'
#'Function to annotate genes
#' @param egIDs
#' @param anotPackage Annotation package.
#' @param toHTML If TRUE a html file will be created.
#' @param outputDir Path of the file created.
#' @param filename Name of the file
#' @param myTitle Title of ???
#' @param specie Specie
#' @param info2show
#' @param linksFile Name of the LinksFile
#' @param maxGenes
#' @importFrom SortableHTMLTables sortable.html.table
#' @importFrom links2File addToLinksFile
#' @return It returns the time that the process least.
#' @examples  
#' \dontrun{
#' genes2annotate <- entrezs[unique(rownames(fitMain$p.value))]
#' genesAnnotated <-BasicP::GeneAnnotation(egIDs = genes2annotate, anotPackage = "org.Hs.eg", 
#' toHTML = TRUE, outputDir = outputDir, filename = "Annotations",  
#' myTitle = "Annotations for all genes analyzed", specie = "homo sapiens", 
#' info2show = c( "Affymetrix", "EntrezGene", "GeneSymbol", "GeneName", "KEGG", "GO"), 
#' linksFile = linksFile, maxGenes = NULL)
#' }
#' @export
#'
GeneAnnotation <- function(egIDs,
                           anotPackage,
                           toHTML = TRUE,
                           outputDir,
                           filename = "annotations",
                           myTitle = "Annotations for all genes analyzed",
                           specie = "Homo_sapiens",
                           info2show =  c("Affymetrix", "EntrezGene"),
                           linksFile,
                           maxGenes = NULL,
                           ...)
{
  ptm <- proc.time()

  NAs <- egIDs[is.na(egIDs)]
  eg <- na.omit(egIDs)

  if(toHTML)
  {
    l.ai <- paste("<span id=\"", names(eg), "\">", names(eg), "</span>", sep="")
  }else{
    l.ai <- names(eg)
  }
  l.eg <- repofun(x = eg, y = "ENTREZ")

  a <- cbind(l.ai, l.eg)
  colnames(a) <- c("AffyID", "EntrezGene")

  if("GeneSymbol" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.gs <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "SYMBOL", toHTML))

    a  <- cbind(a, l.gs)
    colnames(a) <- c(colnames.a, "Symbol")
  }
  if("GeneName"  %in% info2show)
  {
    colnames.a <- colnames(a)

    l.gn <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "GENENAME", toHTML))

    a  <- cbind(a, l.gn)
    colnames(a) <- c(colnames.a, "Description")
  }
  if("KEGG" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.kg <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "PATH", toHTML, specie = specie))

    a  <- cbind(a, l.kg)
    colnames(a) <- c(colnames.a, "KEGG")
  }
  if ("GO" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.go <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "GO", toHTML))
    a  <- cbind(a, l.go)
    colnames(a) <- c(colnames.a, "GO")
  }
  if ("PubMed" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.pm <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "PMID", toHTML))

    a  <- cbind(a, l.pm)
    colnames(a) <- c(colnames.a, "Pubmed")
  }

  if("GenBank" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.gb <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "ACCNUM", toHTML))

    a  <- cbind(a, l.gb)
    colnames(a) <- c(colnames.a, "GenBank")
  }
  if("RefSeq" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.rs <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "REFSEQ", toHTML))

    a  <- cbind(a, l.rs)
    colnames(a) <- c(colnames.a, "RefSeq")
  }
  if ("UniGene" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.ug <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "UNIGENE", toHTML))

    a  <- cbind(a, l.ug)
    colnames(a) <- c(colnames.a, "UniGene")
  }
  if("Ensembl" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.em <- translateIDs(geneIDs = eg, anotPackage, typeID = "ENSEMBL", toHTML, specie = specie)

    a  <- cbind(a, l.em)
    colnames(a) <- c(colnames.a, "Ensembl")
  }
  if("UniProt" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.up <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "UNIPROT", toHTML))

    a  <- cbind(a, l.up)
    colnames(a) <- c(colnames.a, "UniProt")
  }
  if("Cytoband" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.cy <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "MAP", toHTML))

    a  <- cbind(a, l.cy)
    colnames(a) <- c(colnames.a, "Cytoband")
  }
  if("Enzyme" %in% info2show)
  {
    colnames.a <- colnames(a)

    l.ez <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "ENZYME", toHTML))

    a  <- cbind(a, l.ez)
    colnames(a) <- c(colnames.a, "Enzyme")
  }

  if(!is.null(maxGenes))
  {
    if (!toHTML) maxGenes <- nrow(a) # si volem anotacions en .txt enlloc de .html no cal tallar l'output
    k <- ceiling(nrow(a) / maxGenes) # tallarem l'output en k fitxers en funcio del maxim d'entrades per fitxer

    j <- 1 # contador dels "marges" per fer els subsets
    for(i in 1:k)
    {
      genAnnot <- a[j:ceiling((nrow(a)*(i/k))), ]

      outFileName <- paste(filename, paste(i, "of", k ,sep = ""), sep = ".")

      if(toHTML)
      {
        sortable.html.table(df = as.data.frame(genAnnot),
                            output.file = paste0(outFileName,"-sortable.html"),
                            output.directory = outputDir,
                            page.title = myTitle )

      }else{
        write.table(x = genAnnot, file = file.path(outputDir, paste(outFileName, "txt", sep = ".")),
                    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
      }

      anotationsFName <- ifelse(toHTML, paste(outFileName, "html", sep="."), paste(outFileName, "txt", sep="."))
      addToLinksFile(linksFile = linksFile,
                     aFName = anotationsFName,
                     categ = 'ANNOT',
                     desc = "Gene annotations for all genes analyzed")
      j <- j + ceiling((nrow(a) / k)) # actualitzem els marges pel seguent subset
    }
  }else{
    genAnnot <- a

    outFileName <- filename

    if(toHTML)
    {
      sortable.html.table(df = as.data.frame(genAnnot),
                          output.file = paste0(outFileName,"-sortable.html"),
                          output.directory = outputDir,
                          page.title = myTitle )
    }else{
      write.table(x = genAnnot, file = file.path(outputDir, paste(outFileName, "txt", sep = ".")),
                  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }

    anotationsFName <- ifelse(toHTML, paste(outFileName, "html", sep="."), paste(outFileName, "txt", sep="."))
    addToLinksFile(linksFile = linksFile,
                   aFName = anotationsFName,
                   categ = 'ANNOT',
                   desc = "Gene annotations for all genes analyzed")
  }

  return(proc.time() - ptm)
}
