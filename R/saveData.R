#################################################################
expresAndgeneSymbols <- function(expres, expresNames=colnames(expres),
                                 anotPackage = NULL, SYMBOL="SYMBOL",
                                 symbolsVector = NULL)
{
  if (!is.null(anotPackage))
  {
    my_SYMBOL_env <- eval(parse(text = paste(anotPackage, SYMBOL, sep="")))
    geneSymbols <- unlist(mget(rownames(expres), my_SYMBOL_env, ifnotfound=NA))
    expresWithSymbols <- data.frame(geneSymbols, expres)
    colnames(expresWithSymbols) <- c("Symbol", expresNames)
  }else{
    if (!is.null(symbolsVector)){
      geneSymbols <- symbolsVector[rownames(expres)]
      expresWithSymbols <- data.frame(geneSymbols, expres)
      colnames(expresWithSymbols) <- c("Symbol", expresNames)
    }else{
      expresWithSymbols <- expres
    }
  }
  return(expresWithSymbols)
}
###############################################################
write2csv <- function(my.data, fileName, csv = c("csv2", "csv", "txt", "xls"), outputDir)
{
  fileName<- file.path(outputDir, paste(fileName, substr(csv[1], 1, 3) , sep = "."))
  switch (csv[1],
          "csv" = write.csv(my.data, file = fileName, quote = F),
          "csv2" = write.csv2(my.data, file = fileName, quote = F),
          "txt" = write.table(my.data, file = fileName, quote = F))
  # }
}
##################################################################

#'saveData
#'
#'Function that saves a R object in files.
#'
#' @param expres Object to be saved.
#' @param expresNames Names of the samples of expression.
#' @param expres.csv.FileName Name of the csv file.
#' @param csvType Csv type.
#' @param description Description of the file for Linksfile.
#' @param anotPackage Package of annotations.
#' @param SYMBOL ????????
#' @param symbolsVector Name of the symbols table.
#' @param expres.bin.FileName File name of the filtered data.
#' @param linksFile Character string that indicates the path to the txt file.
#' @param outputDir Path of the file created.
#' @importFrom links2File addToLinksFile
#' @examples
#' \dontrun{
#' load("./ResultsDir/normalizedData.Rda")
#' repes <- duplicated(exprs(my.norm), MARGIN=1)
#' exprs(my.norm) <- exprs(my.norm)[!repes,]
#' eset_norm <- my.norm
#' normalized.all.FileName <- "normalized.all"
#' fileType <-"csv2"
#' symbolsTable <- load("./ResultsDir/Symbols.Rda")
#' expres.all.FileName <- "expres.Rda"
#' linksFileName <- "Links.txt"
#' outputDir <- "./ResultsDir"

#' saveData(expres = exprs(eset_norm), expres.csv.FileName = normalized.all.FileName,
#' csvType=fileType, description = "Normalized values for all genes", anotPackage = NULL,
#' symbolsVector = symbolsTable, SYMBOL = "SYMBOL", expres.bin.FileName = expres.all.FileName,
#' linksFile = linksFileName, outputDir = outputDir)}
#' @export

saveData <- function (expres, expresNames=colnames(expres),
                      expres.csv.FileName, csvType, description="Normalized Values",
                      anotPackage, SYMBOL="SYMBOL", symbolsVector = NULL,
                      expres.bin.FileName,
                      linksFile, outputDir )
{
  if (!(is.null(expres.csv.FileName))){
    expres.all <- expresAndgeneSymbols(expres=expres, expresNames=expresNames,
                                       anotPackage=anotPackage, SYMBOL=SYMBOL, symbolsVector=symbolsVector)
    write2csv(expres.all, fileName = expres.csv.FileName, csv = csvType, outputDir = outputDir)
    addToLinksFile(linksFile, paste(expres.csv.FileName, substr(csvType, 1, 3), sep="."), categ = 'DATA', desc= description)
  }

  if (!(is.null(expres.bin.FileName))){
    save(expres, file=file.path(outputDir, expres.bin.FileName))
  }
}

