####################################################################

loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}

#####################################################################
old2db <- function(anot){paste(anot, "db", sep = ".")}
#####################################################################
creaAnotFromChipPackage <- function (chipPackage, field="ENTREZ", cleanNAs=T,
                                     isControl= FALSE, ctlCode=NA, removeControls=FALSE)
{
  if (!require(old2db (chipPackage), character.only=T)){
    stop(paste("Required annotation package", chipPackage," is missing"))
  }else{
    require(old2db (chipPackage), character.only=T) # No deu caldre
  }
  if (isControl){
    cleanNAs<- FALSE
    field <- "ACCNUM"
  }
  name2extract <- paste(chipPackage, field, sep="")
  x<- eval(parse(text= name2extract))
  myAnot <- toTable(x) # Aqui hauria de venir un control d'error
  # per si se li ha donat el nom malament
  myAnotTable <- myAnot[,2]
  names(myAnotTable) <- myAnot[,1]
  if (!cleanNAs){
    name2extract <- paste(chipPackage, "ACCNUM", sep="")
    x<- eval(parse(text= name2extract))
    myAccTable <- toTable(x)
    myNAsNames <- setdiff(myAccTable[,1], myAnot[,1])
    # Compte! : Aquesta diferencia es assimetrica
    myNAs <- rep(NA, length(myNAsNames))
    names(myNAs) <- myNAsNames
    myAnotTable<- c(myAnotTable, myNAs)
  }
  if (isControl){
    if (is.na(ctlCode))
      stop("Control codes must be different from NA")
    controls <-sapply (myAnotTable, function (x) ifelse( regexpr(ctlCode, x) > 0, TRUE, FALSE ))
    if (!removeControls){
      myAnotTable<- myAnotTable[controls]
    }else{
      myAnotTable<- myAnotTable[!controls]
    }
  }
  return(myAnotTable)
}
###############################################################
creaAnotFromPDPackage <- function (dbPackage, field, fieldName=NULL, cleanNAs=T,
                                   multipleIDsSymbol=" /// ", removeMultipleIDs=T,
                                   removeControls=TRUE)
{
  if (!require(dbPackage, character.only=T)){
    stop(paste("Required Platfform design package", dbPackage," is missing"))
  }else{
    require(dbPackage, character.only=T) # No deu caldre
  }
  conn<- db(eval(parse(text=dbPackage)))
  fSetType <- dbGetQuery(conn,
                         paste("SELECT DISTINCT meta_fsetid as transcript_id, type_id",
                               "FROM featureSet, core_mps, type_dict",
                               "WHERE featureSet.fsetid=core_mps.fsetid",
                               "AND featureSet.type=type_dict.type"))
  allExceptControls <- as.character(fSetType$transcript_id[fSetType$type_id=="main"])
  allControls <- as.character(fSetType$transcript_id[fSetType$type_id!="main"])

  # Un cop fet aixo podem mirar de recuperar els ids fent servir getNetAffx
  # getNetAffx necessita un ExpressionSet (suposa que ja ha normalitzat)
  # El fem a ma posant-li els noms dels transcripts en el camp assayData

  if  (removeControls){
    transcriptIds <- allExceptControls
  }else{
    transcriptIds <- allControls
    field <- 1 # Es l'ACCNUM que tambe existeix pels controls
    removeMultipleIDs <- FALSE
  }

  numTranscripts <-length(transcriptIds)
  exp <-matrix(NA, nrow=numTranscripts, ncol=1,
               dimnames=list(transcriptIds,"Empty"))
  rownames(exp) <- transcriptIds
  aD <-assayDataNew(exprs=exp)
  nulEset <-new("ExpressionSet",  assayData =aD)
  featureNames(nulEset) <- transcriptIds
  annotation(nulEset) <-dbPackage

  # Ara invoco getNetAffx per recuperar les anotacions

  featureData(nulEset) <- getNetAffx (nulEset, "transcript")
  geneNames <-pData(featureData(nulEset))$geneassignment

  # Creem la taula d'anotacions

  anotTable <- data.frame(transcriptIds,
                          rep(NA, length(transcriptIds))
  )
  if (is.null(fieldName))
    fieldName<- switch(field,
                       "Accession","GeneSymbol", "Gene Title", "Cytoband","Entrez")
  colnames (anotTable) <- c("transcriptIds", fieldName)

  myAnotTable <-geneNames
  names (myAnotTable)<-transcriptIds

  stopifnot(require(gdata))
  for (i in 1:length(transcriptIds)){
    if (!is.na(geneNames[i])){
      l1<- strsplit(geneNames[i],"///")
      if (runMulticore ==1 || runMulticore ==3) {
        l2 <- mclapply(l1, function(l) strsplit(l, "//"))
        myIds <- unique(unlist(mclapply(l2[[1]], function (x) try(trim(x[[field]])))))
      } else {
        l2 <- lapply(l1, function(l) strsplit(l, "//"))
        myIds <- unique(unlist(lapply(l2[[1]], function (x) try(trim(x[[field]])))))
      }


      anotTable[i,2] <- paste(myIds, collapse="//")
    }
  }
  if(removeMultipleIDs){
    s1<-sapply(anotTable[,2], function(s) strsplit(s, "//"))
    if (runMulticore ==1 || runMulticore ==3) {
      s2<-mclapply (s1, function(s) return(s[1]))
    } else {
      s2<-lapply (s1, function(s) return(s[1]))
    }

    myAnotTable <-unlist(s2)
    names(myAnotTable)<-transcriptIds
  }else{
    myAnotTable <-anotTable[,2]
    names(myAnotTable)<-transcriptIds
  }
  return(myAnotTable)
}

#################################################################

#' createOrLoadAnnotations
#'
#' Function that creates the annotation needed, or if it is already done, loads the file created.
#' The annotation packages used need to be installed previously.
#' @param loadAnnotations FALSE by default. If TRUE the function loads the annotations files.
#' @param chipPackAvailable TRUE if there is a chip package available.
#' @param platformDesignPackAvailable  TRUE if there is a platform design package available.
#' @param chipPackage Name of the chip package. Only if chipPackAvailable is TRUE.
#' @param platformDesignPackage Name of the platform design. Only if platformDesignPackAvailable is TRUE.
#' @param outputDir Path where the annotation will be stored.
#' @param annotationsFileName Name of the file for annotations.
#' @param entrezTableFileName Name of the file for Entrez genes.
#' @param symbolsTableFileName Name of the file for gene symbols ID.
#' @param controlsTableFileName Name of the file for controls.
#' @return A list "anotacions" that contains the annotations, the Entrez genes and the gene symbols ID.
#' @examples
#' \dontrun{
#' loadAnnotations <- FALSE
#' chipPackAvailable <- TRUE
#' platformDesignPackAvailable <- FALSE
#' chipPackage <- "hgu133a2"
#' platformDesignPackage <- NULL
#' outputDir <- "./ResultsDir"
#' annotationsFileName <- "Annotations"
#' entrezTableFileName <-"Entrezs.Rda"
#' symbolsTableFileName <-"Symbols.Rda"
#' controlsTableFileName <- "controls.Rda"
#'
#'
#' anotacions <- createOrLoadAnnotations (loadAnnotations= loadAnnotations,
#' chipPackAvailable = chipPackAvailable, platformDesignPackAvailable = platformDesignPackAvailable,
#' chipPackage = chipPackage, platformDesignPackage = platformDesignPackage,
#' outputDir = outputDir,annotationsFileName = annotationsFileName,
#' entrezTableFileName = entrezTableFileName, symbolsTableFileName = symbolsTableFileName,
#'  controlsTableFileName = controlsTableFileName)}
#' @export

createOrLoadAnnotations <-function (loadAnnotations=FALSE,
                                    chipPackAvailable,
                                    platformDesignPackAvailable,
                                    chipPackage,
                                    platformDesignPackage,
                                    outputDir,
                                    annotationsFileName,
                                    entrezTableFileName,
                                    symbolsTableFileName,
                                    controlsTableFileName)
{
  if(!loadAnnotations){
    if (chipPackAvailable){
      entrezTable  <-creaAnotFromChipPackage (chipPackage=chipPackage,
                                              field='ENTREZID')
      symbolsTable <-creaAnotFromChipPackage (chipPackage=chipPackage,
                                              field='SYMBOL')
      controlsTable <-creaAnotFromChipPackage (chipPackage =chipPackage,
                                               field='ACCNUM',
                                               isControl= TRUE,
                                               ctlCode='AFFX',
                                               removeControls=FALSE)
    }else{
      if (platformDesignPackAvailable){
        entrezTable  <-creaAnotFromPDPackage (dbPackage=platformDesignPackage,
                                              field=5)
        symbolsTable <-creaAnotFromPDPackage (dbPackage=platformDesignPackage,
                                              field=2,
                                              removeMultipleIDs=F)
        controlsTable <-creaAnotFromPDPackage (dbPackage=platformDesignPackage,
                                               field=1,
                                               removeControls=F)
      }else{
        symbolsTable <- NULL          # Si hi ha chipPackage ha de ser NULL
        entrezTable <- NULL           # Si hi ha chipPackage ha de ser NULL
        controlsTable <- NULL
      }
    }
  }else{
    entrezTable <-loadFromFile(file=file.path(outputDir, entrezTableFileName))
    symbolsTable <-loadFromFile(file=file.path(outputDir, symbolsTableFileName))
    controlsTable <- loadFromFile(file=file.path(outputDir, controlsTableFileName))
  }
  if(!is.null(entrezTable)){
    save(entrezTable, file=file.path(outputDir, entrezTableFileName)) }
  if(!is.null(symbolsTable)){
    save(symbolsTable, file=file.path(outputDir, symbolsTableFileName))    }
  if(!is.null(controlsTable)){
    save(controlsTable, file=file.path(outputDir, controlsTableFileName))  }
  if (is.null(annotationsFileName)){
    annotationsFileName <- paste("Annotations", "txt", sep=".")
    write.table(data.frame(Entrez=entrezTable, Symbols=symbolsTable),
                sep="\t", quote=FALSE, file=file.path(outputDir,annotationsFileName))  }
  anotacions=list(Entrez=entrezTable,Symbols=symbolsTable, controls=controlsTable)
  return(anotacions)
}
