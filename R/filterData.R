################################################################
sdf <- function (x, threshold) if (sd(x) < threshold) return(FALSE) else return(TRUE)
#################################################################
filt.by.Signal <- function(x, grups, threshold)
{
  if( max(sapply(split(x, grups), mean), na.rm=TRUE) < threshold) return(FALSE) else return(TRUE)
}

#################################################################
filterBySignal <- function(expres, groups, threshold, sigFun.Name="filt.by.Signal", thr.as.perc=TRUE)
{
  if (thr.as.perc){
    signalThr <- quantile(as.vector(expres), threshold / 100)
  }else{
    signalThr <- threshold
  }
  sigFun=eval(parse(text=sigFun.Name))
  filtered.vals <- apply (expres, 1, sigFun, groups, signalThr)
  expres.filtered <- expres[filtered.vals, ]
  return(expres.filtered)
}
##################################################################
filterByVar <- function (expres, threshold, varFun.Name="sdf", thr.as.perc=TRUE, thr.fun=sd)
{
  if (thr.as.perc){
    SD <- apply(expres, 1, thr.fun)
    variabilityThr <- quantile (SD,  threshold / 100) # SD3Q
  }else{
    variabilityThr <- threshold
  }
  varFun=eval(parse(text=varFun.Name))
  filtered.vals <- apply (expres, 1, varFun, variabilityThr)
  expres.filtered <- expres[filtered.vals, ]
  return(expres.filtered)
}
##################################################################
createFilteringReportFile <- function(filteringReportFName, outputDir)
{
  write(paste(rep("-", 70), collapse = ""), file = file.path(outputDir, filteringReportFName), append = FALSE)
  linia <- paste("Type Of Filtering", "Threshold", "Genes Number", "Number Of Samples", sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  write(paste(rep("-", 70), collapse = ""), file = file.path(outputDir, filteringReportFName), append = TRUE)
}
##################################################################
filteringReport <- function(n.exprs = c(NA ,NA),
                            n.exprs.filtered.Affy = c(NA ,NA),
                            n.exprs.filtered.BySignal,
                            signalThr,
                            n.exprs.filtered.ByVar = c(NA ,NA),
                            varThr,
                            outputDir,
                            filteringReportFName = "FilteringReport.txt")
{
  if(!file.exists(file.path(outputDir, filteringReportFName))) createFilteringReportFile(filteringReportFName, outputDir)

  linia <- paste("Non-filtering", "---", n.exprs[1], n.exprs[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  linia <- paste("By Affymetrix Controls", "---", n.exprs.filtered.Affy[1], n.exprs.filtered.Affy[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  linia <- paste("By Signal", signalThr, n.exprs.filtered.BySignal[1], n.exprs.filtered.BySignal[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  linia <- paste("By Variability", varThr, n.exprs.filtered.ByVar[1], n.exprs.filtered.ByVar[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  write(paste(rep("-", 35), collapse = " "), file = file.path(outputDir, filteringReportFName), append = TRUE)

}
#################################################################

#' filterData
#'
#' Function that filter data and create a report file.
#'
#' @param expres Data set to be filtered.
#' @param controls Vector that indicates wich observations are controls.
#' @param removeNAs If TRUE the Na's observations will be deleted.
#' @param entrezs Name of the file for Entrez genes.
#' @param bySignal If TRUE the data is filtered by signal value.
#' @param signalThr Threshold for filter the data.
#' @param grups Name of the different groups of the data.
#' @param sigFun.Name ????????
#' @param sigThr.as.perc ???????????
#' @param byVar If TRUE the data is filtered by ?????
#' @param variabilityThr Threshold for filter the data.
#' @param varFun.Name ????????
#' @param pairingFun.Name Name of the function to be imported to pair data
#' @param targets Name of the targets file.
#' @param doReport If TRUE a report is created.
#' @param outputDir Path of the file created.
#' @param filteringReportFName Name of the report file.
#' @return Matrix object with te data filtered
#' @examples 
#' \dontrun{
#' load("./ResultsDir/normalizedData.Rda")
#' repes <- duplicated(exprs(my.norm), MARGIN=1)
#' exprs(my.norm) <- exprs(my.norm)[!repes,]
#' eset_norm <- my.norm
#' load("./ResultsDir/controls.Rda")
#' removeNAs <- TRUE
#' load("./ResultsDir/Entrezs.Rda")
#' entrezs <- entrezTable
#' SignalFilter <- TRUE
#' signalThreshold <- 50
#' signalFilter.Function <- "filt.by.Signal"
#' signalThreshold.as.percentage <- TRUE
#' VarFilter <- TRUE 
#' variabilityThreshold <- 50
#' variability.Function <- "sdf"
#' variabilityThreshold.as.percentage <- TRUE
#' pairing.Function <- NULL
#' my.targets <-read.AnnotatedDataFrame("./celfiles/targets.txt", header = TRUE, row.names = 1)
#' doReport <- TRUE
#' outputDir <- "./ResultsDir"
#' FilteringReportFileName <- "FilteringReport.txt"

#' exprs.filtered <- filterData(expres = exprs(eset_norm),controls = names(controlsTable),
#' removeNAs = TRUE, entrezs = entrezs ,bySignal = SignalFilter,signalThr = signalThreshold,
#' grups = pData(eset_norm)$grupo, sigFun.Name = signalFilter.Function, 
#' sigThr.as.perc = signalThreshold.as.percentage, byVar = VarFilter, variabilityThr = variabilityThreshold,
#' varFun.Name = variability.Function, varThr.as.perc = variabilityThreshold.as.percentage,
#' pairingFun.Name = pairing.Function, targets = my.targets, doReport = doReport, outputDir = outputDir, 
#' filteringReportFName = FilteringReportFileName)}
#' @export

filterData <- function (expres,
                        controls = NULL,
                        removeNAs = FALSE,
                        entrezs = NULL,
                        bySignal = TRUE,
                        signalThr,
                        grups = NULL,
                        sigFun.Name = "filt.by.signal",
                        sigThr.as.perc = TRUE,
                        byVar = TRUE,
                        variabilityThr,
                        varFun.Name = "sdf",
                        varThr.as.perc = TRUE,
                        pairingFun.Name = NULL,
                        targets,
                        doReport = NULL,          # Per  defecte executara la crida a la funcio que fa el report
                        outputDir = NULL,                # Directori a on desa el report del filtrat
                        filteringReportFName = NULL)     # Nom per defecte del report

{
  ### Matrix object is needed
  if (is.data.frame(expres))
    expres<-as.matrix(expres)


  ### Remove controls and NAs
  if(!is.null(controls)){
    expresNotControls <- setdiff(rownames(expres), controls)
    if (removeNAs &&(!(is.null(entrezs)))){
      entrezsNotNAs <- names(entrezs[!is.na(entrezs)]) # List that contains notNAs of entrezs
      notControlsNorNAs <- intersect(expresNotControls, entrezsNotNAs) #IDs of the expresion probes
      exprs.filtered.0 <- expres[notControlsNorNAs,] # Filter data
    }else{ # If there is no Entrez or we don't want to remove NAs
      exprs.filtered.0 <- expres[notControls,] # Expresion without controls
    }
  }else { # If there is no controls
    if (removeNAs &&(!(is.null(entrezs)))){
      entrezsNotNAs <- names(entrezs[!is.na(entrezs)])
      exprs.filtered.0 <- expres[entrezsNotNAs,]
    }else{
      exprs.filtered.0 <- expres # No filtered by NA or controls
    }
  }


  ### Filtratge per senyal
  if(bySignal && (!is.null(grups)))
  {
    exprs.filtered.1 <- filterBySignal(exprs.filtered.0, grups, threshold= signalThr, sigFun.Name=sigFun.Name, thr.as.perc=sigThr.as.perc )
  }else{
    exprs.filtered.1 <- exprs.filtered.0
  }

  ### Aparellament. MOLTES vegades aquest pas no es fa. ?????????????
  ### EL control de si es fa o no el basem en la funcio d'aparellament que es passa com parametre
  ### Si es deixa a NULL el pas s'omet
  if (!is.null(pairingFun.Name)){
    exprs.filtered.2 <- do.call(pairingFun.Name, list(exprs.filtered.1, targets))
  }else{
    exprs.filtered.2 <- exprs.filtered.1
  }

  ### Filtratge per variabilitat
  if (byVar)
  {
    exprs.filtered.3 <- filterByVar (exprs.filtered.2, threshold= variabilityThr, varFun.Name=varFun.Name, thr.as.perc=varThr.as.perc)
  }else{
    exprs.filtered.3 <- exprs.filtered.2
  }

  ### Report del proces de filtratge
  if (!is.null(doReport))
  {
    if (doReport)
    {
      if ((!is.null(filteringReportFName)) && (!is.null(outputDir)))
      {
        filteringReport(n.exprs = dim(expres),
                        n.exprs.filtered.Affy = dim(exprs.filtered.0),
                        n.exprs.filtered.BySignal = dim(exprs.filtered.1), signalThr = signalThr,
                        n.exprs.filtered.ByVar = dim(exprs.filtered.3), varThr = variabilityThr,
                        outputDir = outputDir,
                        filteringReportFName = filteringReportFName)
      }
    }
  }

  return(exprs.filtered.3)
}
