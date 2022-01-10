##########################################################################
plotPCA2D <- function(my.data, pc.importance, my.names, my.colors, comp = c(1, 2), posText = 4)
{
  stopifnot(require(rgl))

  scores <- my.data$x
  rownames(scores) <- my.names

  Xlims = c(min(scores[, comp[1]]) - ((max(scores[, comp[1]])-min(scores[, comp[1]]))/6), max(scores[, comp[1]]))
  Ylims = c(min(scores[, comp[2]]), max(scores[, comp[2]]))

  plot(scores,
       main = "Principal Components 2D Plot",
       xlab = paste("PC", comp[1], " ", round(pc.importance[2, comp[1]]*100, 1), "%", sep = ""),
       ylab = paste("PC", comp[2], " ", round(pc.importance[2, comp[2]]*100, 1), "%", sep = ""),
       xlim = Xlims, ylim = Ylims, type = "n")

  for(i in 1:dim(scores)[1])
  {
    points(scores[i, comp[1]],
           scores[i, comp[2]],
           pch = 19, col = my.colors[i])
  }

  text(scores[, comp[1]],
       scores[, comp[2]],
       my.names, cex = 0.6, pos = posText)
}
##########################################################################
plotPCA3D <- function(my.data, pc.importance, my.names, my.colors, comp=c(1, 2, 3))
{
  stopifnot(require(rgl))

  scores <- my.data$x
  rownames(scores) <- my.names

  Xlims = c(min(scores[, comp[1]]), max(scores[, comp[1]]))
  Ylims = c(min(scores[, comp[2]]), max(scores[, comp[2]]))
  Zlims = c(min(scores[, comp[3]]), max(scores[, comp[3]]))

  plot3d(scores,
         main="Principal Components 3D Plot",
         xlab=paste("PC", comp[1], " ", round(pc.importance[2, comp[1]]*100, 1), "%", sep = ""),
         ylab=paste("PC", comp[2], " ", round(pc.importance[2, comp[2]]*100, 1), "%", sep = ""),
         zlab=paste("PC", comp[3], " ", round(pc.importance[2, comp[3]]*100, 1), "%", sep = ""),
         xlim = Xlims, ylim = Ylims, zlim = Zlims,
         type = "n", size = 4)

  for(i in 1:dim(scores)[1])
  {
    points3d(scores[i, comp[1]],
             scores[i, comp[2]],
             scores[i, comp[3]],
             size = 4, pch = 19, col = my.colors[i])
  }

  text3d(scores[, comp[1]],
         scores[, comp[2]],
         scores[, comp[3]],
         my.names, cex = 0.6, pos = 4)
}
##########################################################################
doPCAplot <- function(my.data,
                      sampleNames ,
                      my.colors ,
                      my.groups ,
                      my.cex = 0.8,
                      PCAFFile,
                      outputDir,
                      cor = TRUE,
                      comp = c(1, 2),
                      posText = 4,
                      dim3 = FALSE,
                      csv,
                      my.PCAplot = TRUE,
                      x.coord = -100,
                      y.coord = 100,
                      pch = 19)
{
  if (my.PCAplot)
  {
    pc.my.norm <- prcomp(t(exprs(my.data)))
    pc.importance <- summary(pc.my.norm)$importance
    print(round(pc.importance, 3))

    comp <- c(1, 2)
    plotPCA2D(pc.my.norm, pc.importance, sampleNames, my.colors, comp, posText)

    if (dim3)
    {
      comp <- c(1, 2, 3)
      plotPCA3D(pc.my.norm, pc.importance, sampleNames, my.colors, comp)
    }
  }}


##########################################################################
normplots <- function(my.data, sampleNames, my.colors, my.groups, my.method = "average", my.cex = 0.7, posText = 4, dim3 = FALSE, PCAPlots=TRUE, outputDir, csv=csv)
{

  ### Boxplots

  opt <- par(las=2, cex.axis = my.cex)
  boxplot(exprs(my.data), col = my.colors, names = sampleNames, las = 2 , main = "Normalized (RMA) data", cex.axis = my.cex)
  par(opt)

  ### Dendrograma

  opt <- par(las=2, cex = my.cex)
  clust.euclid.average <- hclust(dist(t(exprs(my.data))), method = my.method)
  plot(clust.euclid.average, main = "Hierarchical clustering of samples",  labels = sampleNames, hang = -1, xlab = " ", sub = " ")
  par(opt)

  ### PCA

  doPCAplot(my.data, sampleNames, my.colors, my.groups, my.cex, PCAFFile, outputDir=outputDir, cor, comp, posText=posText,
            dim3=dim3, csv=csv, my.PCAplot=PCAPlots, x.coord = -100, y.coord = 100, pch = 19)

}
#########################################################
#' addToLinksFile
#' @param linksFile Name of the links file.
#' @param aFName Name of the new object to register in the links file.
#' @param categ Category of the new file.
#' @param desc Description of the file.
#' @param subcateg Subcategory of the file.
#' @export



addToLinksFile <- function(linksFile, aFName, categ = "", desc = "", subcateg = "")
{
  if (!is.null(linksFile))
  {
    write(paste(aFName, categ, subcateg, desc, sep = "\t"), file = linksFile, append = TRUE)
  }
}
##########################################################

#' normplots2File
#'
#' normplots2File plots the result of the normalized data and creates a file.
#'
#' @param my.data Normalized data.
#' @param sampleNames Names of the samples of the normalized data.
#' @param my.colors Colores used in the plots.
#' @param my.groups Name of the groups in the plots
#' @param my.method Method used to make the clusters.
#' @param my.cex Cex used in the plots.
#' @param posText Position of the text.
#' @param dim3 If TRUE the PCA plot is in 3D.
#' @param fileName Name of the file created.
#' @param outputDir Path of the file created.
#' @param PCAPlots If TRUE PCA plots will be created.
#' @param csv Indicates the file type.
#' @param lFile Name of the links file.
#' @importFrom links2File addToLinksFile
#' @examples
#' \dontrun{
#' load("./ResultsDir/normalizedData.Rda")
#' repes <- duplicated(exprs(my.norm), MARGIN=1)
#' exprs(my.norm) <- exprs(my.norm)[!repes,]
#' eset_norm <- my.norm
#' my.colors <- rainbow(length(sampleNames(eset_norm)))
#' my.names <- pData(eset_norm)$ShortName
#' myCex<- 0.8
#' dim3 <- FALSE
#' fileName <- "NormalizedPlots.pdf"
#' outputDir <- "./ResultsDir"
#' PCAPlots <- TRUE
#' csv <- "csv2"

#' normplots2File(my.data = eset_norm, sampleNames = my.names, my.colors = my.colors,
#' my.groups = pData(eset_norm)$Group, my.method = "average",my.cex = myCex ,
#' posText = 2, dim3 = FALSE,fileName = fileName, outputDir = outputDir,PCAPlots = TRUE,
#' csv = fileType)}
#' @export

normplots2File <- function(my.data,
                           sampleNames,
                           my.colors,
                           my.groups,
                           my.method = "average",
                           my.cex = 0.7,
                           posText = 4,
                           dim3 = FALSE,
                           fileName,
                           outputDir,
                           PCAPlots=TRUE,
                           csv,
                           lFile = NULL)
{
  if(!is.null(fileName)) pdf(file.path(outputDir,fileName))

  normplots (my.data=my.data, sampleNames= sampleNames,my.colors= my.colors,
             my.groups= my.groups, my.method = my.method, my.cex = my.cex, posText = posText, dim3 = dim3, PCAPlots=PCAPlots, outputDir=outputDir)

  if(!is.null(fileName)) dev.off() # No entiendo este if.

  addToLinksFile(linksFile=lFile, fileName, categ = 'QC', desc = "Plots of normalized data")
}
