#' checkPairing
#' \code{checkPairing} is a simple function to check if two matrices have the same dimensions and the same row and column names.
#' @param X,Y The Two matrices  to be checked.
checkPairing <- function (X, Y){
  check <- (sum(rownames(X)!=rownames(Y)) == 0) && (sum(colnames(X)!=colnames(Y)) == 0) 
}

#' messageTitle
#' \code{messageTitle} is a simple wrapper for sending messages to console
#' @param aMessage The text of the message.
#' @param underChar Character to use for underlining the message.
#' @example
#' messageTitle("Hello world", "$")
#' 
messageTitle <- function(aMessage, underChar="-"){
  cat (aMessage,"\n")
  cat(paste(rep(underChar, nchar(aMessage)),collapse=""),"\n")
}

#' read2
#' \code{read2} is a simple wrapper to read 2 files with the same format, dimensions and row and column names
#' @param expresFName Name of first file, expected to contain expression values.
#' @param metFName Name of second file, expected to contain methylation values. 
#' @param dataDirectory Name of directory where the files are staore. Defaults to ".".
#' @param sepChar Name of character used to separate matrix columns. Defaults to ";".
#' @param decChar Name of character used as decimal point. Defaults to ".".
#' 
read2 <- function (expresFName, metFName, 
                   dataDirectory=".", sepChar=";", decChar= "."){
  expres<- read.table(file=file.path(dataDirectory, expresFName), header=TRUE, 
             sep=sepChar,dec=decChar, row.names = 1)
  mets <-read.table(file=file.path(dataDirectory, metFName), header=TRUE, 
                    sep=sepChar,dec=decChar, row.names = 1) 
  if (checkData(expres, mets))
    result <- list(expres,mets)
  else
    result <- NULL
  return(result)
}

# Functions \texttt{calcFreqs}, \texttt{numScore}, \texttt{binScore}, \texttt{binScore} 
# and \texttt{plotGene} have been created to respectively count, score or plot points 
# on a scatterplot with a 3x3 grid overimposed.

#' calcFreqs
#' \code{calcFreqs} Given two numeric vectors this function overimposes a grid on the scatterplot defined by YMet~Xmet 
#' #' and returns a 3X3 matrix of counts with the number of points in each cell of the grid 
#' #' for given vertical and horizontal lines.
#' @param xMet, First numeric vector which, in principle, contains the methylation values.
#' @param yExp, Second numeric vector which, in principle, contains the expression values.
#' @param x1, x2, Coordinates of vertical points in the X axis. Because it is expected to contain methylation values that vary between 0 and 1 the default values are 1/3 and 2/3.
#' @param y1, y2, Coordinates of vertical points in the Y axis. Leaving them as NULL assigns them the percentiles of yVec defined by `percY1` and `percY2`.
#' @param percY1, percY2 Values used to act as default for `y1`and `y2` when these are set to `NULL`
#' @example 
#' data(trueGenes)
#' xVec<- as.numeric(myMetilData[trueGene1,])
#' yVec<-as.numeric(myExprData[trueGene1,])
#' titolT <- trueGene1
#' plotGenSel(xMet=xVec, yExp=yVec, titleText=titolT, x1=1/3, x2=2/3)
#' messageTitle(paste("Cell counts for gene",trueGene1))
#' calcFreqs(xMet=xVec, yExp=yVec, x1=1/3, x2=2/3, y1=NULL, y2=NULL, percY1=1/3, percY2=2/3)
#'
#' xVec<- as.numeric(myMetilData[falseGene1,])
#' yVec<-as.numeric(myExprData[falseGene1,])
#' titolF <- falseGene1
#' messageTitle(paste("Cell counts for gene",falseGene1))
#' plotGenSel(xMet=xVec, yExp=yVec, titleText=titolF, x1=1/3, x2=2/3)
#' calcFreqs(xMet=xVec, yExp=yVec, x1=1/3, x2=2/3, y1=NULL, y2=NULL, percY1=1/3, percY2=2/3)
#' 
calcFreqs <- function (xMet, yExp, x1, x2, y1=NULL, y2=NULL, 
                      percY1=1/3, percY2=2/3)  
{
  xMet<-as.numeric(xMet)
  yExp<-as.numeric(yExp)
  freqsMat <-matrix(0, nrow=3, ncol=3)
  xVals <- c(x1, x2)
  minExp<-min(yExp); maxExp <- max(yExp); delta<- maxExp-minExp
  if (is.null(y1)) y1<- minExp + percY1*delta
  if (is.null(y2)) y2<- minExp + percY2*delta
  yVals <- c(y1, y2)
  condX <- c("(xMet<=x1)", "((xMet>x1) & (xMet<=x2))", "(xMet>x2)")
  condY <- c("(yExp>y2)", "((yExp<=y2) & (yExp>y1))", "(yExp<=y1)") 
  for (i in 1:3){
    for (j in 1:3){
      condij <- paste(condX[j], condY[i], sep="&")
      freqsMat [i,j] <- sum(eval(parse(text=condij)))
    }
  }
  return(freqsMat)
}

#' toReqMat
#' \code{toReqMat} can be used to turn a matrix of required percentages into 
#' a matrix of required counts to facilitate its scoring
#' @param numPoints Number of points in a scatterplot. Used to turn the required percentages into required counts.
#' @param aReqPercentMat Matrix of required percentages
toReqMat <- function (numPoints, aReqPercentMat){
  return(round(aReqPercentMat*numPoints/100,0)) #  return(round(aReqPercentMat*sum(numPoints)/100,0))
}

#' binScore
#' \code{binScore} can be used to score scatterplots by directly comparing 
#' the sample counts with a matrix of minimal or maximal percentages/counts
#' to be found in each cell. It implements the three bands rule implicitly
#' by setting threshold values
#' @param aGrid A matrix of counts as computed by `calcFreqs` function
#' @param aReq A matrix of minimum or maximum counts to be found in each cell 
#' if L-shape is TRUE
#' @examples 
#' reqPercentages <- matrix (c(15, 5, 0, 0, 5, 5, 10, 10, 15), nrow=3, byrow=TRUE)
#' (countsRnd    <- matrix(floor(runif(9)*10)+1, nrow=3, ncol=3
#' (reqRnd <- toReqMat (sum(countsRnd), reqPercentages)) 
#' binScore(countsRnd, reqRnd)
#' (countsTrueL1 <- matrix (c(20, 3, 0, 10, 2, 2, 20, 10, 20), nrow=3, byrow=TRUE))
#' (reqTrueL <- toReqMat (sum(countsTrueL1), reqPercentages))
#' binScore(countsTrueL1, reqTrueL)
binScore <- function(aGrid, aReq){
  # show(aGrid)
  #   show(aReq)
  #  cat((aGrid[1,1]>= aReq[1,1]), "\t",(aGrid[1,2]<= aReq[1,2]), "\t",(aGrid[1,3]<= aReq[1,3]), "\n", 
  #      (aGrid[2,1]>=aReq[2,1]), "\t",(aGrid[2,2]<=aReq[2,2]), "\t",(aGrid[2,3]<=aReq[2,3]), "\n",
  #      (aGrid[3,1]>=aReq[3,1]), "\t",(aGrid[3,2]>=aReq[3,2]), "\t",(aGrid[3,3]>=aReq[3,3]), "\n")
  
  comp <- (aGrid[1,1]>= aReq[1,1])&&(aGrid[1,2]<= aReq[1,2])&&(aGrid[1,3]<= aReq[1,3]) && 
    (aGrid[2,1]>=aReq[2,1]) &&(aGrid[2,2]<=aReq[2,2]) &&(aGrid[2,3]<=aReq[2,3]) &&
    (aGrid[3,1]>=aReq[3,1])&&(aGrid[3,2]>=aReq[3,2])&&(aGrid[3,3]>=aReq[3,3])
  return(comp)
}

#' numScore
#' \code{numScore} can be used to score scatterplot using a weight matrix. 
#' The scoring does not incorporate logical conditions such as "if xij < C ..."
#' @param aGrid A matrix of counts as computed by `calcFreqs` function
#' @param LShaped A boolean value indicating if this scatterplot can be seen as LShaped. 
#' @param aWeightMifL A matrix of weights to score the previous counts if the scatterplot has been classified as L.
#' @param aWeightMifNonL A matrix of weights to score the previous counts if the scatterplot has been classified as non-L
#' @example 
#' xVecT <- as.numeric(trueLMet[1,]); yVecT<- as.numeric(trueLExpr[1,])
#' xVecF <- as.numeric(falseLMet[1,]); yVecF <- as.numeric(falseLExpr[1,])
#' trueFreq <- calcFreqs(xMet=xVecT, yExp=yVecT, x1=1/3, x2=2/3)
#' falseFreq <- calcFreqs(xMet=xVecF, yExp=yVecF, x1=1/3, x2=2/3)
#' weightsIfL    <- matrix (c(2,-1,-99,1,0,-1,1,1,2), nrow=3, byrow=TRUE)
#' weightsIfNonL <- matrix (c(0,-1,-99,0,0,-1,0,0,0), nrow=3, byrow=TRUE)
#' numScore (trueFreq, weights)
numScore <- function (aGrid, LShaped, aWeightMifL, aWeightMifNonL){
  if (!LShaped){
    scoresM <- aGrid * aWeightMifNonL
  }else{
    scoresM <- aGrid * aWeightMifL
    }
  return(sum(scoresM))
}

#' scoreGenesMat
#' \code{scoreGenesMat} scores scatterplot using a binary and a numeric schemes on a row-wise basis 
#' @param mets Matrix of methylation values 
#' @param expres Matrix of expression values
#' @param aReqPercentsMat Matrix of minimum maximum percentage of counts to have in a given cell
#' @param aWeightMifL A matrix of weights to score the previous counts if the scatterplot has been classified as L.
#' @param aWeightMifNonL A matrix of weights to score the previous counts if the scatterplot has been classified as non-L
#' PARAMETROS QUE HEREDA (O NO) DE CALCFREQS
#' @param x1, x2, Coordinates of vertical points in the X axis. Because it is expected to contain methylation values that vary between 0 and 1 the default values are 1/3 and 2/3.
#' @param y1, y2, Coordinates of vertical points in the Y axis. Leaving them as NULL assigns them the percentiles of yVec defined by `percY1` and `percY2`.
#' @param percY1, percY2 Values used to act as default for `y1`and `y2` when these are set to `NULL`
#' @param trace Set to TRUE to print count matrices as they are computed
standardWeightMifL    <- matrix (c(2,-2,-100,1, 0,-2, 0, 0, 0), nrow=3, byrow=TRUE)
standardWeightMifNonL <- matrix (c(0,-2,-100,0, 0,-2, 0, 0, 0), nrow=3, byrow=TRUE)
standardPercentsMat <- matrix (c(10, 20, 0, 5, 0, 20, 0, 5, 10), nrow=3, byrow=TRUE)

scoreGenesMat <- function(mets, expres,
                           x1=1/3, x2=2/3, 
                           y1=NULL, y2=NULL, percY1=1/3, percY2=2/3,
                           aReqPercentsMat, # =standardPercentsMat,
                           aWeightMifL, # =standardWeightM,
                           aWeightMifNonL) # , trace=FALSE)
{
  N <- dim(mets)[2]
  Ngenes <-nrow(mets)
  scores <- data.frame(logicSc=rep(FALSE, Ngenes), numericSc=rep(0,Ngenes))
  rownames(scores)<- rownames(mets)
  minmaxCounts <- toReqMat (N, aReqPercentMat=aReqPercentsMat)
  for (gene in 1:Ngenes){
    theGene <- rownames(expres)[gene]
    xVec<- mets[theGene,]
    yVec<- expres[theGene,]
    geneGrid <- calcFreqs(xMet=xVec, yExp=yVec, x1=x1, x2=x2,
                          y1=y1, y2=y2, percY1=percY1, percY2=percY2)
    # TRACE
    # if (trace) show(geneGrid)
    # TRACE
    binSc <-  binScore (geneGrid, minmaxCounts)
    scores[gene, "logicSc"] <- binSc
    numSc <- numScore (geneGrid, LShaped=binSc,
                       aWeightMifL=aWeightMifL, 
                       aWeightMifNonL=aWeightMifNonL )
    scores[gene, "numericSc"] <- numSc
  }
  return (scores)
}

#' plotGeneSel
#' \code{PlotGene} plot points on a scatterplot with a 3x3 grid overimposed
#' @example
#' data(trueLGenes)
#' trueGene1 <-rownames(trueLGenes)[1]
#' xVec<- as.numeric(myMetilData[trueGene1,])
#' yVec<-as.numeric(myExprData[trueGene1,])
#' titolT <- paste (trueGene1, "(May be GRM)")
#' plotGeneSel(xMet=xVec, yExp=yVec, titleText=titolT, x1=1/3, x2=2/3)
#'
#'
plotGenSel <- function(xMet, yExp, titleText, 
                         x1=1/3, x2=2/3, y1=NULL, y2=NULL, 
                         percY1=1/3, percY2=2/3, plotGrid=TRUE)  
{
    minExp<-min(yExp); maxExp <- max(yExp); delta<- maxExp-minExp
    plot(xMet,yExp,  xlim=c(0,1), ylim=c(minExp, maxExp), main=titleText)
    if (plotGrid){
      if (is.null(y1)) y1<- minExp + percY1*delta
      if (is.null(y2)) y2<- minExp + percY2*delta
      abline(v=x1);  abline(v=x2)
      abline(h=y1);  abline(h=y2)
    }
}

#' plotGenesMat
#' \code{plotGenesMat} is a simple wrapper for plotting the scatterplots associated with two matrices.
#' @examples
#' plotGenesMat (mets=trueLMet, expres=trueLExpr, fileName="trueGenesPlots.pdf")
#' plotGenesMat (mets=falseLMet, expres=falseLExpr, fileName="falseGenesPlots.pdf")
#'
plotGenesMat <- function(mets, expres, fileName, text4Title=NULL, 
                         x1=1/3, x2=2/3, 
                         y1=NULL, y2=NULL, 
                         percY1=1/3, percY2=2/3,  
                         plotGrid=TRUE, logicSc = NULL){
  if (!is.null(fileName))
    pdf(fileName)
  if (!is.null(text4Title)){
    text4Title <- paste(rownames(expres),text4Title, sep=", ")
  }else{
    if (is.null(logicSc)){
      text4Title<- rownames(expres)
    }else{
      text4Title<- paste(rownames(expres), "\n L-shaped = ", logicSc, sep = " ") #text4Title<- rownames(expres)
    }
  }
#  opt<-par(mfrow=c(2,2))
  for (gene in 1:nrow(expres)){
    xVec<- as.numeric(mets[gene,])
    yVec<- as.numeric(expres[gene,])
    plotGenSel(xMet=xVec, yExp=yVec, titleText=text4Title[gene], 
               x1=x1, x2=x2, percY1=percY1, percY2=percY2, plotGrid=plotGrid) #x1=1/3, x2=2/3, plotGrid=plotGrid)
  }
#  par(opt)
  if (!is.null(fileName))
    dev.off()  
}

#' plotGeneByName
#' \code{plotGeneByName} plot points on a scatterplot with a 3x3 grid overimposed. 
#' The name of a the gene is provided jointly with the matrix and used to select the row to be plotted
#' @examples
#' plotGeneByName (gene="HOOK1", mets=falseLMet, expres=falseLExpr, fileName=NULL)
#'
plotGeneByName <- function(geneName, mets, expres, fileName, text4Title=NULL, 
                         plotGrid=TRUE, figs=c(2,2)){
  if (!is.null(fileName))
    pdf(fileName)
  if (!is.null(text4Title)){
    text4Title <- paste(geneName, text4Title, sep=", ")
  }else{
    text4Title<- geneName
  }
  if (geneName %in% rownames(expres)){
    genePos <- which(rownames(expres)==geneName)
  }else{
    genePos <- NULL
  }
  if(!(is.null(genePos))){
    xVec<- as.numeric(mets[genePos,])
    yVec<- as.numeric(expres[genePos,])
    plotGenSel(xMet=xVec, yExp=yVec, titleText=text4Title, 
               x1=1/3, x2=2/3, plotGrid=plotGrid)
  }
  if (!is.null(fileName))
    dev.off()  
}



