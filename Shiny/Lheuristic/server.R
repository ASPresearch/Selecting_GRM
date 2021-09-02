#---------------------------------------------------------------
# Versio 1.1 Mars-2018
#---------------------------------------------------------------

library(shiny)
# library(lpattern)

# We will use the set of functions in gridFunctions.R
source("gridFunctions.R")
# source("https://raw.githubusercontent.com/alexsanchezpla/scripts/master/Selecting_Genes_Regulated_by_Methylation/Rcode/gridFunctions.R")

#Set working directory
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")

# Set memory size limit to 30 Mb. Intended to allow uploading big files

options(shiny.maxRequestSize=30*1024^2)

#server
shinyServer(function(input, output, clientData, session) {
  
  Ngenes <- reactive({
    if(input$allgenes){
      Ngenes <- NULL
    }else{
    Ngenes <- input$Ngenes
    }
    })
  
  #Load user methylation data or demo methylation data
  metFile <- reactive({
    
    metinFile <- input$metfile
    
    if (!is.null(metinFile)){
      metFile <- read.csv2(metinFile$datapath, header = TRUE,
                       sep = input$sep1, dec = input$dec1, quote = input$quote1, 
                       row.names = 1)
    }else{
      dataset <- input$dataset





      if(dataset=="GEOAll"){
        metFile <- read.table(file.path(dataDir, "GEOMethData.csv"), header = TRUE,
                          sep = ";", dec = ".", quote = '"', 
                          row.names = 1)
      }else{
        metFile <- read.table(file.path(dataDir, "GEOTrueLMetilacion.csv"), header = TRUE,
                          sep = ";", dec = ".", quote = '"', 
                          row.names = 1)
      }
     }





    return(metFile)
    })
  
  metData <- reactive({
    
    Ngenes <- Ngenes()
    metFile <- metFile()
    
    if(!is.null(Ngenes)){
      Ngenes <- ifelse(Ngenes < nrow(metFile), Ngenes, nrow(metFile))
      metData <- metFile[1:Ngenes,]
    }else{
      metData <- metFile
    }
    
    return(metData)
    })
  
  #Load user expression data or demo expression data
  exprFile <- reactive({
    
    exprinFile <- input$exprfile
    Ngenes <- Ngenes()
    
    if(!is.null(exprinFile)){
      exprFile <- read.csv2(exprinFile$datapath, header = TRUE,
                     sep = input$sep2, dec = input$dec2, quote = input$quote2, 
                     row.names = 1)
    }else{
      
      dataset <- input$dataset









      if(dataset=="GEOAll"){
        exprFile <- read.csv2(file.path(dataDir, "GEOExpData.csv"), header = TRUE,
                         sep = ";", dec = ".", quote = '"', 
                         row.names = 1)
      }else{
        exprFile <- read.csv2(file.path(dataDir, "GEOTrueLExpression.csv"), header = TRUE,
                         sep = ";", dec = ".", quote = '"', 
                         row.names = 1)
      }







    }
    return(exprFile)
    })
  
  exprData <- reactive({
    
    Ngenes <- Ngenes()
    exprFile <- exprFile()
    
    if(!is.null(Ngenes)){
      Ngenes <- ifelse(Ngenes < nrow(exprFile), Ngenes, nrow(exprFile))
      exprData <- exprFile[1:Ngenes,]
    }else{
      exprData <- exprFile
    }
    
    return(exprData)
  })
  
  percentMat <- eventReactive(input$submit, {
    
    n11input <- input$n11; n12input <- input$n12; n13input <- input$n13
    n21input <- input$n21; n22input <- input$n22; n23input <- input$n23
    n31input <- input$n31; n32input <- input$n32; n33input <- input$n33
    
    aReqPercentsMat <- matrix(c(n11input,n12input,n13input,
                                n21input,n22input,n23input,
                                n31input,n32input,n33input), 
                              nrow=3, byrow=TRUE)
    N <- dim(metData())[2]
    toReqMat(N, aReqPercentsMat)
  })
  
  output$matrix <- renderTable({
    m <- percentMat()
    matrix(as.integer(m), nrow=3)
    }, colnames = FALSE)

  
  #Build a table with columns, logicSc and numericSc, and gene symbol as row name. 
  metTable <- eventReactive(input$submit, {


























      n11input <- input$n11; n12input <- input$n12; n13input <- input$n13
      n21input <- input$n21; n22input <- input$n22; n23input <- input$n23
      n31input <- input$n31; n32input <- input$n32; n33input <- input$n33
      
      aReqPercentsMat <- matrix(c(n11input, n12input, n13input,
                                  n21input, n22input, n23input,
                                  n31input, n32input, n33input), 
                                  nrow=3, byrow=TRUE)
      
      w11input <- input$w11; w12input <- input$w12; w13input <- input$w13
      w21input <- input$w21; w22input <- input$w22; w23input <- input$w23
      w31input <- input$w31; w32input <- input$w32; w33input <- input$w33
      
      neWeightM <- matrix(c(w11input, w12input, w13input,
                            w21input, w22input, w23input,
                            w31input, w32input, w33input), 
                          nrow=3, byrow=TRUE)
      
      wf11input <- input$wf11; wf12input <- input$wf12; wf13input <- input$wf13
      wf21input <- input$wf21; wf22input <- input$wf22; wf23input <- input$wf23
      wf31input <- input$wf31; wf32input <- input$wf32; wf33input <- input$wf33
      
      neFWeightM <- matrix(c(wf11input, wf12input, wf13input,
                             wf21input, wf22input, wf23input,
                             wf31input, wf32input, wf33input), 
                             nrow=3, byrow=TRUE)
      
      trueScores <- scoreGenesMat (mets=metData(), expres = exprData(),
                                   x1=input$xaxis[1], x2=input$xaxis[2],
                                   y1=NULL, y2=NULL, percY1 = input$yaxis[1], percY2 = input$yaxis[2],
                                   aReqPercentsMat=aReqPercentsMat,
                                   aWeightMifL=neWeightM,
                                   aWeightMifNonL=neFWeightM)


    trueScores
    })
  
  #Subset the true L-shaped genes from "metTable"
  seltable <- reactive({
    allGenes <- metTable()
    selGenes <- subset(allGenes, allGenes$logicSc == TRUE)
    selGenes
  })
  
  #Subset the false L-shaped genes from "metTable"
  unseltable <- reactive({
    allGenes <- metTable()
    unselGenes <- subset(allGenes, allGenes$logicSc == FALSE)
    unselGenes
  })
  
  #Count and output the number of genes classified as true L-shaped genes
  output$nselect1 <- output$nselect2 <- renderText({
    
    selGenes <- rownames(seltable())
    
    paste("L-shaped Genes: ", length(selGenes))
  })
  
  #Count and output the number of genes classified as false L-shaped genes
  output$nunselect1 <- output$nunselect2 <- renderText({
    
    unselGenes <- rownames(unseltable())
    
    paste("Non L-shaped Genes: ", length(unselGenes))
  })
  
  #Write table to a .csv file, with all results, to be downloaded by the user
  output$downloadAll <- downloadHandler(
    filename = paste("AllGenes_", Sys.Date(), ".csv", sep=''),
    content = function(file){write.csv(metTable(), file)}
    )
  
  #Write table to a .csv file, with true L-shaped genes results, to be downloaded 
  #by the user
  output$downloadSelect <- downloadHandler(
    filename = function() {paste("LshapedGenes_", Sys.Date(), ".csv", sep="")},
    content = function(file) {write.csv(seltable(), file)}
    )
  
  #Write table to a .csv file, with false L-shaped genes results, to be downloaded 
  #by the user
  output$downloadUnsel <- downloadHandler(
    filename = paste("NonLshapedGenes_", Sys.Date(), ".csv", sep=''),
    content = function(file){write.csv(unseltable(), file)}
    )
  
  #Display the n first rows of the full table of genes. User can choose the n
  #of rows to be displayed
  output$alltable <- renderTable({
    
    n <- as.numeric(input$Nall)
    alltable <- metTable()
    if(n > length(alltable[,1])){
      n <- length(alltable[,1])
    }
    alltable[1:n, ]
    
    }, colnames = TRUE, rownames = TRUE)
  
  #Plot the scatterplots of all genes. These graphics are saved as a PDF file that can be 
  #downloaded by the user
  allplots <- eventReactive(input$submit,{
    
    metData <- metData()
    alltable <- metTable()
    allGenes <- rownames(alltable)
    
    mets <- subset(metData, rownames(metData) %in% allGenes)
    expres <- subset(exprData(), rownames(exprData()) %in% allGenes)
    n <- ifelse(length(mets[,1])<4, length(mets[,1]), 4)
    opt<- par(mfrow=c(2,2))

    plotGenesMat(mets=mets[1:n,], 
                 expres=expres[1:n,], 
                 x1=input$xaxis[1], x2=input$xaxis[2], logicSc = alltable$logicSc,
                 percY1=input$yaxis[1], percY2=input$yaxis[2], fileName = NULL)
    par(opt)

    pdf("results/temp/allplots.pdf")
    opt<- par(mfrow=c(2,2))
    plotGenesMat(mets=subset(metData, rownames(metData) %in% allGenes), 
                 expres=subset(exprData(), rownames(exprData()) %in% allGenes), 
                 x1=input$xaxis[1], x2=input$xaxis[2], 
                 y1=NULL, y2=NULL, percY1=input$yaxis[1], percY2=input$yaxis[2], 
                 logicSc = alltable$logicSc,
                 fileName = NULL)
    par(opt)
    dev.off()
  })
  
  output$allplots <- renderPlot({allplots()})
  
  output$pdflinkAll <- downloadHandler(
    
    filename = paste("Allplots_", Sys.Date(),".pdf", sep=""),
    content <- function(file) {
      file.copy("results/temp/allplots.pdf", file)
    })
  
  #Display the n first rows of the true L-shaped genes' table. User can choose the n
  #of rows to be displayed
  output$selectable <- renderTable({
    
    n <- as.numeric(input$Nsel)
    seltable <- seltable()
    if(n > length(seltable[,1])){
      n <- length(seltable[,1])
    }
    seltable[1:n, ]
    
    }, colnames = TRUE, rownames = TRUE)
  
  #Plot the scatterplots of the true L-shaped genes. These graphics are saved as a PDF file that 
  #can be downloaded by the user
  selectplots <- eventReactive(input$submit,{
    
    metData <- metData()
    seltable <- seltable()
    selGenes <- row.names(seltable())
    
    mets <- subset(metData, rownames(metData) %in% selGenes)
    expres <- subset(exprData(), rownames(exprData()) %in% selGenes)
    n <- ifelse(length(mets[,1])<4, length(mets[,1]), 4)
    opt<- par(mfrow=c(2,2))

    plotGenesMat(mets=mets[1:n,], 
                 expres=expres[1:n,], 
                 x1=input$xaxis[1], x2=input$xaxis[2], logicSc = seltable$logicSc,
                 y1=NULL, y2=NULL, percY1=input$yaxis[1], percY2=input$yaxis[2], fileName = NULL)
    par(opt)

    pdf("results/temp/selectplots.pdf")
    opt<- par(mfrow=c(2,2))
    plotGenesMat(mets=subset(metData, rownames(metData) %in% selGenes), 
                 expres=subset(exprData(), rownames(exprData()) %in% selGenes), 
                 x1=input$xaxis[1], x2=input$xaxis[2], logicSc = seltable$logicSc,
                 y1=NULL, y2=NULL, percY1=input$yaxis[1], percY2=input$yaxis[2], fileName = NULL)
    par(opt)
    dev.off()
  })
  
  output$selectplots <- renderPlot({selectplots()})
  
  output$pdflinkSelect <- downloadHandler(
    
    filename = paste("Lshapedplots_", Sys.Date(),".pdf", sep=""),
    content <- function(file) {
      file.copy("results/temp/selectplots.pdf", file)
    })
  
  #Display the n first rows of the false L-shaped genes' table. User can choose the n
  #of rows to be displayed
  output$unseltable <- renderTable({
    
  n <- as.numeric(input$Nunsel)
  unseltable <- unseltable()
  if(n > length(unseltable[,1])){
    n <- length(unseltable[,1])
    }
  unseltable[1:n, ]
    
    }, colnames = TRUE, rownames = TRUE)
  
  #Plot the scatterplots of the false L-shaped genes. These graphics are saved as a PDF file that 
  #can be downloaded by the user
  
  unselplots <- eventReactive(input$submit,{
    
    metData <- metData()
    unseltable <- unseltable()
    unselGenes <- row.names(unseltable())
    
    mets <- subset(metData, rownames(metData) %in% unselGenes)
    expres <- subset(exprData(), rownames(exprData()) %in% unselGenes)
    n <- ifelse(length(mets[,1])<4, length(mets[,1]), 4)
    opt<- par(mfrow=c(2,2))

    plotGenesMat(mets=mets[1:n,], 
                 expres=expres[1:n,],
                 x1=input$xaxis[1], x2=input$xaxis[2], logicSc = unseltable$logicSc,
                 y1=NULL, y2=NULL, percY1=input$yaxis[1], percY2=input$yaxis[2], fileName = NULL)
    par(opt)

    pdf("results/temp/unselplots.pdf")
    opt<- par(mfrow=c(2,2))
    plotGenesMat(mets=subset(metData, rownames(metData) %in% unselGenes), 
                 expres=subset(exprData(), rownames(exprData()) %in% unselGenes), 
                 x1=input$xaxis[1], x2=input$xaxis[2], logicSc = unseltable$logicSc,
                 y1=NULL, y2=NULL, percY1=input$yaxis[1], percY2=input$yaxis[2], fileName = NULL)
    par(opt)
    dev.off()
  })
  
  output$unselplots <- renderPlot({unselplots()})
  
  output$pdflinkUnsel <- downloadHandler(
    
    filename = paste("NonLshapedplots_", Sys.Date(),".pdf", sep=""),
    content <- function(file) {
      file.copy("results/temp/unselplots.pdf", file)
    })
  
  #Select all genes or a set of genes
  observe({
    allgenes <- input$allgenes
    if(allgenes){
      updateNumericInput(session, "Ngenes", value = "")
    }else{
      updateNumericInput(session, "Ngenes", value = 200)
    }
  })
  
  #Changing the value of input$n31 from the server, so that the total sum of cells does not 
  #exceed 100%
  observe({
    n11 <- input$n11; n12 <- input$n12; n13 <- input$n13
    n21 <- input$n21; n22 <- input$n22; n23 <- input$n23
    n31 <- input$n31; n32 <- input$n32; n33 <- input$n33
    
#    updateNumericInput(session, "n11", max = 100-(n12+n13+n21+n22+n23+n31+n32+n33))
#    updateNumericInput(session, "n12", max = 100-(n11+n13+n21+n22+n23+n31+n32+n33))
#    updateNumericInput(session, "n13", max = 100-(n11+n12+n21+n22+n23+n31+n32+n33))
#    updateNumericInput(session, "n21", max = 100-(n11+n12+n13+n22+n23+n31+n32+n33))
#    updateNumericInput(session, "n22", max = 100-(n11+n12+n13+n21+n23+n31+n32+n33))
#    updateNumericInput(session, "n23", max = 100-(n11+n12+n13+n21+n22+n31+n32+n33))
#    updateNumericInput(session, "n31", max = 100-(n11+n12+n13+n21+n22+n23+n32+n33))
#    updateNumericInput(session, "n32", max = 100-(n11+n12+n13+n21+n22+n23+n31+n33))
#    updateNumericInput(session, "n33", max = 100-(n11+n12+n13+n21+n22+n23+n31+n32))
  })
  
  #Reset parameters to default
  observeEvent(input$reset, {
    updateNumericInput(session, "Ngenes", value = 200)
    updateCheckboxInput(session, "allgenes", value=FALSE)
    updateSliderInput(session, "xaxis", value = c(1/3, 2/3))
    updateSliderInput(session, "yaxis", value = c(1/3, 2/3))
    updateNumericInput(session, "n11", value = 10)
    updateNumericInput(session, "n12", value = 20)
    updateNumericInput(session, "n13", value = 1)
    updateNumericInput(session, "n21", value = 5)
    updateNumericInput(session, "n22", value = 40)
    updateNumericInput(session, "n23", value = 20)
    updateNumericInput(session, "n31", value = 0)
    updateNumericInput(session, "n32", value = 5)
    updateNumericInput(session, "n33", value = 10)
    updateNumericInput(session, "w11", value = 2)
    updateNumericInput(session, "w12", value = -2)
    updateNumericInput(session, "w13", value = -25)
    updateNumericInput(session, "w21", value = 1)
    updateNumericInput(session, "w22", value = 0)
    updateNumericInput(session, "w23", value = -2)
    updateNumericInput(session, "w31", value = 1)
    updateNumericInput(session, "w32", value = 1)
    updateNumericInput(session, "w33", value = 2)
    updateNumericInput(session, "wf11", value = 0)
    updateNumericInput(session, "wf12", value = -2)
    updateNumericInput(session, "wf13", value = -25)
    updateNumericInput(session, "wf21", value = 1)
    updateNumericInput(session, "wf22", value = 0)
    updateNumericInput(session, "wf23", value = -2)
    updateNumericInput(session, "wf31", value = 1)
    updateNumericInput(session, "wf32", value = 1)
    updateNumericInput(session, "wf33", value = 2)
  })
})