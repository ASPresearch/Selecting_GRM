library(shiny)
source("load_and_create_data_Latest.R")
#Set working directory
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")

#server
shinyServer(function(input, output, clientData, session) {

  #Load user methylation data or demo methylation data
  metFile <- reactive({
    
    metinFile <- input$metfile
    
    if (!is.null(metinFile)){
      metFile <- read.csv2(metinFile$datapath, header = TRUE,
                           sep = input$sep1, dec = input$dec1, quote = input$quote1, 
                           row.names = 1, stringsAsFactors = F)
    }
    return(metFile)
  })
  
  #Load user expression data or demo expression data
  exprFile <- reactive({
    
    exprinFile <- input$exprfile
    #Ngenes <- Ngenes()
    
    if(!is.null(exprinFile)){
      exprFile <- read.csv2(exprinFile$datapath, header = TRUE,
                            sep = input$sep2, dec = input$dec2, quote = input$quote2, 
                            row.names = 1, stringsAsFactors = F)
    }
    return(exprFile)
  })  

  ### --------------------------------- ###
# My own functions
# You can access the value of the widget with input$text, e.g.
  output$value <- renderPrint({ 
    value.to.show = paste("Number of Exp Samples set:",
                          ncol(exprFile()), 
                          "Number of Meth Samples set:",
                          ncol(metFile()),"",sep=" ")
  })
  

# # Lectura de los archivos  
#   geo.read <- reactive({
#       geo.val <- getGEO.anno(input$text)
#       # Return Matrix GEO annotated
#       return(geo.val)  
#   })  
  
# output$selectable <- renderDataTable({
#   show.table <- geo.read()
# })

ExpReadTable <- reactive({
  # Get same name input
  exp = exprFile()
  meth = metFile()
  expr.gene = rownames(exp)
  meth.gene = rownames(meth)
  f.exp.gene= exp[expr.gene %in% meth.gene, ]
  f.meth.gene= exp[meth.gene %in% expr.gene, ]
  o.exp = f.exp.gene[order(rownames(f.exp.gene)), ]
  o.meth = f.meth.gene[order(rownames(f.meth.gene)), ]
  return(o.exp)
})  
  
MethReadTable <- reactive({
  # Get same name input
  exp = exprFile()
  meth = metFile()
  expr.gene = rownames(exp)
  meth.gene = rownames(meth)
  f.exp.gene= exp[expr.gene %in% meth.gene, ]
  f.meth.gene= meth[meth.gene %in% expr.gene, ]
  o.exp = f.exp.gene[order(rownames(f.exp.gene)), ]
  o.meth = f.meth.gene[order(rownames(f.meth.gene)), ]
  return(o.meth)
})  

# shape de los valores
scagnostics <- eventReactive(input$action, {

    # Get Shape Files
  # Add filters to reduce the size of the set. 
  shape.val = lapply(rownames(ExpReadTable()), function(x) shape.param(exp = ExpReadTable()[x, ], meth=MethReadTable()[x, ] ) )
  
#   # Obtain Values
#   # Clean shape table
   shape.table <-  table(NA,NA,NA)
   for (i in shape.val) {  shape.table = rbind(shape.table, i ) }
   is.null.vector = as.vector(unlist(lapply(1:length(shape.val), function(x) is.null(shape.val[[x]])  )))
   shape.tabl <-  as.data.frame(shape.table)
   rownames(shape.tabl) = rownames(ExpReadTable())
   print(shape.tabl)
   # Filter by setting selecction
   # L-shape Default
   param.Convex = filt.fun(min = input$Convex[1], max = input$Convex[2], param = "Convex",set=shape.tabl)
   param.Skinny = filt.fun(min = input$Skinny[1], max = input$Skinny[2], param = "Skinny",set=shape.tabl)
   param.Stringy = filt.fun(min = input$Stringy[1], max = input$Stringy[2], param = "Stringy",set=shape.tabl)
   
   # Other Elements
#    param.Outlying = filt.fun(min = input$Outlying[1], max = input$Outlying[2], param = "Outlying",set=shape.tabl)
#    param.Skewed = filt.fun(min = input$Skewed[1], max = input$Skewed[2], param = "Skewed",set=shape.tabl)
#    param.Clumpy = filt.fun(min = input$Clumpy[1], max = input$Clumpy[2], param = "Clumpy",set=shape.tabl)
#    param.Sparse = filt.fun(min = input$Sparse[1], max = input$Sparse[2], param = "Sparse",set=shape.tabl)
#    param.Striated = filt.fun(min = input$Striated[1], max = input$Striated[2], param = "Striated",set=shape.tabl)
#    param.Monotonic = filt.fun(min = input$Monotonic[1], max = input$Monotonic[2], param = "Monotonic",set=shape.tabl)
#    
   # Selection
   # all.param = which(1:ncol(shape.tabl) %in% param.Outlying & 1:ncol(shape.tabl) %in% param.Skewed & 1:ncol(shape.tabl) %in% param.Clumpy & 1:ncol(shape.tabl) %in% param.Sparse & 1:ncol(shape.tabl) %in% param.Striated & 1:ncol(shape.tabl) %in% param.Convex & 1:ncol(shape.tabl) %in% param.Skinny & 1:ncol(shape.tabl) %in% param.Stringy & 1:ncol(shape.tabl) %in% param.Monotonic )
   
   all.param = which(1:ncol(shape.tabl) %in% param.Convex & 1:ncol(shape.tabl) %in% param.Skinny & 1:ncol(shape.tabl) %in% param.Stringy)
    shape.tabl = shape.tabl[all.param, ]
  # Return table
  return(shape.tabl)
  
})

# Show L-shape defined user
output$selectable <- renderTable({
      seltable <- scagnostics()
    }, rownames = T)

# Some plot example
someplots <- eventReactive(input$showChart, {
  gene.name = rownames(scagnostics()[input$geneName,  ])
  scatter.expMeth(as.numeric(ExpReadTable()[gene.name,]), as.numeric(MethReadTable()[gene.name,] ), gen.name = gene.name )
})

output$someplots <- renderPlot({
    someplots() 
  })

#Write table to a .csv file, with true L-shaped genes results, to be downloaded 
  #by the user
  output$downloadSelect <- downloadHandler(
    filename = function() {paste("LshapedGenes_", Sys.Date(), ".csv", sep="")},
    content = function(file) {write.csv(scagnostics(), file)}
  )
})