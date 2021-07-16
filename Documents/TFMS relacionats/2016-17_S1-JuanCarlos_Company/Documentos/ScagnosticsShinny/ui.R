library(shiny)

shinyUI(navbarPage("L-Shape",
                   tabPanel("Scagnostics",
                            sidebarLayout(
                              sidebarPanel(
                                tabsetPanel(
                                  tabPanel("Upload data",
                                           #br(),h4("Add GEO id"),
                                           #textInput("text", label = ""),
                                           hr(), h4("Choose input files"),
                                           fileInput('metfile', 'Upload your methylation array',
                                                     accept = c('text/csv',
                                                                '.csv')),
                                           fileInput('exprfile', 'Upload your expression microarray or RNAseq',
                                                     accept = c('text/csv','.csv'), ),
                                           p(strong('Set format parameters of your',
                                                    'methylation data file')),
                                           fluidRow(
                                             column(4,radioButtons('sep1', 'Separator',
                                                                   c(Tab='\t',
                                                                     Comma=',',
                                                                     Semicolon=';'),
                                                                   ';')),
                                             column(4,radioButtons('dec1', 'Decimal',
                                                                   c(Point='.',
                                                                     Comma=','),
                                                                   ',')),
                                             column(4,radioButtons('quote1', 'Quote',
                                                                   c(None='',
                                                                     'Double'='"',
                                                                     'Single'="'"),
                                                                   '"'))),
                                           p(strong('Set format parameters of your',
                                                    'expression data file')),
                                           fluidRow(
                                             column(4,radioButtons('sep2', 'Separator',
                                                                   c(Tab='\t',
                                                                     Comma=',',
                                                                     Semicolon=';'),
                                                                   ';')),
                                             column(4,radioButtons('dec2', 'Decimal',
                                                                   c(Point='.',
                                                                     Comma=','),
                                                                   ',')),
                                             column(4,radioButtons('quote2', 'Quote',
                                                                   c(None='',
                                                                     'Double'='"',
                                                                     'Single'="'"),
                                                                   '"')
                                             )
                                           )
                                  ),
                                  
                                  tabPanel("L-shape Settings",
                                           h3("Parameters", align = "left"),
                                           fluidRow(column(8,
                                                           actionButton("action", label = "Compute"))),br(),
                                           sliderInput("Convex", "convex",
                                                       min = 0, max = 1, value = c(0, 0.1)),
                                           sliderInput("Skinny", "skinny",
                                                       min = 0, max = 1, value = c(0.2, 0.8)),
                                           sliderInput("Stringy", "stringy",
                                                       min = 0, max = 1, value = c(0.6, 1))
#                                            sliderInput("Outlying", "outlying",
#                                                        min = 0, max = 2, value = c(0, 1)),
#                                            sliderInput("Skewed", "skewed",
#                                                        min = 0, max = 1, value = c(0,1)),
#                                            sliderInput("Clumpy", "clumpy",
#                                                        min = 0, max = 1, value = c(0, 1)),
#                                            sliderInput("Sparse", "sparse",
#                                                        min = 0, max = 1, value = c(0, 1)),
#                                            sliderInput("Striated", "striated",
#                                                        min = 0, max = 1, value = c(0, 1)),
#                                            sliderInput("Sonotonic", "monotonic",
#                                                        min = 0, max = 1, value = c(0, 1))
                                  ))),
                              
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("L-shaped Genes", br(),
                                           column(4, downloadButton("downloadSelect", "Download table")),
                                           fluidRow(tableOutput("selectable")) ,
                                           hr(),
                                           p(strong('Graphic')),
                                           fluidRow(column(4,textInput("geneName", label = "")),
                                           column(8, actionButton("showChart", label = "Show Graphic"))),
                                           imageOutput("someplots")
                                  
                                  )
                                  
                                )
                              )
                            )
                   ),
                   tabPanel("Help")
)
)
