library(shiny)
library(DT)
library(shinythemes)
library(shinycssloaders)
library(plotly)

fluidPage(theme = shinytheme("cosmo"),includeCSS("www/styles.css"),
    #titlePanel("Differential Expression Analysis"),
    navbarPage(title="Logo",
               tabPanel("Hierarchical Clustering & Correlation Analysis",
                        sidebarLayout(sidebarPanel(
                          fileInput("file1", "Upload file in CSV", accept = ".csv", 
                                    placeholder = "Enter Count table in csv"),
                          selectInput("dmethod","Dissimilarity measures", 
                                      choices = c("euclidean","maximum","manhattan",
                                                  "canberra","binary","minkowski"),
                                      selected="euclidean"),
                          selectInput("hmethod", "Linkage methods", 
                                      choices = c("ward.D", "ward.D2", "single", "complete", 
                                                  "average", "mcquitty", "median", "centroid"),
                                      selected = "complete"),
                          fluidRow(column(9,actionButton("run_hclust","Run")),
                                   column(3,p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>"))))
                        ),
                        mainPanel(
                          selectInput("res", "Results", c("Distance Matrix","Dendogram", "PCA", "Boxplot", "Correlation Plot"), 
                                      selected = "Distance Matrix"),
                          conditionalPanel(condition = "input.res == 'Distance Matrix'",
                                           withSpinner(DTOutput(outputId = 'table2'))),
                          conditionalPanel(condition = "input.res == 'Dendogram'",
                                           withSpinner(plotlyOutput("dendoplot"))),
                          conditionalPanel(condition = "input.res == 'PCA'",
                                           withSpinner(plotlyOutput("pca"))),
                          conditionalPanel(condition = "input.res == 'Boxplot'",
                                           withSpinner(plotlyOutput("boxplot"))),
                          conditionalPanel(condition = "input.res == 'Correlation Plot'",
                                           withSpinner(plotlyOutput("cor_plot", height = "100%")))),
                        )),
               tabPanel("Differenetial Gene Expression Analysis",
                             sidebarLayout(
                               sidebarPanel(
                                 p(HTML("<b>DGE Parameters</b>")),
                                 radioButtons("tool", "Tool", c("DESeq2", "edgeR")),
                                 selectInput("control", "Control", character(0), multiple = TRUE),
                                 selectInput("case", "Case", character(0), multiple = TRUE),
                                 textInput("comp","Enter Comparision Name","Control_vs_Case"),
                                 p(HTML("<b>Data filteration Parameter</b>")),
                                 numericInput("pvalue","P-Value", value = 0.05),
                                 numericInput("log2fc","log2FC", value = 1),
                                 numericInput("fdr","FDR", value = NULL),
                                 actionButton("submit","Run")
                               ),
                               mainPanel(
                                 tabsetPanel(tabPanel("Result table", tags$br(),
                                                      div(align="right", 
                                                          downloadButton("downloadData", 
                                                                         "Download", 
                                                                         icon = icon("download"))),
                                                      tags$br(),
                                                      withSpinner(DTOutput('table'))),
                                             tabPanel("Filtered Table", tags$br(),
                                                      div(align="right", 
                                                          downloadButton("downloadData2", 
                                                                         "Download", 
                                                                         icon = icon("download"))),
                                                      tags$br(),
                                                      withSpinner(DTOutput('table3'))),
                                             tabPanel("Result Summary", tags$br(),
                                                      div(align="right", 
                                                          downloadButton("downloadData3", 
                                                                         "Download", 
                                                                         icon = icon("download"))),
                                                      tags$br(),
                                                      withSpinner(plotlyOutput('sum', height = "100%"))
                                                      )
                                             )
                             ))),
               tabPanel("Volcano Plot",
                        fluidRow(column(3,
                                        tags$br(),actionButton("plotvol", "Generate Plot", width = "80%")
                        ),
                        column(3,selectInput("vformat","Select File Format",
                                             choices = c("png","jpeg","pdf"), selected = "png")
                        ),
                        column(3,numericInput("vol_width", "Width in px",13),
                               numericInput("vol_height", "Height in px",7)),
                        column(3,tags$br(), downloadButton("downloadPlotvol", "Download",
                                                icon = icon("download")))
                        ),tags$br(),
                        withSpinner(
                          plotlyOutput("volcano", height = "100%", width = "80%")
                        )),
               tabPanel("Heatmaps",
                        sidebarLayout(sidebarPanel(
                          radioButtons("strategy", "Strategy", c("Top n", "Selected Genes")),
                          conditionalPanel(condition = "input.strategy == 'Top n'",
                                           numericInput("topn", "N number of genes", value = 25)),
                          conditionalPanel(condition = "input.strategy == 'Selected Genes'",
                                           p(HTML("<b style=\"color:black\">Enter list of genes</b>")),
                                           tags$style(type = "text/css", "textarea {width:100%}"),
                                           tags$textarea(
                                             id = "genes_list", placeholder = "Just paste a list of genes. One genes in one line.",
                                             rows = 8, ""
                                           )),
                          actionButton("submit","Run") 
                                           #textAreaInput("input_text", "Enter list of genes"))
                        ),
                        mainPanel(
                          withSpinner(plotlyOutput('heatmaps', height = "1000px"))
                        )
                        )),
               tabPanel("About")
               )
    )
 
