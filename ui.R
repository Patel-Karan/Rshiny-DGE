library(shiny)
library(DT)
library(shinythemes)
library(shinycssloaders)
library(plotly)

fluidPage(theme = shinytheme("cosmo"),includeCSS("www/styles.css"),
    #titlePanel("Differential Expression Analysis"),
    navbarPage(title="Logo",tabPanel("Analysis",
                             sidebarLayout(
                               sidebarPanel(
                                 fileInput("file1", "Upload file in CSV", accept = ".csv", 
                                           placeholder = "Enter Count table in csv"),
                                 radioButtons("tool", "Tool", c("DESeq2", "edgeR")),
                                 selectInput("control", "Control", character(0), multiple = TRUE),
                                 selectInput("case", "Case", character(0), multiple = TRUE),
                                 actionButton("submit","Run"),
                               ),
                               mainPanel(
                                 fluidRow(
                                   column(4,
                                          selectInput("fformat","Select File Format",
                                                      choices = c("csv","xlsx","txt"), selected = "csv")
                                   ),
                                   column(4,textInput("comp","Enter Comparision Name","Control_vs_Case")),
                                   column(4, tags$br(),downloadButton("downloadData", "Download", 
                                                            icon = icon("download")))
                                 ),
                                 withSpinner(
                                   DTOutput(outputId = 'table')
                                 ))
                             )
                             ),
               tabPanel("Hierarchical Clustering & Correlation Analysis",
                        sidebarLayout(sidebarPanel(
                          selectInput("dmethod","Select distance method", 
                                      choices = c("euclidean","maximum","manhattan",
                                                  "canberra","binary","minkowski"),
                                      selected="euclidean"),
                          selectInput("hmethod", "Select hclust methods", 
                                      choices = c("ward.D", "ward.D2", "single", "complete", 
                                                  "average", "mcquitty", "median", "centroid"),
                                      selected = "complete"),
                          selectInput("hclust_Samples","Select Samples for Hclust"
                                      ,character(0), multiple = TRUE),
                          actionButton("run_hclust","Run")
                        ),
                        mainPanel(
                          navbarPage(NULL,
                                     tabPanel("Distance Matrix",
                                              withSpinner(
                                                DTOutput(outputId = 'table2')
                                              )
                                              ),
                                     tabPanel("Dendogram",
                                              withSpinner(
                                                plotlyOutput("dendoplot")
                                              )),
                                     tabPanel("PCA",
                                              withSpinner(
                                       plotlyOutput("pca")
                                     )),
                                     tabPanel("BoxPlot",
                                              withSpinner(
                                                plotlyOutput("boxplot")
                                              )),
                                     tabPanel("Correlation"))
                        ))),
               tabPanel("Volcano Plot",
                        fluidRow(column(3,
                                        tags$br(),actionButton("plotvol", "Generate Plot", width = "80%"),
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
               tabPanel("Heatmaps"),
               tabPanel("About")
               )
    )
