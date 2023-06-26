library(shiny)
library(DT)
library(shinythemes)
library(shinycssloaders)
library(plotly)

# Define UI for application that draws a histogram
fluidPage(theme = shinytheme("united"),
    titlePanel("Differential Expression Analysis"),
    navbarPage(NULL,tabPanel("Analysis",
                             sidebarLayout(
                               sidebarPanel(
                                 fileInput("file1", "Upload file in CSV"),
                                 radioButtons("tool", "Tool", c("DESeq2", "edgeR")),
                                 selectInput("control", "Control", character(0), multiple = TRUE),
                                 selectInput("case", "Case", character(0), multiple = TRUE),
                                 textInput("comp","Comparision","Control_vs_Case"),
                                 actionButton("submit","Run"),
                                 downloadButton("downloadData", "Download")
                               ),
                               mainPanel(
                                 withSpinner(
                                   DTOutput(outputId = 'table')
                                 ))
                             )
                             ),
               tabPanel("Hierarchical Clustering & Correlation Analysis"),
               tabPanel("Volcano Plot",
                        withSpinner(
                          plotlyOutput("volcano")
                        )),
               tabPanel("Heatmaps"),
               tabPanel("About")
               )
    )
