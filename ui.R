library(shiny)
library(DT)
library(shinythemes)

# Define UI for application that draws a histogram
fluidPage(theme = shinytheme("united"),

    # Application title
    titlePanel("Differential Expression Analysis"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          fileInput("file1", "Upload files in CSV"),
          radioButtons("tool", "Tool", c("DESeq2", "edgeR")),
          selectInput("control", "Control", character(0), multiple = TRUE),
          selectInput("case", "Case", character(0), multiple = TRUE),
          textInput("comp","Comparision","Control_vs_Case"),
          actionButton("submit","Run")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          DTOutput(outputId = 'table')
        )
    )
)
