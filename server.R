library(shiny)
library(dplyr)
library(tibble)
library(DT)
library(plotly)

options(shiny.maxRequestSize=30*1024^2)

function(input, output, session) {
  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    tbl <- read.table(inFile$datapath, sep=",", header = TRUE, 
                      row.names = 1,stringsAsFactors = FALSE)
    return(tbl)
  })
  
  observeEvent(data(),{
    updateSelectInput(session, "control",choices = names(data()))
    })
  
  observeEvent(data(),{
    updateSelectInput(session, "case",choices = names(data())) 
  })
  
  df <- eventReactive(input$submit,{
      d_filter <- data() %>% select(c(input$control,input$case))
      nC <- length(input$control)
      nT <- length(input$case)
      source("Scripts/dexp_deseq2.R", local = TRUE)
      dexp <- d_exp(d_filter,nC,nT)
      dexp <- dexp %>% rownames_to_column(var = "GeneID")
      return(dexp)
      })
  
  output$table <- renderDT({datatable(df() , filter = "top", 
                                              options = list(scrollX=TRUE),
                                              rownames = FALSE)})
  
  output$downloadData <- downloadHandler(filename = function() {
      paste(input$comp, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(df(), file, row.names = FALSE)
    }
  )
  
  output$volcano <- renderPlotly({
    plot_ly(data=df()%>%filter(padj != "NA"), x=~log2FoldChange, y=~(-log2(padj)), mode="markers",type = "scatter") %>% 
      layout(title = "Volcano Plot", xaxis = list(title = 'Fold Change'), 
             yaxis = list(title = '-log(FDR)'))
  }) 
  
  }
