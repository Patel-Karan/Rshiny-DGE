library(shiny)
library(dplyr)
library(DT)

options(shiny.maxRequestSize=30*1024^2)

function(input, output, session) {
  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    tbl <- read.table(inFile$datapath, sep=",", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
    return(tbl)
  })
  
  observeEvent(data(),{
    updateSelectInput(session, "control",choices = names(data()))
    })
  
  observeEvent(data(),{
    updateSelectInput(session, "case",choices = names(data()))
  })
  
  observeEvent(input$submit,{
    output$table <- renderDT({
      d_filter <- datatable(data() %>% select(c(input$control,input$case))
                            ,filter = "top", options = list(scrollX=TRUE))
      return(d_filter)
      })
  })
}
