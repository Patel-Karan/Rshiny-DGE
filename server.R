library(shiny)
library(dplyr)
library(tidyr)
library(reshape2)
library(tibble)
library(DT)
library(plotly)
library(ggplot2)
library(scales)
library(ggdendro)
library(ggcorrplot)
library(heatmaply)

options(shiny.maxRequestSize=30*1024^2)
options(java.parameters = "-Xmx30g")


function(input, output, session){

### Data Reading.
  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    tbl <- read.table(inFile$datapath, sep=",", header = TRUE, 
                      row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
    return(tbl)
  })

### Hierarchical Clustering & Correlation Analysis
  df_hclust <- eventReactive(input$run_hclust,{
    source("Scripts/Then.R")
    validate(
      need(try(input$file1), "Please upload the Count data in CSV format in Analysis Section!")
    )
    if (length(names(data()))/2 == 0){
      nC <- length(names(data()))/2
      nT <- length(names(data()))-nC
      source("Scripts/dexp_hclust.R", local = TRUE)
      dexp_hclust <- d_exp(data(),nC,nT,input$hmethod,input$dmethod)
    }else{
      nC <- round(length(names(data()))/2)
      nT <- length(names(data()))-nC
      source("Scripts/dexp_hclust.R", local = TRUE)
      dexp_hclust <- d_exp(data(),nC,nT,input$hmethod,input$dmethod)
    }
    return(dexp_hclust)
  })
  
### Distance Matrix table rendering
  output$table2 <- renderDT({
    datatable( df_hclust()$distance, filter = "top", options = list(scrollX=TRUE, rownames = FALSE))
  })
  
### Dendogram plot rendering.
  output$dendoplot <- renderPlotly({
    dendop <- ggdendrogram(df_hclust()$clustering, rotate = FALSE, leaf_labels = FALSE, size = 2)
    ggplotly(dendop)
  })
  
### PCA Plot
  output$pca <- renderPlotly(({
    source("Scripts/PCA_Plot.R")
    ggplotly(plt_pca(df_hclust()$pca), tooltip = c("text"))
  }))
  
### BoxPlot
  output$boxplot <- renderPlotly({
    source("Scripts/BoxPlot.R")
    ggplotly(plt_box(df_hclust()$result_df)) 
  })
  
### Correlation Plot
  output$cor_plot <- renderPlotly({
    # corp <- ggcorrplot(df_hclust()$corr$correlation,p.mat = df_hclust()$corr$pval, 
    #                    sig.level = 0.01, insig = "blank",
    #                    colors = c("red", "white", "blue"),)
    # ggplotly(corp)
    heatmaply_cor(df_hclust()$corr$correlation, colors = c("red", "white", "blue"), show_dendrogram = c(FALSE,FALSE))
  })
  
### Differential Gene Expression Analysis
  observeEvent(data(),{
    updateSelectInput(session, "control",choices = names(data()))
  })
  observeEvent(data(),{
    updateSelectInput(session, "case",choices = names(data())) 
  })
  df <- eventReactive(input$submit,{
    source("Scripts/Then.R")
    validate(
      need(try(input$file1), "Please upload the Count data in CSV format!") %then%
        need(try(input$control), "Please select the Control samples!") %then%
        need(try(input$case), "Please select the Treated/Case Samples!")
    )
    source("Scripts/dexp_deseq2.R", local = TRUE)
    dexp <- diff_exp_deseq2(data(),input$control,input$case)
    return(dexp)
  })
  
### DGE Analysis table filtration using different parameters
  filter_df <- reactive({
    source("Scripts/Filter.R")
    return(Filtered(df(),input$pvalue,input$fdr,input$log2fc))
  })
  
### DGE Summary
  output$sum <- renderPlotly({
    source("Scripts/Summary_PieChart.R", local = TRUE)
    plt_pie(df())
  })
  
### RAW DGE Analysis table rendering
  output$table <- renderDT({datatable(df() , filter = "top",
                                      options = list(scrollX=TRUE),
                                      rownames = FALSE)})
  
### Filtered DGE Analysis table rendering
  output$table3 <- renderDT({
    validate(need(try(filter_df()),
                  "Please input either P-Value/FDR/All \n\n NOTE: Only fold change filteration is not permitted."))
    datatable(filter_df(), filter = "top", options = list(scrollX=TRUE), rownames = FALSE)
  })

### Download RAW DGE Output table in CSV format  
  output$downloadData <- downloadHandler(filename = function() {
    paste(input$comp,".csv", sep = "")
  },
  content = function(file) {
    write.csv(df(), file, row.names = FALSE)
  }
  )
  
### Download DGE filtered table in CSV format 
  output$downloadData2 <- downloadHandler(filename = function() {
    paste(input$comp,"_Filtered.csv", sep = "")
  },
  content = function(file) {
    write.csv(df(), file, row.names = FALSE)
  }
  )
  
### Volcano Plot 
  pvol <- eventReactive(input$plotvol,{
    source("Scripts/Volcano_Plot.R", local = TRUE)
    plotv <- pltvol(df())
    return(plotv)
  })
  
### Volcano Plot rendering
  output$volcano <- renderPlotly({ggplotly(pvol(),tooltip = c("text"))}) 
 
### Download Volcano Plot 
  output$downloadPlotvol <- downloadHandler(filename = function() {
    paste(input$comp,"_Volcano",".",input$vformat, sep = "")
  },
  content = function(file) {
    if (input$vformat == "png"){
      ggsave(file, pvol(),device = input$vformat,
             width = input$vol_width, height = input$vol_height, 
             limitsize = FALSE, bg= "white", dpi = 250)
    }else if (input$vformat == "jpeg"){
      ggsave(file, pvol(),device = input$vformat,
             width = input$vol_width, height = input$vol_height, 
             limitsize = FALSE, bg= "white", dpi = 250)
    }else if (input$vformat == "pdf"){
      ggsave(file, pvol(),device = input$vformat,
             width = input$vol_width, height = input$vol_height, 
             limitsize = FALSE, bg= "white", dpi = 250)
    }
  })
  
  ### Heatmaps
  output$heatmaps <- renderPlotly({
    if (input$strategy == "Top n"){
      source("Scripts/Heatmap.R")
      plt_heatmap(filter_df(),input$topn,input$log2fc)
    }
  })
}