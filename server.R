library(shiny)
library(dplyr)
library(tidyr)
library(reshape2)
library(tibble)
library(DT)
#library(xlsx)
library(plotly)
library(ggplot2)
library(ggdendro)

options(shiny.maxRequestSize=30*1024^2)
options(java.parameters = "-Xmx30g")



function(input, output, session) {
  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    tbl <- read.table(inFile$datapath, sep=",", header = TRUE, 
                      row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
    return(tbl)
  })
  
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
      d_filter <- data() %>% select(c(input$control,input$case))
      nC <- length(input$control)
      nT <- length(input$case)
      source("Scripts/dexp_deseq2.R", local = TRUE)
      dexp <- d_exp(d_filter,nC,nT)
      dexp <- dexp %>% rownames_to_column(var = "GeneID")
      return(dexp)
      })
  
  output$table <- renderDT({datatable(df() , filter = "top", 
                                              options = list(scrollX=TRUE, 
                                                             extensions = 'Buttons',
                                                             buttons = c('copy', 'csv', 'excel', 'pdf')),
                                              rownames = FALSE)})
  
  output$downloadData <- downloadHandler(filename = function() {
      paste(input$comp,".",input$fformat, sep = "")
    },
    content = function(file) {
      if (input$fformat == "csv"){
        write.csv(df(), file, row.names = FALSE)
      }else if (input$fformat == "xlsx"){
        write.xlsx(df(),file, sheetNames="Sheet1", row.names = FALSE, append = FALSE)
      }else if (input$fformat == "txt"){
        write.table(df(), file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  )
  
  pvol <- eventReactive(input$plotvol,{
    df2 <- df() %>% filter(log2FoldChange != "NA")
    df2$Significant <- "Not-Significant"
    df2$Significant[df2$log2FoldChange > 1 & df2$pvalue < 0.05] <- "Upregulated"
    df2$Significant[df2$log2FoldChange < -1 & df2$pvalue < 0.05] <- "Downregulated"
    source("Scripts/Volcano_Plot.R", local = TRUE)
    plotv <- pltvol(df2)
    return(plotv)
  })
  
  output$volcano <- renderPlotly({ggplotly(pvol(),tooltip = c("text"))}) 
  
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
  
  observeEvent(data(),{
    updateSelectInput(session, "hclust_Samples",choices = append("All",names(data()))) 
  })
  
  df_hclust <- eventReactive(input$run_hclust,{
    source("Scripts/Then.R")
    validate(
      need(try(input$file1), "Please upload the Count data in CSV format in Analysis Section!") %then%
      need(try(input$hclust_Samples), "Please select either All or Samples !")
    )
    if ((length(input$hclust_Samples)==1) & ("All" %in% input$hclust_Samples)){
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
    }else if((length(input$hclust_Samples)>1) & ("All" %in% input$hclust_Samples)){
      dexp_hclust <- NULL
    }else{
      d_filter <- data() %>% select(c(input$hclust_Samples))
      if (length(names(d_filter))/2 == 0){
        nC <- length(names(d_filter))/2
        nT <- length(names(d_filter))-nC
        source("Scripts/dexp_hclust.R", local = TRUE)
        dexp_hclust <- d_exp(d_filter,nC,nT,input$hmethod,input$dmethod)
      }else{
        nC <- round(length(names(d_filter))/2, digits = 0)
        nT <- length(names(d_filter))-nC
        source("Scripts/dexp_hclust.R", local = TRUE)
        dexp_hclust <- d_exp(d_filter,nC,nT,input$hmethod,input$dmethod)
      }
    }
    return(dexp_hclust)
  })
  
  output$table2 <- renderDT({
    validate(
      need(df_hclust()$d, 'Please select either All or Samples in multiples!'),
    )
    datatable( df_hclust()$distance, filter = "top", options = list(scrollX=TRUE, rownames = FALSE))
    })
  
  output$dendoplot <- renderPlotly({
    validate(
      need(df_hclust()$d, 'Please select either All or Samples in multiples!'),
    )
    dendop <- ggdendrogram(df_hclust()$model$m, rotate = FALSE, leaf_labels = FALSE, size = 2)
    ggplotly(dendop)
  })
  
  output$pca <- renderPlotly(({
    validate(
      need(df_hclust()$d, 'Please select either All or Samples in multiples!'),
    )
    scores2 = df_hclust()$pca$pca_score
    pca_plot <- plot_ly(scores2, x =scores2$PC1 , y = scores2$PC2, 
                        text = rownames(scores2),mode = "markers",marker = list(size = 11))
    pca_plot <- layout(pca_plot,
                       xaxis = list(title=paste("PC1 (", df_hclust()$pca$p_percent[1], "%)", sep = "")),
                       yaxis = list(title=paste("PC2 (", df_hclust()$pca$p_percent[2], "%)", sep = "")))
    ggplotly(pca_plot, tooltip = c("text"))
  }))
  
  output$boxplot <- renderPlotly({
    if (input$hclust_Samples == "All"){
      box_plot <- df_hclust()$model$df %>% 
        pivot_longer(1:length(names(data())),names_to="Sample",values_to = "RNC") %>%
        ggplot(aes(x=Sample,y=-log10(RNC),fill=Sample))+
        geom_boxplot(alpha=0.3)+
        theme_bw()+
        labs(x="Samples",y="Raw Normalized Counts")+
        theme(legend.title = element_blank(),plot.title = element_text(size=20, face="bold"))
      ggplotly(box_plot)
    }else{
      box_plot <- df_hclust()$model$df %>% 
        pivot_longer(1:length(names(input$hclust_Samples)),names_to="Sample",values_to = "RNC") %>%
        ggplot(aes(x=Sample,y=-log10(RNC),fill=Sample))+
        geom_boxplot(alpha=0.3)+
        theme_bw()+
        labs(x="Samples",y="Raw Normalized Counts")+
        theme(legend.title = element_blank(),plot.title = element_text(size=20, face="bold"))
      ggplotly(box_plot)
    }
    
  })
}