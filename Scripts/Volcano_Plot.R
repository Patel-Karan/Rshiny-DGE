pltvol <- function(df){
  p <- df%>%ggplot(aes(x = log2FoldChange , y = -log10(pvalue),
                   text = paste('Gene-ID: ',GeneID,
                                '\nLog2FoldChange: ',log2FoldChange,
                                '\n-log10(pvalue) :',-log10(pvalue)),
                   col = Significant))+
    geom_point(size=1) +
    guides(color = guide_legend(override.aes = list(size = 3)))+
    labs(color = NULL,x = "log2 fold change", y = "-log10 p-value") +
    theme_minimal() +
    theme(text = element_text(face="bold"),
          axis.ticks = element_line(linewidth = 1),
          axis.text = element_text(colour = "black"),
          panel.border = element_rect(linetype = "solid", color = "black", size = 1, fill = NA),
          legend.justification = c("right", "top"),
          legend.background = element_rect(fill=NA,size = 0.25, linetype="solid", colour ="black")) +
    ggtitle("Volcano Plot") +
    geom_vline(xintercept = c(-0.6, 0.6), col = "#5B464B", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "#5B464B", linetype = 'dashed') +
    scale_color_manual(values = c("#B31B21", "#9897A9", "#1465AC"))
  return(p)
}
