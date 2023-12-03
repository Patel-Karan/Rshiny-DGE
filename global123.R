required_packages = c("shiny",
                      "DT",
                      "shinythemes", 
                      "shinycssloaders", 
                      "plotly", 
                      "dplyr", 
                      "tidyr", 
                      "reshape2", 
                      "tibble", 
                      "ggplot2", 
                      "ggdendro" ,
                      "deseq2")

missing_packages =
  required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages)) {
  cat(paste(
    "Missing packages:",
    paste(missing_packages, collapse = ";"),
    "\nAttempting to install them."
  ))
  install.packages(missing_packages)
}

source("Scripts/dexp_deseq2.R",local = TRUE)
source("Scripts/dexp_hclust.R", local = TRUE)
source("Scripts/Then.R", local = TRUE)
source("Scripts/Volcano_Plot.R", local = TRUE)