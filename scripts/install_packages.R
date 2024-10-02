
#Library block from the Rmd - can be run separately

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager", repos='https://ftp.acc.umu.se/mirror/CRAN/')
#BiocManager::install()

#Matti styled definition:
packages <- c("data.table","knitr", "ggplot2",
              "phyloseq", "microbiome","vegan","biomformat","mia",
              "survival","survminer",
              "cmprsk", "glmnet","DESeq2","RColorBrewer","cowplot","ggpubr","ggrepel",
              "scales","kableExtra","tableone", "tidyverse","ANCOMBC")

is.installed <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, ask=F)
  }
  sapply(pkg, require, character.only = TRUE)
}
is.installed(packages)
