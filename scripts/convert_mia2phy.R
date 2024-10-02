#Convert mia to phyloseq... File is generated for Oleg

#Run this at cluster with following commands:
#grun.py -n convert_mia2phy -q highmem.q -c "source $HOME/.bashrc; conda activate finriskR_env; Rscript scripts/convert_mia2phy.R"

#Libraries
library(data.table)
suppressPackageStartupMessages(library(tidyverse))
library(phyloseq)
library(microbiome)


phy_out <- "data_mg/phy_gg2_MGS.rds"

#Paths from external file
path.list <- fread("paths.txt") %>% 
  group_by(path_name) %>%
  transmute(named_vec = list(path)) %>%
  deframe()

phy <- readRDS(path.list$phymg_file) %>% mia::makePhyloseqFromTreeSummarizedExperiment()

saveRDS(phy, phy_out)

