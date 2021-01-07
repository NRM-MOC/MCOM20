# This script generates a pdf appendix based on a list of MoCiS objects
#

library(MoCiS)
library(MoCiS.tools)
library(tidyverse)
library(readxl)
library(knitr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)
library(ggspatial)
library(scales)

#setwd("Y:/MFO_Privat/MG/Milj?gifts?vervakning/MCOM/MCOM2020/R/MCOM20_report/")
setwd("Y:/MFO_Privat/MG/MOCIS/MCOM/R/MCOM20_report/")

source("functions/functions.R")

# loads "mocis_all" a list of MoCiS objects
load("../MCOM20_MoCiS/mocis_all.Rdata")

# Generate plots and tables for the objects in mocis_all
# Figures saved in subdirectory "figs" and tables in "tables"
walk(mocis_all, plot_ts)
walk(mocis_all, plot_mean_maps)
walk(mocis_all, plot_rate_maps)
walk(mocis_all, save_tables)

# Construct a table of contents for the appendix
contents <- map_df(mocis_all, mocis_contents) %>% 
  distinct() %>% 
  left_join(variables()) %>% 
  filter(!((group == "Biological data") & gen == "MYTI")) %>% # Exclude biological data for mussels
  arrange(group_order, var_order) %>% 
  select(gen, var, var_name, group, group_order, var_order) %>%
  distinct() %>% 
  na.omit()

# Generates an Rmd
file_name <- "appendix_draft_2.Rmd"
appendix(contents, file_name, 
         exclude_maps = c("CR", "AHCH", "LINDA", "PFHXA", "PFHPA", "LPFBS", "LPFDS", "PFDS"))
#renders to pdf
rmarkdown::render(file_name)
