#this script downloads the key marker genes of L. major specified in the Dataset S01 in the Supporting Material of Catta-Preta et al 2024.
#it combines them in a single data frame to make it easier to use in downstream analysis.
#link: https://www.pnas.org/doi/10.1073/pnas.2406776121#supplementary-materials

#clean the environment####
rm(list = ls(all.names = TRUE))
gc()

#libraries####
library(tidyverse)
library(readxl)

#set links and paths
markers_link <- "https://www.pnas.org/doi/suppl/10.1073/pnas.2406776121/suppl_file/pnas.2406776121.sd01.xlsx"
dest_file <- "../inputs/external_data/pnas.2406776121.sd01.xlsx"

#download the list
#download.file(markers_link, destfile = dest_file) can't download directly (error 403: forbidden). Download manually instead.

#import it
sheets <- excel_sheets(dest_file)
marker_tables <- list()
for(sheet in sheets){
  if(sheet == "Cluster 1 - Procyclic"){
    marker_tables[[sheet]] <- read_excel(dest_file, sheet = sheet, col_names = TRUE, skip = 2)  # first sheet contains two extra rows that must be ignored 
  }else{
    marker_tables[[sheet]] <- read_excel(dest_file, sheet = sheet, col_names = TRUE)  
  }
  
}

#merge it into a single data frame
markers_table <- bind_rows(marker_tables, .id = "cluster")

#export the table
write_tsv(markers_table, "../inputs/external_data/markers_from_Lmj_atlas_formatted.tsv")
