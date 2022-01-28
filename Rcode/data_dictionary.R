###################################################################
# File: data_dictionary.R                             
#                                                                 
# Purpose: Creating data dictionary of features used in analysis                                 
#                                                                 
#                                                                 
# Author: A.Oliver				                                        
# Date: 1/20/21						                                        
#                                                                 
# Inputs (1):                                                     
# (1) Be able to source merging_features_new_data script          
#      successfully 
# (2) Data_dictionary_raw_for_supplement.xlsx
#                                                                 
#                     #####################                       
#                                                                 
# Outputs (1):                                                        
# (1) supplemental_data_dictionary.txt - data dictionary for the
#      features used in this analysis
#                                                                 
# Usage: Run the entire script without changes.                   
###################################################################

###################################
## Data Dictionary for supplemental
###################################

## load libraries
library(tidyverse)
library(readxl)
library(janitor)
## set wd
setwd("/home/datasets/new_datasets/")

## source all the features used in the analysis if needed
#source(file = "/home/scripts/merging_features_new_data.R")

## load up raw data dictionary with all features
data_dictionary <- read_excel(path = "/home/datasets/from_yasmine/Data_dictionary_raw_for_supplement.xlsx")
data_dictionary$`Variable / Field Name` <- make_clean_names(data_dictionary$`Variable / Field Name`)

## select out features that are used in ML
paper_dictionary <- subset(data_dictionary, data_dictionary$`Variable / Field Name` %in% colnames(for_directed_hypothesis_testing))
paper_dictionary_clean <- paper_dictionary %>% 
  rename(., "feature" = "Variable / Field Name") %>% 
  select(., feature, `Field Label`) %>% 
  clean_names()

## write to file
write.table(x = paper_dictionary_clean, 
            file = "/home/Manuscript/supplemental_data_dictionary.txt", 
            sep = "\t", quote = F, row.names = F)

## what are we missing
#missing_features <- base::setdiff(colnames(for_directed_hypothesis_testing), paper_dictionary$feature)
#View(as.data.frame(missing_features))

