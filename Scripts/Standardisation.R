#!/usr/bin/env Rscript

#################################################################
##                      Loading in Packages                    ##
#################################################################

dependencies <- c('dplyr', 'readr', 'vroom', 'here', 'forcats', 'stringr')
for (package in dependencies) {
  if (!package %in% installed.packages()) {
    install.packages(package)
  }
}

suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(vroom))
suppressMessages(library(here))

main_folder <- paste0(str_remove(here(), 'Scripts'), '/')

# This retrieves RScript commandline arguments from the bash script and saves them as a character vector.
args = commandArgs(trailingOnly=TRUE)

#################################################################
##                        Loading in Data                      ##
#################################################################

prediction_result <- vroom(paste0(main_folder, args[2]), col_names=T, col_types = "fcccdd", altrep = F)
human_result <- vroom(paste0(main_folder, args[1]), col_names=T, col_types = "fcccdd")

#rename colnames
colnames(prediction_result) <- c("HLA", "Peptide", "Binding_Core", "Identity", "Rank" , "Affinity")
colnames(human_result) <- c("HLA", "Peptide", "Binding_Core", "Identity", "Rank" , "Affinity")

#################################################################
##                  Selecting only needed HLA's                ##
#################################################################

#Select only those HLA's from the reference that are present in the data you want to standardize, saves memory
human_result <- filter(human_result, HLA %in% levels(prediction_result$HLA))
#drop unused factors
human_result$HLA <- forcats::fct_drop(human_result$HLA)

##################################################################
##                    Normalizing the scores                    ##
##################################################################

# Only select epitopes with %Rank lower than 2%
prediction_result <- filter(prediction_result, Rank <= 2)
#human_result <- filter(human_result, Rank <= 2)

# Only select unique epitopes for a certain HLA. The same epitope can be derived from another protein, however.
prediction_result <- distinct(prediction_result, HLA, Peptide, Identity, .keep_all = T)

#create empty list to save data in for each HLA
allele_list <- vector(mode = "list", length = length(levels(prediction_result$HLA)))

#Loop over all alleles in the dataset, calculate the minimum score and then calculate the relative affinity for each epitope
#Save it all in the list for each HLA separately
for (Allele in 1:length(levels(prediction_result$HLA))) {
  min_score <- filter(human_result, HLA == levels(prediction_result$HLA)[Allele]) %>%
    slice(which.min(Affinity)) %>%
    select(Affinity) %>%
    pull #change affinity

  allele_specific_predictions <- filter(prediction_result, HLA == levels(prediction_result$HLA)[Allele])
  allele_specific_predictions$Affinity <- (min_score / allele_specific_predictions$Affinity)
  allele_list[[Allele]] <- allele_specific_predictions
  rm(allele_specific_predictions)
}

rm(prediction_result, human_result)

#Bind the individual dataframes containing the separate HLA data together to a single dataframe
prediction_result <- bind_rows(allele_list)

rm(allele_list)

#Write the standardized results to the output file provided in the bash script
vroom_write(prediction_result, paste0(main_folder, args[2]), delim = " ", col_names = TRUE)


