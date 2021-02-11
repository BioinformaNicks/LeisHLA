
#################################################################
##                      Loading in Packages                    ##
#################################################################
dependencies <- c('tidyr', 'tibble', 'dplyr', 'readr', 'stringr', 'vroom',
                  'here', 'lmerTest', 'forcats')
for (package in dependencies) {
  if (!package %in% installed.packages()) {
    install.packages(package)
  }
}

suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(vroom))
suppressMessages(library(lmerTest))
suppressMessages(library(here))

#################################################################
##                    Setting folder settings                   #
#################################################################

main_folder <- paste0(str_remove(here(), 'Scripts'), '/')
data_folder <- paste0(str_remove(here(), 'Scripts'), '/Data/')
output_folder <- paste0(str_remove(here(), 'Scripts'), '/Output/')

#################################################################
##                      Set files to open                       #
#################################################################

protective_files <- c("LBraziliensis_Protective_StrongBinders.txt", "LDonovani_Protective_StrongBinders.txt",
                      "LInfantum_Protective_StrongBinders.txt", "LMexicana_Protective_StrongBinders.txt")

risk_files <- c("LBraziliensis_Risk_StrongBinders.txt", "LDonovani_Risk_StrongBinders.txt",
                "LInfantum_Risk_StrongBinders.txt", "LMexicana_Risk_StrongBinders.txt")

#################################################################
##                        Loading in data                       #
#################################################################

strong_binder_loading <- function(filename) {
  #' This function loads in the strong-binding epitope data
  #' 
  #' @description Loads in the strong-binding epitope data and returns it to environment
  #'
  #' @param filename character vector containing filenames you wish to load
  #' 
  #' @usage strong_binder_loading(filename)
  #' @return returns a list containing dataframes of the strong-binding epitope data in the different files

  strong_binders <- vector(mode = "list", length = 4)
  for (file in 1:length(strong_binders)) {
    strong_binders[[file]] <- vroom(paste0(main_folder, filename[file]), col_names = T, col_types = "fcccddc")
    names(strong_binders)[file] <- paste0(str_split(filename[file], "_", simplify = TRUE)[1], "_", str_split(filename[file], "_", simplify = TRUE)[2])
  }
  return(strong_binders)
}

protective <- strong_binder_loading(protective_files)

risk <- strong_binder_loading(risk_files)

# suppressMessages(list2env(protective,envir=.GlobalEnv))
# suppressMessages(list2env(risk,envir=.GlobalEnv))

#This unlists the strong-binding epitope dataframes to the object environment
LBraziliensis_Protective <- enframe(protective) %>% unnest(cols = c(value)) %>% filter(., name == 'LBraziliensis_Protective') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))
LBraziliensis_Risk <- enframe(risk) %>% unnest(cols = c(value)) %>% filter(., name == 'LBraziliensis_Risk') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))

LDonovani_Protective <- enframe(protective) %>% unnest(cols = c(value)) %>% filter(., name == 'LDonovani_Protective') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))
LDonovani_Risk <- enframe(risk) %>% unnest(cols = c(value)) %>% filter(., name == 'LDonovani_Risk') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))

LInfantum_Protective <- enframe(protective) %>% unnest(cols = c(value)) %>% filter(., name == 'LInfantum_Protective') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))
LInfantum_Risk <- enframe(risk) %>% unnest(cols = c(value)) %>% filter(., name == 'LInfantum_Risk') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))

LMexicana_Protective <- enframe(protective) %>% unnest(cols = c(value)) %>% filter(., name == 'LMexicana_Protective') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))
LMexicana_Risk <- enframe(risk) %>% unnest(cols = c(value)) %>% filter(., name == 'LMexicana_Risk') %>% select(!name) %>% mutate(HLA = forcats::fct_drop(HLA))

##################################################################
##                      Duplicate cleaning                       #
##################################################################

duplicate_epitope_filtering <- function(protective_data, risk_data) {
  #' This function filters out any duplicate epitopes, keeping only unique ones with the highest affinity
  #' additionally, it filters out epitopes common between the protective and risk data
  #'
  #' @description This function filters out any duplicate epitopes, keeping only unique ones
  #' 
  #' @param protective_data Dataframe containing the data you wish to rid of duplicates
  #' @param risk_data Dataframe containing data you wish to filter out of protective_data
  #' 
  #' @usage duplicate_epitope_filtering(protective_data, risk_data)
  #' @return a dataframe containing filtered data to be used in conservancy across alleles check

   filtered_data <- anti_join(protective_data, risk_data, by = 'Peptide') %>%
     group_by(HLA, Peptide) %>%
     slice(which.max(Affinity)) %>%
     ungroup()

   return(filtered_data)

}

uniq_LBraziliensis_Protective <- duplicate_epitope_filtering(LBraziliensis_Protective, LBraziliensis_Risk)

uniq_LBraziliensis_Risk <- duplicate_epitope_filtering(LBraziliensis_Risk, LBraziliensis_Protective) %>%
  filter(., HLA %in% c("HLA-DQA10602-DQB10302", "HLA-DQA10509-DQB10302")) %>%
  mutate(HLA = forcats::fct_drop(HLA)) %>%
  mutate(HLA = as.factor(case_when(str_starts(HLA, 'HLA-DQA1') ~ 'DQB1_0302'))) %>%
  group_by(HLA, Peptide) %>%
  slice(which.max(Affinity)) %>%
  ungroup()

uniq_LDonovani_Protective <- duplicate_epitope_filtering(LDonovani_Protective, LDonovani_Risk) %>%
  filter(., HLA %in% c('DRB1_0101', 'DRB1_1602')) %>%
  mutate(HLA = forcats::fct_drop(HLA))

uniq_LDonovani_Risk <- duplicate_epitope_filtering(LDonovani_Risk, LDonovani_Protective) %>%
  filter(., HLA %in% c('DRB1_1101', 'DRB1_1104', 'DRB1_1301', 'DRB1_1302', 'DRB1_1404')) %>%
  mutate(HLA = forcats::fct_drop(HLA)) %>%
  mutate(HLA = as.factor(case_when(str_starts(HLA, 'DRB1_11') ~ 'DRB1_11',
                                   str_starts(HLA, 'DRB1_1301') ~ 'DRB1_1301',
                                   str_starts(HLA, 'DRB1_1302') ~ 'DRB1_1302',
                                   str_starts(HLA, 'DRB1_1404') ~ 'DRB1_1404'))) %>%
  group_by(HLA, Peptide) %>%
  slice(which.max(Affinity)) %>%
  ungroup()

uniq_LInfantum_Protective <- duplicate_epitope_filtering(LInfantum_Protective, LInfantum_Risk) %>%
  filter(., HLA %in% c('DRB1_0101', 'DRB1_1501', 'DRB1_1502', 'DRB1_1602')) %>%
  mutate(HLA = forcats::fct_drop(HLA))

uniq_LInfantum_Risk <- duplicate_epitope_filtering(LInfantum_Risk, LInfantum_Protective) %>%
  filter(., HLA %in% c('DRB1_1101', 'DRB1_1104', 'DRB1_1301', 'DRB1_1302', 'DRB1_1404')) %>%
  mutate(HLA = forcats::fct_drop(HLA)) %>%
  mutate(HLA = as.factor(case_when(str_starts(HLA, 'DRB1_11') ~ 'DRB1_11',
                                   str_starts(HLA, 'DRB1_1301') ~ 'DRB1_1301',
                                   str_starts(HLA, 'DRB1_1302') ~ 'DRB1_1302',
                                   str_starts(HLA, 'DRB1_1404') ~ 'DRB1_1404'))) %>%
  group_by(HLA, Peptide) %>%
  slice(which.max(Affinity)) %>%
  ungroup()

uniq_LMexicana_Protective <- duplicate_epitope_filtering(LMexicana_Protective, LMexicana_Risk)
uniq_LMexicana_Risk <- duplicate_epitope_filtering(LMexicana_Risk, LMexicana_Protective)

##################################################################
##          Selecting only top 100 epitopes per allele           #
##################################################################

top_n_selector <- function(data_to_select_from, n_epitopes) {
  #' This function selects only the top n epitopes per allele with the strongest affinity for that allele
  #' 
  #' @description selects only the top n epitopes per allele with the strongest affinity for that allele
  #' 
  #' @param data_to_select_from the strong-binding epitope data to select from
  #' @param n_epitopes the number n of epitopes you wish to select
  #' 
  #' @usage top_n_selector(data_to_select_from, n_epitopes)
  #' @return a dataframe containing the top n epitopes for each allele
  
  individual_allele_list <- vector(mode = "list", length = length(levels(data_to_select_from$HLA)))
  
  for (allele in 1:length(levels(data_to_select_from$HLA))) {
    individual_allele_list[[allele]] <- filter(data_to_select_from, HLA == levels(data_to_select_from$HLA)[allele]) %>% slice_max(., order_by = Affinity, n = n_epitopes, with_ties = FALSE)
  }
  
  top_n_per_allele <- bind_rows(individual_allele_list)
  return(top_n_per_allele)
}

LBraziliensis_Protective_top_100 <- top_n_selector(uniq_LBraziliensis_Protective, 100)
LBraziliensis_Risk_top_100 <- top_n_selector(uniq_LBraziliensis_Risk, 100)

LDonovani_Protective_top_100 <- top_n_selector(uniq_LDonovani_Protective, 100)
LDonovani_Risk_top_100 <- top_n_selector(uniq_LDonovani_Risk, 100)

LInfantum_Protective_top_100 <- top_n_selector(uniq_LInfantum_Protective, 100)
LInfantum_Risk_top_100 <- top_n_selector(uniq_LInfantum_Risk, 100)

LMexicana_Protective_top_100 <- top_n_selector(uniq_LMexicana_Protective, 100)
LMexicana_Risk_top_100 <- top_n_selector(uniq_LMexicana_Risk, 100)

#################################################################
##              Creation of Mixed Model-ready table             #
#################################################################

mm_ready_table_maker <- function(protective_data, risk_data) {
  #' this function creates a table that is ready for linear mixed effect modelling
  #' 
  #' @description creates a table that is ready for linear mixed effect modelling
  #' 
  #' @param protective_data the protective-associated strong-binding epitope data
  #' @param risk_data the protective-associated strong-binding epitope data
  #'
  #' @usage mm_ready_table_maker(protective_data, risk_data)
  #' @return a table that is ready for linear mixed effect modelling

  protective_ready <- mutate(protective_data, Status = 'Protective')
  risk_ready <- mutate(risk_data, Status = 'Risk')
  
  mm_table <- bind_rows(protective_ready, risk_ready) %>% select(HLA, Peptide, Affinity, Status) %>% mutate(Status = as.factor(Status))

  return(mm_table)
}

#here, only the top 100 epitopes were used as an example, but all epitopes were used in the publication

LBraziliensis_MM_table <- mm_ready_table_maker(LBraziliensis_Protective_top_100, LBraziliensis_Risk_top_100)
#scaling can be needed or tried out, as an example
LBraziliensis_MM_table[,3] <- scale(LBraziliensis_MM_table[,3])
LDonovani_MM_table <- mm_ready_table_maker(LDonovani_Protective_top_100, LDonovani_Risk_top_100)
LInfantum_MM_table <- mm_ready_table_maker(LInfantum_Protective_top_100, LInfantum_Risk_top_100)
LMexicana_MM_table <- mm_ready_table_maker(LMexicana_Protective_top_100, LMexicana_Risk_top_100)

across_species <- bind_rows(LBraziliensis_MM_table, LDonovani_MM_table, LInfantum_MM_table, LMexicana_MM_table)

#setting REML to T or F impacts the singularity of results (thus confidence in result), here set to true as an example
MM_LBraziliensis <- lmer(Affinity ~ Status + (1|HLA), data = LBraziliensis_MM_table)
MM_LDonovani <- lmer(Affinity ~ Status + (1|HLA), data = LDonovani_MM_table, REML = F)
MM_LInfantum <- lmer(Affinity ~ Status + (1|HLA), data = LInfantum_MM_table, REML = F)
MM_LMexicana <- lmer(Affinity ~ Status + (1|HLA), data = LMexicana_MM_table, REML = F)

