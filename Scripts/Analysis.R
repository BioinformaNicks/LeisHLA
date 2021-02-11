#!/usr/bin/env Rscript

#################################################################
##                      Loading in Packages                    ##
#################################################################

dependencies <- c('dplyr', 'readr', 'stringr', 'vroom', 'ggplot2', 'ggsci', 'ggExtra', 'ggseqlogo', 'reticulate', 'here')
for (package in dependencies) {
  if (!package %in% installed.packages()) {
    install.packages(package)
  }
}
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(vroom))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(ggExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(ggseqlogo))
suppressMessages(library(reticulate))
suppressMessages(library(here))

main_folder <- paste0(str_remove(here(), 'Scripts'), '/')
output_folder <- paste0(str_remove(here(), 'Scripts'), '/Output/')

# This retrieves RScript commandline arguments from the bash script and saves them as a character vector.
args = commandArgs(trailingOnly=TRUE)
threshold = if (suppressWarnings(is.na(as.numeric(args[1])))) 0 else as.numeric(args[1])
protective_file <- args[2]
risk_file <- args[3]
stat_output_file <- paste0(output_folder, args[4])
plot_output_file <- args[5]

##################################################################
##                  Extracting data from input                  ##
##################################################################

acquire_strong_binders <- function(f_threshold, filename) {
  #' This function loads in the standardized epitope prediction data and then extracts the number of total predicted epitopes per HLA
  #' This function then extracts strong binding epitopes.
  #' 
  #' @description Loads in standardized epitope prediction data, extracts # of total predicted epitopes/HLA and the strong binders.
  #' 
  #' @param f_threshold An integer or float representing the relative affinity threshold with which to identify strong binding epitopes.
  #' @param filename A character string that contains the filename of the standardized epitope prediction data.
  #' 
  #' @usage acquire_strong_binders(f_threshold, filename)
  #' @return A list containing a vector with the number of total predicted epitopes/HLA and a tibble dataframe with the strong binding epitopes.
  
  total_result <- vroom(paste0(main_folder, filename), col_names = T, col_types = "fcccdd")
  colnames(total_result) <- c("HLA", "Peptide", "Binding_Core", "Identity", "Rank", "Affinity")
  
  epitope_length_vector <- c()
  strong_binders_results <- vector(mode = "list", length = length(levels(total_result$HLA)))
  strong_binding_number_vec <- c()
  
  for (Allele in 1:length(levels(total_result$HLA))) {
    epitopes <- filter(total_result, HLA == levels(total_result$HLA)[Allele])
    epitope_length_vector <- append(epitope_length_vector, nrow(epitopes))
    strong_binders_results[[Allele]] <- filter(epitopes, Affinity > f_threshold)
    strong_binding_number_vec <- append(strong_binding_number_vec, nrow(strong_binders_results[[Allele]]))
  }
  
  return_list <- list(epitope_length_vector, strong_binders_results, strong_binding_number_vec)
  return(return_list)
}

# Run the function for the protective data, saving the output in a temporary list
temp_protective_results <- acquire_strong_binders(threshold, protective_file)
# Extracting the contents of the temporary list: the number of total predicted epitopes, and a tibble df containing strong binders
protective_total_epitopes <- temp_protective_results[[1]]
protective_strong_binders_per_allele <- temp_protective_results[[2]]
protective_strong_binders_df <- bind_rows(protective_strong_binders_per_allele)
protective_num_strongbinders <- temp_protective_results[[3]]
rm(temp_protective_results)

# Run the function for the risk data, saving the output in a temporary list
temp_risk_results <- acquire_strong_binders(threshold, risk_file)
# Extracting the contents of the temporary list: the number of total predicted epitopes, and a tibble df containing strong binders
risk_total_epitopes <- temp_risk_results[[1]]
risk_strong_binders_per_allele <- temp_risk_results[[2]]
risk_strong_binders_df <- bind_rows(risk_strong_binders_per_allele)
risk_num_strongbinders <- temp_risk_results[[3]]
rm(temp_risk_results)

##################################################################
##        Annotating StrongBinders with Protein function        ##
##################################################################

#This code chunk and python script is only neccessary if running script on Leishmania data
source_python(here('Protein_Function_Annotator.py'))
protective_strong_binders_df <- protective_annotator(protective_strong_binders_df)
risk_strong_binders_df <- risk_annotator(risk_strong_binders_df)
rm(lbrazil_fasta_dict, ldonovani_fasta_dict, linfantum_fasta_dict, lmajor_fasta_dict, lmexicana_fasta_dict)

##################################################################
##      Identify common epitopes between Protective and Risk    ##
##################################################################

#Join together protective and risk data to identify common epitopes between them
joined_strongbinders_df <- inner_join(protective_strong_binders_df, risk_strong_binders_df, by = 'Peptide')
joined_strongbinders_df <- select(joined_strongbinders_df, 1, 3, 4, 5, 6, 7, 2, 8, 9, 10, 11, 12, 13)
colnames(joined_strongbinders_df) <- c("HLA_Protective", "BindingCore_Protective", "Identity_Protective", "Rank_Protective", "Affinity_Protective", "Protein_Protective",
                                       "Peptide",
                                       "HLA_Risk", "BindingCore_Risk", "Identity_Risk", "Rank_Risk", "Affinity_Risk", "Protein_Risk")

#filter out common epitopes between protective and risk
uniq_protective_strong_binders_df <- anti_join(protective_strong_binders_df, risk_strong_binders_df, by = 'Peptide')
uniq_risk_strong_binders_df <- anti_join(risk_strong_binders_df, protective_strong_binders_df, by = 'Peptide')

# Save the protective strong-binders to a file
vroom_write(protective_strong_binders_df,
            paste0(output_folder, str_split(protective_file, "_", simplify = TRUE)[1], '_Protective_StrongBinders.txt'),
            delim = " ", col_names = TRUE)

# Save the risk strong-binders to a file
vroom_write(risk_strong_binders_df,
            paste0(output_folder, str_split(risk_file, "_", simplify = TRUE)[1], '_Risk_StrongBinders.txt'),
            delim = " ", col_names = TRUE)

# Save the common strong-binders to a file
vroom_write(joined_strongbinders_df,
            paste0(output_folder, str_split(risk_file, "_", simplify = TRUE)[1], '_Common_StrongBinders.txt'),
            delim = " ", col_names = TRUE)

#################################################################
##        Identifying conserved epitopes across alleles        ##
#################################################################

conservancy_across_alleles <- function(uniq_strong_binders) {
  #'  This function identifies epitopes that bind across multiple alleles
  #' 
  #' @description Identifies epitopes that binds across multiple alleles
  #' 
  #' @param uniq_strong_binders Dataframe containing the epitope data
  #' 
  #' @usage conservancy_across_alleles(uniq_strong_binders)
  #' @return A dataframe containing the epitopes that bind across multiple alleles
  
  between_allele_conservancy_df <- uniq_strong_binders %>% 
    group_by(HLA, Peptide) %>%
    slice(which.max(Affinity)) %>%
    ungroup() %>%
    mutate(HLA = case_when(str_starts(HLA, 'DRB1_15') ~ 'DRB1_15',
                           str_starts(HLA, 'DRB1_16') ~ 'DRB1_16',
                           str_starts(HLA, 'DRB1_01') ~ 'DRB1_01',
                           str_starts(HLA, 'HLA-DQA10509') ~ "HLA-DQA10509-DQB103",
                           str_starts(HLA, 'HLA-DQA10602') ~ "HLA-DQA10602-DQB10302",
                           str_starts(HLA, 'HLA-DPA10103') ~ "HLA-DPA10103-DPB10401",
                           str_starts(HLA, 'HLA-DPA10401') ~ "HLA-DPA10401-DPB10101",
                           str_starts(HLA, 'DRB1_11') ~ 'DRB1_11',
                           str_starts(HLA, 'DRB1_13') ~ 'DRB1_13',
                           str_starts(HLA, 'DRB1_14') ~ 'DRB1_14')) %>% 
    group_by(HLA, Peptide) %>%
    slice(which.max(Affinity)) %>%
    ungroup() %>% 
    add_count(Peptide) %>%
    filter(n > 1)
  
  return(between_allele_conservancy_df)
}

protective_conserved_epitopes <- conservancy_across_alleles(uniq_protective_strong_binders_df)
risk_conserved_epitopes <- conservancy_across_alleles(uniq_risk_strong_binders_df)


##################################################################
##                      Statistical Output                      ##
##################################################################

# Results are written to a .txt log file
prot_allele_names <- levels(protective_strong_binders_df$HLA)
risk_allele_names <- levels(risk_strong_binders_df$HLA)

write_lines(paste("Total number of strong-binding Protective epitopes:", sum(protective_num_strongbinders)), stat_output_file)
write_lines(paste("Total number of strong-binding Risk epitopes:", sum(risk_num_strongbinders)), stat_output_file, append=T)
write_lines('', stat_output_file, append=T)
write_lines(c("Number of strong-binding epitopes for each Protective-associated allele:", paste(prot_allele_names, protective_num_strongbinders)), stat_output_file, append=T)
write_lines('', stat_output_file, append=T)
write_lines(c("Number of strong-binding epitopes for each Risk-associated allele:", paste(risk_allele_names, risk_num_strongbinders)), stat_output_file, append=T)
write_lines('', stat_output_file, append=T)
write_lines(c("Number of strong-binding epitopes common to Protective and Risk: ", nrow(joined_strongbinders_df)), stat_output_file, append=T)
write_lines('', stat_output_file, append=T)
write_lines("List of protective epitopes conserved across alleles: ", stat_output_file, append=T)
vroom_write(protective_conserved_epitopes, stat_output_file, delim = " ", col_names = TRUE, append = TRUE)
write_lines('', stat_output_file, append=T)
write_lines("List of risk epitopes conserved across alleles: ", stat_output_file, append=T)
vroom_write(risk_conserved_epitopes, stat_output_file, delim = " ", col_names = TRUE, append = TRUE)

##################################################################
##              Relative Affinity Distribution Graph            ##
##################################################################

# transform data to be easily converted to a long format for ease of plotting and statistics
protective_affinities <- select(uniq_protective_strong_binders_df, "Affinity") %>% mutate(., Status = 'Protective')
risk_affinities <- select(uniq_risk_strong_binders_df, "Affinity") %>% mutate(., Status = 'Risk')

# Status needs to be a factor for statistics and plotting
protective_affinities$Status <- as.factor(protective_affinities$Status)
risk_affinities$Status <- as.factor(risk_affinities$Status)

# bind the groups together in a long format
affinity_comparison <- bind_rows(protective_affinities, risk_affinities)

# Output of the relative affinity distribution comparison

if (!is.na(plot_output_file)) {
  h <- ggplot(affinity_comparison, aes(x=Affinity, y = after_stat(density), fill=Status)) +
    geom_histogram(binwidth = 0.2, color = 'black', position=position_dodge(0.16))
  
  h_plotdata <- ggplot_build(h)$data[[1]]
  h_plotdata$Status <- as.factor(h_plotdata$group)
  levels(h_plotdata$Status) <- levels(affinity_comparison$Status)
  
  distplot <- ggplot(h_plotdata, aes(x=x, y=y, fill = Status)) +
    geom_bar(color = 'black', stat = "identity") +
    geom_text(data = h_plotdata, aes(label = count), nudge_y = 0.125, fontface = 'bold') +
    ggtitle('NetMHCIIpan Relative Affinity Distribution') + 
    labs(y = 'Density', x='Relative Affinity') +
    theme_bw() + 
    theme(axis.text.y = element_text(color='black', size=11),
          axis.title.y = element_text(color='black', size=14),
          axis.text.x = element_text(color='black', size=11),
          axis.title.x = element_text(color='black', size=14),
          plot.title = element_text(size=16, hjust=0.5),
          legend.background = element_rect(colour='black', fill='white', linetype='solid'),
          legend.title = element_text(size=14, face='bold'),
          legend.text = element_text(size=11),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_continuous(breaks=seq(1, round(max(affinity_comparison$Affinity)/2,1)*2, by = 0.2)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0,NA)) +
    scale_fill_npg()
  
  #save to file
  ggsave(paste0(output_folder, plot_output_file), plot = distplot, width = 14)
}

#################################################################
##             SeqLogo plots RvP over all alleles              ##
#################################################################

#create named list
sequence_list <- list(Protective = uniq_protective_strong_binders_df$Binding_Core,
                      Risk = uniq_risk_strong_binders_df$Binding_Core)

#create seqlogo plot
RvP_seqlogo <- ggseqlogo::ggseqlogo(sequence_list, ncol=2) + 
  ggtitle(paste(str_split(plot_output_file, "_", simplify = TRUE)[1], "Protective vs Risk SeqLogo")) +
  labs(y = 'Bits', x = 'Position', fill = 'Amino Acid Properties') +
  guides(fill=guide_legend(title="Amino Acid Properties")) +
  theme(strip.text.x = element_text(size=14, face='bold'),
        strip.background.x =  element_blank(),
        axis.text.y = element_text(color='black', size=12),
        axis.title.y = element_text(color='black', size=16),
        axis.text.x = element_text(color='black', size=12, face = 'bold'),
        axis.title.x = element_text(color='black', size=16),
        plot.title = element_text(size=16, hjust=0.5),
        legend.background = element_rect(colour='black', fill='white'),
        legend.title = element_text(size=14, face='bold'),
        legend.text = element_text(size=11))

#save to file
ggsave(paste0(output_folder, str_split(plot_output_file, "[.]", simplify = TRUE)[1], "_SeqLogo.png"), plot = RvP_seqlogo, width = 14)

#################################################################
##            SeqLogo plot each allele individually            ##
#################################################################

seqlogo_allele_plotter <- function(strongbinder_df, status, numcol) {

  individual_allele_list <- vector(mode = "list", length = length(levels(strongbinder_df$HLA)))
  for (Allele in 1:length(levels(strongbinder_df$HLA))) {
    individual_allele_list[[Allele]] <- filter(strongbinder_df, HLA == levels(strongbinder_df$HLA)[Allele]) %>%
      select(Binding_Core) %>% pull
  }
  names(individual_allele_list) <- levels(strongbinder_df$HLA)
  
  #create seqlogo plot
  allele_seqlogo <- ggseqlogo(individual_allele_list, ncol=numcol) + 
    ggtitle(paste(str_split(plot_output_file, "_", simplify = TRUE)[1], status, "SeqLogo")) +
    labs(y = 'Bits', x = 'Position', fill = 'Amino Acid Properties') +
    guides(fill=guide_legend(title="Amino Acid Properties")) +
    theme(strip.text.x = element_text(size=14, face='bold'),
          strip.background.x =  element_blank(),
          axis.text.y = element_text(color='black', size=12),
          axis.title.y = element_text(color='black', size=16),
          axis.text.x = element_text(color='black', size=12, face = 'bold'),
          axis.title.x = element_text(color='black', size=16),
          plot.title = element_text(size=16, hjust=0.5),
          legend.background = element_rect(colour='black', fill='white'),
          legend.title = element_text(size=14, face='bold'),
          legend.text = element_text(size=11))
  
  ggsave(paste0(output_folder, str_split(plot_output_file, "_", simplify = TRUE)[1], "_", status, "_SeqLogo.png"), plot = allele_seqlogo, width = 14)
}

if (length(levels(uniq_protective_strong_binders_df$HLA)) > 4) {
  seqlogo_allele_plotter(uniq_protective_strong_binders_df,
                         "Protective",
                         numcol = length(levels(uniq_protective_strong_binders_df$HLA))/2)
} else {
  seqlogo_allele_plotter(uniq_protective_strong_binders_df,
                         "Protective",
                         numcol = length(levels(uniq_protective_strong_binders_df$HLA)))
}

if (length(levels(uniq_risk_strong_binders_df$HLA)) > 4) {
  seqlogo_allele_plotter(uniq_risk_strong_binders_df,
                         "Risk",
                         numcol = length(levels(uniq_risk_strong_binders_df$HLA))/2)
} else {
  seqlogo_allele_plotter(uniq_risk_strong_binders_df,
                         "Risk",
                         numcol = length(levels(uniq_risk_strong_binders_df$HLA)))
}

