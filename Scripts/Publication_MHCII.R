
#################################################################
##                      Loading in Packages                    ##
#################################################################

dependencies <- c('dplyr', 'readr', 'stringr', 'vroom', 'ggplot2', 'ggsci', 'ggExtra', 'grid', 'ggseqlogo', 'here', 'UpSetR')
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
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(ggseqlogo))
suppressMessages(library(UpSetR))
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
  
  strong_binders <- vector(mode = "list", length = length(filename))
  for (file in 1:length(strong_binders)) {
    strong_binders[[file]] <- vroom(paste0(main_folder, filename[file]), col_names = T, col_types = "fcccddc")
    names(strong_binders)[file] <- paste0(str_split(filename[file], "_", simplify = TRUE)[1], "_", str_split(filename[file], "_", simplify = TRUE)[2])
  }
  return(strong_binders)
}

protective <- strong_binder_loading(protective_files)

risk <- strong_binder_loading(risk_files)

suppressMessages(list2env(protective,envir=.GlobalEnv))
suppressMessages(list2env(risk,envir=.GlobalEnv))

rm(protective, risk)

#################################################################
##                  Duplicate epitope filtering                 #
#################################################################

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
  
  uniq_epitope_df <- anti_join(protective_data, risk_data, by = 'Peptide') %>% 
    group_by(HLA, Peptide) %>%
    slice(which.max(Affinity)) %>%
    ungroup()
  
  return(uniq_epitope_df)
}

#################################################################
##          Create panel of affinity distribution plots         #
#################################################################

RvP_affinities <- function(protective_data, risk_data, plotname) {
  #' this function creates a barplot with the relative affinity distribution of the protective and risk strong-binding epitopes
  #' 
  #' @description creates a barplot with the relative affinity distribution of the protective and risk strong-binding epitopes
  #' 
  #' @param protective_data dataframe containing the protective strong-binding epitopes
  #' @param risk_data dataframe containing the risk strong-binding epitopes
  #' @param plotname title of plot
  #' 
  #' @usage RvP_affinities(protective_data, risk_data, plotname)
  #' @return a plot object containing a barplot of the relative affinity distribution
  
  protective_affinities <- select(protective_data, "Affinity") %>% mutate(., Status = 'Protective')
  risk_affinities <- select(risk_data, "Affinity") %>% mutate(., Status = 'Risk')
  
  # Status needs to be a factor for statistics and plotting
  protective_affinities$Status <- as.factor(protective_affinities$Status)
  risk_affinities$Status <- as.factor(risk_affinities$Status)
  
  # bind the groups together in a long format
  affinity_comparison <- bind_rows(protective_affinities, risk_affinities)
  
  # Output of the relative affinity distribution comparison
  h <- ggplot(affinity_comparison, aes(x=Affinity, y = after_stat(density), fill=Status)) +
    geom_histogram(binwidth = 0.2, color = 'black', position=position_dodge(0.16))
  
  h_plotdata <- ggplot_build(h)$data[[1]]
  h_plotdata$Status <- as.factor(h_plotdata$group)
  levels(h_plotdata$Status) <- levels(affinity_comparison$Status)
  
  h_plotdata <- as.data.frame(h_plotdata %>% group_by(Status) %>% mutate(MaxCount = sum(count)))
  h_plotdata$y <- h_plotdata$count / h_plotdata$MaxCount
  
  distplot <- ggplot(h_plotdata, aes(x=x, y=y, fill = Status)) +
    geom_col(color = 'black') +
    geom_text(data = h_plotdata, aes(label = count), nudge_y=0.0250) +
    ggtitle(plotname) + 
    labs(y = 'Proportion of epitopes', x='Relative Affinity') +
    theme_bw() + 
    theme(axis.text.y = element_text(color='black', size=12),
          axis.title.y = element_text(color='black', size=14),
          axis.text.x = element_text(color='black', size=12),
          axis.title.x = element_text(color='black', size=14),
          plot.title = element_text(size=16, hjust=0.5),
          legend.background = element_rect(colour='black', fill='white', linetype='solid'),
          legend.title = element_text(size=14, face='bold'),
          legend.text = element_text(size=11, face ='bold'),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_continuous(breaks=seq(1, round(max(affinity_comparison$Affinity)/2,1)*2, by = 0.2)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0,0.6)) +
    scale_fill_manual(values = c("#00A087FF", "#E64B35FF"))
  
  return(distplot)
  
}

LBraziliensis_RvP <- RvP_affinities(LBraziliensis_Protective, LBraziliensis_Risk, "L. braziliensis MHC-II Relative Affinity Distribution")
LDonovani_RvP <- RvP_affinities(LDonovani_Protective, LDonovani_Risk, "L. donovani MHC-II Relative Affinity Distribution")
LInfantum_RvP <- RvP_affinities(LInfantum_Protective, LInfantum_Risk, "L. infantum MHC-II Relative Affinity Distribution")
LMexicana_RvP <- RvP_affinities(LMexicana_Protective, LMexicana_Risk, "L. mexicana MHC-II Relative Affinity Distribution")

#create a panel of plots for publication
affinity_panel <- ggarrange(LBraziliensis_RvP, LDonovani_RvP, LInfantum_RvP, LMexicana_RvP, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

#ggsave(paste0(output_folder, 'Publication_ready_affinities.png'), affinity_panel, width = 22, height = 11)

#################################################################
##              Create panel of affinity boxplots               #
#################################################################

uniq_LBraziliensis_Protective <- duplicate_epitope_filtering(LBraziliensis_Protective, LBraziliensis_Risk)
uniq_LBraziliensis_Risk <- duplicate_epitope_filtering(LBraziliensis_Risk, LBraziliensis_Protective)

uniq_LDonovani_Protective <- duplicate_epitope_filtering(LDonovani_Protective, LDonovani_Risk)
uniq_LDonovani_Risk <- duplicate_epitope_filtering(LDonovani_Risk, LDonovani_Protective)

uniq_LInfantum_Protective <- duplicate_epitope_filtering(LInfantum_Protective, LInfantum_Risk)
uniq_LInfantum_Risk <- duplicate_epitope_filtering(LInfantum_Risk, LInfantum_Protective)

uniq_LMexicana_Protective <- duplicate_epitope_filtering(LMexicana_Protective, LMexicana_Risk)
uniq_LMexicana_Risk <- duplicate_epitope_filtering(LMexicana_Risk, LMexicana_Protective)

RvP_boxplot_affinities <- function(protective_data, risk_data, plotname, xaxis_angle, xaxis_vjust) {
  #' this function creates a boxplot with the relative affinity of the protective and risk strong-binding epitopes for each allele
  #' 
  #' @description creates a boxplot with the relative affinity of the protective and risk strong-binding epitopes for each allele
  #' 
  #' @param protective_data dataframe containing the protective strong-binding epitopes
  #' @param risk_data dataframe containing the risk strong-binding epitopes
  #' @param plotname title of plot
  #' @param xaxis_angle change the angle of the x-axis labels
  #' @param xaxis_vjust change the vertical position of the x-axis labels
  #' 
  #' @usage RvP_boxplot_affinities(protective_data, risk_data, plotname, xaxis_angle, xaxis_vjust)
  #' @return a plot object containing a boxplot of the relative affinity distribution across each allele

  protective_affinities <- protective_data %>% mutate(., Status = 'Protective')
  risk_affinities <- risk_data %>% mutate(., Status = 'Risk')
  
  # Status needs to be a factor for statistics and plotting
  protective_affinities$Status <- as.factor(protective_affinities$Status)
  risk_affinities$Status <- as.factor(risk_affinities$Status)
  # bind the groups together in a long format
  affinity_comparison <- bind_rows(protective_affinities, risk_affinities)
  
  # Output of the relative affinity distribution per allele
  distplot <- ggplot(affinity_comparison, aes(x=HLA, y = Affinity, fill=Status)) +
    geom_boxplot() +
    ggtitle(plotname) + 
    labs(y = 'Affinity', x='HLA') +
    theme_bw() + 
    theme(axis.text.y = element_text(color='black', size=11),
          axis.title.y = element_text(color='black', size=14),
          axis.text.x = element_text(color='black', size=11, angle = xaxis_angle, vjust=xaxis_vjust),
          axis.title.x = element_text(color='black', size=14),
          plot.title = element_text(size=16, hjust=0.5),
          legend.background = element_rect(colour='black', fill='white', linetype='solid'),
          legend.title = element_text(size=14, face='bold'),
          legend.text = element_text(size=11, face ='bold'),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color="#dadfe6", size=0.25)) +
    scale_fill_manual(values = c("#00A087FF", "#E64B35FF"))
  
  return(distplot)
  
}

LBraziliensis_affinity_boxplot <- RvP_boxplot_affinities(uniq_LBraziliensis_Protective, uniq_LBraziliensis_Risk, 'L. braziliensis Relative Affinity across HLA class II Alleles', 20, 0.6)
LDonovani_affinity_boxplot <- RvP_boxplot_affinities(uniq_LDonovani_Protective, uniq_LDonovani_Risk, 'L. donovani Relative Affinity across HLA class II Alleles', 0, 1)
LInfantum_affinity_boxplot <- RvP_boxplot_affinities(uniq_LInfantum_Protective, uniq_LInfantum_Risk, 'L. infantum Relative Affinity across HLA class II Alleles', 0, 1)
LMexicana_affinity_boxplot <- RvP_boxplot_affinities(uniq_LMexicana_Protective, uniq_LMexicana_Risk, 'L. mexicana Relative Affinity across HLA class II Alleles', 20, 0.6)

#create a panel of boxplots for publication
affinity_boxplot_panel <- ggarrange(LBraziliensis_affinity_boxplot, LMexicana_affinity_boxplot,
                                    LInfantum_affinity_boxplot, LDonovani_affinity_boxplot,
                                    ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

#ggsave(paste0(output_folder, 'Publication_ready_affinity_boxplot.png'), affinity_boxplot_panel, width = 22, height = 11)

#################################################################
##                Proteins unique to Protective                 #
#################################################################

proteins_unique_to_protective <- function(protective_data, risk_data) {
  #' This function filters out any duplicate epitopes, keeping only unique ones
  #' additionally, it filters out proteins common between the protective and risk data
  #'
  #' @description This function filters out any duplicate epitopes, keeping only unique ones
  #' 
  #' @param protective_data Dataframe containing the data you wish to rid of duplicates
  #' @param risk_data Dataframe containing data you wish to filter out of protective_data
  #' 
  #' @usage proteins_unique_to_protective(protective_data, risk_data)
  #' @return a dataframe containing filtered data to find proteins unique to protective
  
  uniq_protective <- anti_join(protective_data, risk_data, by = 'Identity') %>%
    distinct(HLA, Peptide, .keep_all = T) %>% distinct(Identity, .keep_all = T)
  
  return(uniq_protective)
}

#Total proteins unique to protective
uniq_proteins_LBraziliensis_Protective <- proteins_unique_to_protective(LBraziliensis_Protective, LBraziliensis_Risk)
uniq_proteins_LDonovani_Protective <- proteins_unique_to_protective(LDonovani_Protective, LDonovani_Risk)
uniq_proteins_LInfantum_Protective <- proteins_unique_to_protective(LInfantum_Protective, LInfantum_Risk)
uniq_proteins_LMexicana_Protective <- proteins_unique_to_protective(LMexicana_Protective, LMexicana_Risk)

uniq_proteins_LBraziliensis_Risk <- proteins_unique_to_protective(LBraziliensis_Risk, LBraziliensis_Protective)
uniq_proteins_LDonovani_Risk <- proteins_unique_to_protective(LDonovani_Risk, LDonovani_Protective)
uniq_proteins_LInfantum_Risk <- proteins_unique_to_protective(LInfantum_Risk, LInfantum_Protective)
uniq_proteins_LMexicana_Risk <- proteins_unique_to_protective(LMexicana_Risk, LMexicana_Protective)

#################################################################
##                Create panel of SeqLogo plots                 #
#################################################################

#Prep L. braziliensis to merge DQB1*03 alleles
merged_LBraziliensis_Risk <- LBraziliensis_Risk %>%
  mutate(HLA = as.factor(case_when(str_starts(HLA, 'HLA-DQA1') ~ 'DQB1_03'))) %>%
  group_by(HLA, Peptide) %>%
  slice(which.max(Affinity)) %>%
  ungroup()

#Prep L. donovani to merge DRB1*11
merged_LDonovani_Risk <- LDonovani_Risk %>%
  mutate(HLA = as.factor(case_when(str_starts(HLA, 'DRB1_11') ~ 'DRB1_11',
                                   str_starts(HLA, 'DRB1_1301') ~ 'DRB1_1301',
                                   str_starts(HLA, 'DRB1_1302') ~ 'DRB1_1302',
                                   str_starts(HLA, 'DRB1_1404') ~ 'DRB1_1404'))) %>%
  group_by(HLA, Peptide) %>%
  slice(which.max(Affinity)) %>%
  ungroup()

#Prep L. infantum to merge DRB1*11
merged_LInfantum_Risk <- LInfantum_Risk %>%
  mutate(HLA = as.factor(case_when(str_starts(HLA, 'DRB1_11') ~ 'DRB1_11',
                                   str_starts(HLA, 'DRB1_1301') ~ 'DRB1_1301',
                                   str_starts(HLA, 'DRB1_1302') ~ 'DRB1_1302',
                                   str_starts(HLA, 'DRB1_1404') ~ 'DRB1_1404'))) %>%
  group_by(HLA, Peptide) %>%
  slice(which.max(Affinity)) %>%
  ungroup()

#first aggregate data and filter out duplicates for each allele

Protective_aggregated <- bind_rows(LBraziliensis_Protective, LDonovani_Protective, LInfantum_Protective, LMexicana_Protective) %>%
  distinct(HLA, Peptide, .keep_all = T)
Risk_aggregated <- bind_rows(merged_LBraziliensis_Risk, merged_LDonovani_Risk, merged_LInfantum_Risk, LMexicana_Risk) %>%
  distinct(HLA, Peptide, .keep_all = T)

rm(merged_LBraziliensis_Risk)

seqlogo_allele_plotter <- function(aggregated_data, status, numcol) {
  #' This function creates consensus sequence logos of the epitope binding cores for each allele
  #' 
  #' @description creates consensus sequence logos of the epitope binding cores for each allele
  #' 
  #' @param aggegrated_data a dataframe containing the strong-binding epitope data of each species aggregated
  #' @param status the association group: protective or risk
  #' @param numcol the number of columns you wish to use for creation of an allele matrix
  #' 
  #' @usage seqlogo_allele_plotter(aggregated_data, status, numcol)
  #' @return a plot object containing the consensus sequence logos
  
  individual_allele_list <- vector(mode = "list", length = length(levels(aggregated_data$HLA)))
  
  for (Allele in 1:length(levels(aggregated_data$HLA))) {
    individual_allele_list[[Allele]] <- filter(aggregated_data, HLA == levels(aggregated_data$HLA)[Allele]) %>%
      select(Binding_Core) %>% pull
  }
  
  names(individual_allele_list) <- levels(aggregated_data$HLA)
  
  #create seqlogo plot
  allele_seqlogo <- ggseqlogo(individual_allele_list, ncol = numcol) + 
    ggtitle(status) +
    labs(y = 'Bits', x = 'Position', fill = 'Amino Acid Properties') +
    guides(fill=guide_legend(title="Amino Acid Properties")) +
    theme(strip.text.x = element_text(size=14, face='bold'),
          strip.background.x =  element_blank(),
          axis.text.y = element_text(color='black', size=12),
          axis.title.y = element_text(color='black', size=16),
          axis.text.x = element_text(color='black', size=12, face = 'bold'),
          axis.title.x = element_text(color='black', size=16),
          plot.title = element_text(size=18, hjust=0.5),
          legend.background = element_rect(colour='black', fill='white'),
          legend.title = element_text(size=14, face='bold'),
          legend.text = element_text(size=11))
  
  return(allele_seqlogo)
}


protective_seqlogo <- seqlogo_allele_plotter(Protective_aggregated, 'Protective', 2)
risk_seqlogo <- seqlogo_allele_plotter(Risk_aggregated, 'Risk', 2)

seqlogo_panel <- ggarrange(protective_seqlogo, risk_seqlogo, ncol=2, nrow=1, common.legend = TRUE, legend="bottom") %>% annotate_figure(text_grob("Sequence logos across HLA class II alleles", size = 20, face = 'bold'))

#ggsave(paste0(output_folder, 'Publication_ready_Seqlogo_Panel.png'), seqlogo_panel, width = 22, height = 11)



