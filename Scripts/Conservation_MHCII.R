
#################################################################
##                      Loading in Packages                    ##
#################################################################

dependencies <- c('dplyr', 'readr', 'stringr', 'vroom', 'ggplot2', 'ggsci', 'ggExtra', 'grid', 'gridExtra', 'ggseqlogo', 'here', 'UpSetR')
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
  
  strong_binders <- vector(mode = "list", length = 4)
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
  #' This function filters out any duplicate epitopes, keeping only unique ones
  #' additionally, it filters out epitopes common between the protective and risk data
  #'
  #' @description This function filters out any duplicate epitopes, keeping only unique ones
  #' 
  #' @param protective_data Dataframe containing the data you wish to rid of duplicates
  #' @param risk_data Dataframe containing data you wish to filter out of protective_data
  #' 
  #' @usage duplicate_epitope_filtering(protective_data, risk_data)
  #' @return a dataframe containing filtered data to be used in conservancy across alleles check
  
  between_allele_conservancy_df <- anti_join(protective_data, risk_data, by = 'Peptide') %>% 
    distinct(HLA, Peptide, .keep_all = T)
  
  return(between_allele_conservancy_df)
}

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
  
}

#Total proteins unique to protective
uniq_proteins_LBraziliensis_Protective <- proteins_unique_to_protective(LBraziliensis_Protective, LBraziliensis_Risk)
uniq_proteins_LDonovani_Protective <- proteins_unique_to_protective(LDonovani_Protective, LDonovani_Risk)
uniq_proteins_LInfantum_Protective <- proteins_unique_to_protective(LInfantum_Protective, LInfantum_Risk)
uniq_proteins_LMexicana_Protective <- proteins_unique_to_protective(LMexicana_Protective, LMexicana_Risk)

##################################################################
##                  Conservancy across Alleles                   #
##################################################################

#Check the binding of epitopes across the different alleles associated with the species
Between_Alleles_LBraziliensis_Protective <- duplicate_epitope_filtering(LBraziliensis_Protective, LBraziliensis_Risk) %>% 
  add_count(Peptide) %>% filter(n > 1) %>% select(!HLA) %>%
  inner_join(., LBraziliensis_Protective, by = c('Peptide', 'Binding_Core', 'Identity', 'Rank', 'Affinity', 'Protein')) %>%
  select(8, 1, 2, 3, 4, 5, 6, 7) #%>% filter(., Identity %in% uniq_proteins_LBraziliensis_Protective$Identity)

Between_Alleles_LDonovani_Protective <- duplicate_epitope_filtering(LDonovani_Protective, LDonovani_Risk) %>% 
  add_count(Peptide) %>% filter(n > 1) %>% select(!HLA) %>%
  inner_join(., LDonovani_Protective, by = c('Peptide', 'Binding_Core', 'Identity', 'Rank', 'Affinity', 'Protein')) %>%
  select(8, 1, 2, 3, 4, 5, 6, 7) #%>% filter(., Identity %in% uniq_proteins_LDonovani_Protective$Identity)

Between_Alleles_LInfantum_Protective <- duplicate_epitope_filtering(LInfantum_Protective, LInfantum_Risk) %>% 
  add_count(Peptide) %>% filter(n > 1) %>% select(!HLA) %>%
  inner_join(., LInfantum_Protective, by = c('Peptide', 'Binding_Core', 'Identity', 'Rank', 'Affinity', 'Protein')) %>%
  select(8, 1, 2, 3, 4, 5, 6, 7) #%>% filter(., Identity %in% uniq_proteins_LInfantum_Protective$Identity)

Between_Alleles_LMexicana_Protective <- duplicate_epitope_filtering(LMexicana_Protective, LMexicana_Risk) %>%
  add_count(Peptide) %>% filter(n > 1) %>% select(!HLA) %>%
  inner_join(., LMexicana_Protective, by = c('Peptide', 'Binding_Core', 'Identity', 'Rank', 'Affinity', 'Protein')) %>%
  select(8, 1, 2, 3, 4, 5, 6, 7) #%>% filter(., Identity %in% uniq_proteins_LMexicana_Protective$Identity)

##################################################################
##                  Conservancy across Species                   #
##################################################################

#check the conservation of epitopes across species
#first filter out duplicates and common epitopes between protective and risk
pre_between_Species_LBraziliensis_Protective <- duplicate_epitope_filtering(LBraziliensis_Protective, LBraziliensis_Risk) #%>% filter(., Identity %in% uniq_proteins_LBraziliensis_Protective$Identity)
pre_between_Species_LDonovani_Protective <- duplicate_epitope_filtering(LDonovani_Protective, LDonovani_Risk) #%>% filter(., Identity %in% uniq_proteins_LDonovani_Protective$Identity)
pre_between_Species_LInfantum_Protective <- duplicate_epitope_filtering(LInfantum_Protective, LInfantum_Risk) #%>% filter(., Identity %in% uniq_proteins_LInfantum_Protective$Identity)
pre_between_Species_LMexicana_Protective <- duplicate_epitope_filtering(LMexicana_Protective, LMexicana_Risk) #%>% filter(., Identity %in% uniq_proteins_LMexicana_Protective$Identity)

joining_species_data <- function(species_1, species_2) {
  #' Join the dataframes of two species together by common epitope
  #' 
  #' @description Join the dataframes of two species together by common epitope
  #' 
  #' @param species_1 the first species you wish to join
  #' @param species_2 the second species you wish to join
  #' 
  #' @usage joining_species_data(species_1, species_2)
  #' @return a dataframe containing the two species dataframes joined together by common epitope
  
  joined_species_data <- inner_join(species_1, species_2, by = 'Peptide') %>%
    distinct(Peptide, .keep_all = T) %>% select(2, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
  return(joined_species_data)
  
}

#join species for every combination possible and rename columns
LBraziliensis_with_LDonovani <- joining_species_data(pre_between_Species_LBraziliensis_Protective, pre_between_Species_LDonovani_Protective)
colnames(LBraziliensis_with_LDonovani) <- c("Peptide",
                                            "HLA_LBraziliensis", "Binding_Core_LBrazilienzis", "Identity_LBraziliensis", "Rank_LBraziliensis", "Affinity_LBraziliensis", "Protein_LBraziliensis",
                                            "HLA_LDonovani", "Binding_Core_LDonovani", "Identity_LDonovani", "Rank_LDonovani", "Affinity_LDonovani", "Protein_LDonovani")
LBraziliensis_with_LInfantum <- joining_species_data(pre_between_Species_LBraziliensis_Protective, pre_between_Species_LInfantum_Protective)
colnames(LBraziliensis_with_LInfantum) <- c("Peptide",
                                            "HLA_LBraziliensis", "Binding_Core_LBrazilienzis", "Identity_LBraziliensis", "Rank_LBraziliensis", "Affinity_LBraziliensis", "Protein_LBraziliensis",
                                            "HLA_LInfantum", "Binding_Core_LInfantum", "Identity_LInfantum", "Rank_LInfantum", "Affinity_LInfantum", "Protein_LInfantum")
LBraziliensis_with_LMexicana <- joining_species_data(pre_between_Species_LBraziliensis_Protective, pre_between_Species_LMexicana_Protective)
colnames(LBraziliensis_with_LMexicana) <- c("Peptide",
                                            "HLA_LBraziliensis", "Binding_Core_LBrazilienzis", "Identity_LBraziliensis", "Rank_LBraziliensis", "Affinity_LBraziliensis", "Protein_LBraziliensis",
                                            "HLA_LMexicana", "Binding_Core_LMexicana", "Identity_LMexicana", "Rank_LMexicana", "Affinity_LMexicana", "Protein_LMexicana")

LDonovani_with_LInfantum <- joining_species_data(pre_between_Species_LDonovani_Protective, pre_between_Species_LInfantum_Protective)
colnames(LDonovani_with_LInfantum) <- c("Peptide",
                                        "HLA_LDonovani", "Binding_Core_LDonovani", "Identity_LDonovani", "Rank_LDonovani", "Affinity_LDonovani", "Protein_LDonovani",
                                        "HLA_LInfantum", "Binding_Core_LInfantum", "Identity_LInfantum", "Rank_LInfantum", "Affinity_LInfantum", "Protein_LInfantum")
LDonovani_with_LMexicana <- joining_species_data(pre_between_Species_LDonovani_Protective, pre_between_Species_LMexicana_Protective)
colnames(LDonovani_with_LMexicana) <- c("Peptide",
                                        "HLA_LDonovani", "Binding_Core_LDonovani", "Identity_LDonovani", "Rank_LDonovani", "Affinity_LDonovani", "Protein_LDonovani",
                                        "HLA_LMexicana", "Binding_Core_LMexicana", "Identity_LMexicana", "Rank_LMexicana", "Affinity_LMexicana", "Protein_LMexicana")

LInfantum_with_LMexicana <- joining_species_data(pre_between_Species_LInfantum_Protective, pre_between_Species_LMexicana_Protective)
colnames(LInfantum_with_LMexicana) <- c("Peptide",
                                        "HLA_LInfantum", "Binding_Core_LInfantum", "Identity_LInfantum", "Rank_LInfantum", "Affinity_LInfantum", "Protein_LInfantum",
                                        "HLA_LMexicana", "Binding_Core_LMexicana", "Identity_LMexicana", "Rank_LMexicana", "Affinity_LMexicana", "Protein_LMexicana")

conserved_between_all_species <- inner_join(LBraziliensis_with_LDonovani, LDonovani_with_LInfantum, by = 'Peptide') %>% select(!contains('.x'))
colnames(conserved_between_all_species) <- str_remove(colnames(conserved_between_all_species), '[.]y')
conserved_between_all_species <- inner_join(conserved_between_all_species, LInfantum_with_LMexicana, by = 'Peptide') %>% select(!contains('.x'))
colnames(conserved_between_all_species) <- str_remove(colnames(conserved_between_all_species), '[.]y')

LB_LD_LI <- inner_join(LBraziliensis_with_LDonovani, LDonovani_with_LInfantum, by = 'Peptide') %>% select(!contains('.x'))
colnames(LB_LD_LI) <- str_remove(colnames(LB_LD_LI), '[.]y')

LB_LD_LM <- inner_join(LBraziliensis_with_LDonovani, LDonovani_with_LMexicana, by = 'Peptide') %>% select(!contains('.x'))
colnames(LB_LD_LM) <- str_remove(colnames(LB_LD_LM), '[.]y')

LB_LI_LM <- inner_join(LBraziliensis_with_LInfantum, LInfantum_with_LMexicana, by = 'Peptide') %>% select(!contains('.x'))
colnames(LB_LI_LM) <- str_remove(colnames(LB_LI_LM), '[.]y')

LD_LI_LM <- inner_join(LDonovani_with_LInfantum, LInfantum_with_LMexicana, by = 'Peptide') %>% select(!contains('.x'))
colnames(LD_LI_LM) <- str_remove(colnames(LD_LI_LM), '[.]y')

rm(pre_between_Species_LBraziliensis_Protective, pre_between_Species_LDonovani_Protective, pre_between_Species_LInfantum_Protective, pre_between_Species_LMexicana_Protective)

#################################################################
##              Conserved across species Upset Plot            ##
#################################################################

intersections_conserved_species <- c(
  L.braziliensis = 1697,
  L.donovani = 587,
  L.infantum = 1437,
  L.mexicana = 1567,
  "L.braziliensis&L.donovani" = 68,
  "L.braziliensis&L.infantum" = 170,
  "L.braziliensis&L.mexicana" = 147,
  "L.donovani&L.infantum" = 388,
  "L.donovani&L.mexicana" = 135,
  "L.infantum&L.mexicana" = 396,
  "L.braziliensis&L.donovani&L.infantum" = 66,
  "L.braziliensis&L.donovani&L.mexicana" = 44,
  "L.braziliensis&L.infantum&L.mexicana" = 118,
  "L.donovani&L.infantum&L.mexicana" = 130,
  "L.braziliensis&L.donovani&L.infantum&L.mexicana" = 44
)

#png(file = 'UpsetPlot_BetweenSpecies_Conservation.png', width = 1400, height = 700)
  upset(fromExpression(intersections_conserved_species),
        nintersects = 15, 
        nsets = 4, 
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.6, 0.4),
        number.angles = 0, 
        text.scale = 1.4, 
        point.size = 6, 
        line.size = 1,
        sets.bar.color = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")
  )
  grid.text('MHC-II Epitope conservation between Leishmania spp.', x = 0.65, y=0.97, gp=gpar(fontsize=14))
#dev.off()

##################################################################
##  Unique to Protective & Conserved across species Upset Plot   #
##################################################################

intersections_conserved_species <- c(
  L.braziliensis = 1697,
  L.donovani = 587,
  L.infantum = 1437,
  L.mexicana = 1567,
  "L.braziliensis&L.donovani" = 44,
  "L.braziliensis&L.infantum" = 126,
  "L.braziliensis&L.mexicana" = 116,
  "L.donovani&L.infantum" = 318,
  "L.donovani&L.mexicana" = 85,
  "L.infantum&L.mexicana" = 306,
  "L.braziliensis&L.donovani&L.infantum" = 42,
  "L.braziliensis&L.donovani&L.mexicana" = 17,
  "L.braziliensis&L.infantum&L.mexicana" = 69,
  "L.donovani&L.infantum&L.mexicana" = 81,
  "L.braziliensis&L.donovani&L.infantum&L.mexicana" = 17
)

#png(file = 'UpsetPlot_ProtUniq_BetweenSpecies_Conservation.png', width = 1400, height = 700)
  upset(fromExpression(intersections_conserved_species),
        nintersects = 15, 
        nsets = 4, 
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.6, 0.4),
        number.angles = 0, 
        text.scale = 1.4, 
        point.size = 6, 
        line.size = 1,
        sets.bar.color = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")
  )
  grid.text('MHC-II Epitope conservation between Leishmania spp.', x = 0.65, y=0.97, gp=gpar(fontsize=14))
#dev.off()

#################################################################
##              Binding across alleles Upset Plot               #
#################################################################

intersections_binding_lbraziliensis <- c(
  DRB1_1501 = 526,
  DRB1_1502 = 329,
  DRB1_1601 = 489,
  DRB1_1602 = 353,
  "DRB1_1501&DRB1_1502" = 124,
  "DRB1_1601&DRB1_1602" = 214,
  "DRB1_1501&DRB1_1601" = 8,
  "DRB1_1501&DRB1_1602" = 5,
  "DRB1_1502&DRB1_1601" = 62,
  "DRB1_1502&DRB1_1602" = 78
)

intersections_binding_ldonovani <- c(
  DRB1_0101 = 151,
  DRB1_1602 = 436,
  "DRB1_0101&DRB1_1602" = 73
)

intersections_binding_linfantum <- c(
  DRB1_0101 = 143,
  DRB1_1501 = 540,
  DRB1_1502 = 345,
  DRB1_1602 = 409,
  "DRB1_0101&DRB1_1501" = 0,
  "DRB1_0101&DRB1_1502" = 15,
  "DRB1_0101&DRB1_1602" = 64,
  "DRB1_1501&DRB1_1502" = 115,
  "DRB1_1501&DRB1_1602" = 10,
  "DRB1_1502&DRB1_1602" = 87,
  "DRB1_0101&DRB1_1501&DRB1_1502" = 0,
  "DRB1_0101&DRB1_1501&DRB1_1602" = 0,
  "DRB1_0101&DRB1_1502&DRB1_1602" = 15,
  "DRB1_1501&DRB1_1502&DRB1_1602" = 9
)

intersections_binding_lmexicana <- c(
  DRB1_1501 = 542,
  DRB1_1502 = 396,
  DRB1_1602 = 425,
  HLA_DPA10103_DPB10401 = 204,
  "DRB1_1501&DRB1_1502" = 118,
  "DRB1_1501&DRB1_1602" = 9,
  "DRB1_1502&DRB1_1602" = 110
)

upsetplot_lbraziliensis <- upset(fromExpression(intersections_binding_lbraziliensis),
                                 nintersects = 10,  nsets = 4, order.by = "freq", decreasing = T,
                                 mb.ratio = c(0.6, 0.4), number.angles = 0, line.size = 1,
                                 text.scale = 1.6, point.size = 6, 
                                 sets.bar.color = c("#E64B35FF","#F39B7FFF","#8491B4FF","#4DBBD5FF"))

upsetplot_ldonovani <- upset(fromExpression(intersections_binding_ldonovani),
                             nintersects = 3, nsets = 2, order.by = "freq",  decreasing = T,
                             mb.ratio = c(0.6, 0.4), number.angles = 0, line.size = 1,
                             text.scale = 1.6, point.size = 6, 
                             sets.bar.color = c("#8491B4FF","#00A087FF"))

upsetplot_linfantum <- upset(fromExpression(intersections_binding_linfantum),
                             nintersects = 14,  nsets = 4, order.by = "freq", decreasing = T,
                             mb.ratio = c(0.6, 0.4), number.angles = 0, line.size = 1,
                             text.scale = 1.8, point.size = 6, 
                             sets.bar.color = c("#F39B7FFF","#8491B4FF","#4DBBD5FF","#00A087FF"))

upsetplot_lmexicana <- upset(fromExpression(intersections_binding_lmexicana),
                             nintersects = 7,  nsets = 4, order.by = "freq", decreasing = T,
                             mb.ratio = c(0.6, 0.4), number.angles = 0, line.size = 1,
                             text.scale = 1.8, point.size = 5, 
                             sets.bar.color = c("#F39B7FFF","#4DBBD5FF", "#8491B4FF", "#91D1C2FF"))

plotnames <- c('Binding of L. donovani epitopes across protective MHC-II', 'Binding of L. infantum epitopes across protective MHC-II')
plotnames2 <- c('Binding of L. braziliensis epitopes across protective MHC-II', 'Binding of L. mexicana epitopes across protective MHC-II')

upset_plotter_across_alleles <- function(species_plots_list, name_of_plots, filename, plot_width) {
  #' this function combines the separate plots of different species in one plot and saves it to a file
  #' 
  #' @description combines the separate plots of different species in one plot and saves it to a file
  #' 
  #' @param species_plots_list a list containing the plot objects of the different species
  #' @param name_of_plots a character vector with the title you wish to give each plot
  #' @param filename the filename you wish to save the entire plot in
  #' @param plot_width the width of the plot for optimal resolution
  #' 
  #' @usage upset_plotter_across_alleles(species_plots_list, name_of_plots, filename, plot_width)
  #' @return nothing, saves a file with filename param
  
  
  for (plot_num in 1:length(name_of_plots)) {
    print(species_plots_list[[plot_num]])
    grid.text(name_of_plots[plot_num], x = 0.65, y=0.97, gp=gpar(fontsize=14))
    grid.edit('arrange', name = name_of_plots[plot_num])
    vp <- grid.grab()
    species_plots_list[[plot_num]] <- vp
  }
  
  png(file = filename, width = plot_width, height = 700)
    gridExtra::grid.arrange(grobs = species_plots_list, ncol = 2)
  dev.off()
  
}

upset_plotter_across_alleles(list(upsetplot_ldonovani, upsetplot_linfantum),
                             plotnames,
                             'UpSet_AcrossAlleles_LDonovaniLInfantum.png', 1400)

upset_plotter_across_alleles(list(upsetplot_lbraziliensis, upsetplot_lmexicana),
                             plotnames2,
                             'UpSet_AcrossAlleles_LBraziliensisLMexicana.png', 1800)

