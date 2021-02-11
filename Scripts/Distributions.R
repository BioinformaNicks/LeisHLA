#!/usr/bin/env Rscript

#################################################################
##                      Loading in Packages                    ##
#################################################################

dependencies <- c('dplyr', 'readr', 'stringr', 'vroom', 'ggplot2', 'ggsci', 'ggExtra', 'ggseqlogo', 'here', 'scales')
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
suppressMessages(library(here))

main_folder <- paste0(str_remove(here(), 'Scripts'), '/')
output_folder <- paste0(str_remove(here(), 'Scripts'), '/Output/')

# This retrieves RScript commandline arguments from the bash script and saves them as a character vector.
args = commandArgs(trailingOnly=TRUE)

#################################################################
##                        Loading in Data                      ##
#################################################################

prediction_result <- vroom(paste0(main_folder, args[1]), col_names=T, col_types = "fcccdd")
prediction_result <- select(prediction_result, c("HLA", "Affinity"))
human_result <- vroom(paste0(main_folder, args[2]), col_names=T, col_types = "fcccdd")
human_result <- select(human_result, c("HLA", "Affinity"))

#################################################################
##                  Selecting only needed HLA's                ##
#################################################################

#Select only those HLA's from the reference that are present in the data you want to standardize, saves memory
human_result <- filter(human_result, HLA %in% levels(prediction_result$HLA))

#Rename the HLA's to include the species
prediction_result$HLA <- paste(str_split(args[1], "_", simplify = TRUE)[1], prediction_result$HLA)
prediction_result$HLA <- as.factor(prediction_result$HLA)
human_result$HLA <- paste("Human", human_result$HLA)
human_result$HLA <- as.factor(human_result$HLA)

#################################################################
##                        Merging of data                      ##
#################################################################

graph_data <- bind_rows(human_result, prediction_result)
colnames(graph_data) <- c("Allele", "Affinity")
rm(human_result, prediction_result)

##################################################################
##                  Creating distribution graphs                ##
##################################################################

# this is the round_any function from the plyr package;
#to avoid conflicts with dplyr, plyr is not loaded
round_to <- function(value, accuracy, f=round) {
  f(value / accuracy) * accuracy
}

step_width <- round_to(max(graph_data$Affinity), 10000)/10

#h object is created with a geom_histogram and the needed data
h <- ggplot(graph_data, aes(x=Affinity, y = after_stat(density), fill=Allele)) +
        geom_histogram(binwidth = step_width, color = 'black', position=position_dodge(step_width*0.8))

#plotdata is taken from the h object, to later create a geom_bar object
h_plotdata <- layer_data(h, i = 1L) #h_plotdata <- ggplot_build(h)$data[[1]]
rm(h)
h_plotdata$Allele <- as.factor(h_plotdata$group)
levels(h_plotdata$Allele) <- levels(graph_data$Allele)

if (length(levels(graph_data$Allele)) > 2) {
  color_palette <- scale_fill_manual(values = rep(pal_npg("nrc")(length(levels(graph_data$Allele))/2),2))
} else {
  color_palette <- scale_fill_manual(values = rep(pal_npg("nrc")(length(levels(graph_data$Allele))),2))
}

#create a geom_bar with the plotdata from geom_hist
distplot <- ggplot(h_plotdata, aes(x=x, y=y, fill = Allele)) +
  geom_bar(color = 'black', stat = "identity") +
  ggtitle('NetMHCIIpan Affinity Distribution') + 
  labs(y = 'Density', x='Affinity (in nM)') +
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
  scale_x_continuous(labels=scales::comma, breaks=seq(0, round_to(max(graph_data$Affinity), 1000), by = step_width)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0,NA)) +
  color_palette

#save to file
ggsave(paste0(output_folder, args[3]), plot = distplot, width = 14)