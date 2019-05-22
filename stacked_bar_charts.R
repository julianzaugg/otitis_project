# Stacked bar graphs for abundances
# Showing taxa abundances for each group
# Requires cowplot!!


library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(cowplot)

# ------------------------------------------------------------------------------------------
# Various colour palletes
my_colour_pallete <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
my_colour_pallete_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
my_colour_pallete_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
my_colour_pallete_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
my_colour_pallete_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
my_colour_pallete_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
my_colour_pallete_32_distinct <- c("#ea7e00","#ca0074","#d1c69b","#474007","#bb00ad","#9c80ff","#be3300","#542e72","#00b9f5","#09436b","#8b0036","#9ac8e6","#ff1059","#959eff","#154a11","#0290f4","#ff7762","#7dbf00","#ff8194","#834c00","#006e73","#f9bb5d","#d6c943","#017229","#00d3a8","#732427","#36e191","#6a8200","#efb3ea","#3227bb","#ff90e1","#e92a12")
# lesion_pallete_7 <- c("#8558d6","#6ee268","#d247ad","#c9d743","#d7453e","#59a237","#d78f2a")
# patient_pallete_45 <- c("#d64530","#585fb1","#795d97","#9e4773","#3f6921","#71692c","#a2b93c","#d571cc","#9b3e97","#33947a","#98ad66","#448a4e","#869ae0","#5ce7af","#e085a3","#dfdc87","#d19be2","#5cb735","#e38269","#3db6c0","#50b565","#50902c","#a98a2c","#dde84a","#db3d76","#5fe485","#7c8329","#b3e791","#6fe965","#5ebce9","#3c86c1","#2a6a45","#65b688","#6651d1","#af4ed3","#df872f","#56e4db","#737cea","#ac464b","#dd37b5","#995b2b","#daac6f","#92e2be","#a2e24b","#e0be3a")
my_colour_pallete_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
my_colour_pallete_10_soft <- c("#9E788F","#4C5B61","#678D58","#AD5233","#A0A083","#4D456A","#588578","#D0AC4C","#2A7BA0","#931621")

my_colour_pallete_12_soft <-c("#9E788F","#4C5B61","#678D58","#AD5233","#A0A083","#4D456A","#588578","#D0AC4C","#2A7BA0","#931621", "#c75a93", "#7c7731")
  
  


# ------------------------------------------------------------------------------------------

common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
  panel.background = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white", size = 1),
  legend.key=element_blank(),
  legend.direction="vertical",
  legend.background = element_rect(colour ="white", size = .3),
  legend.text.align = 0,
  legend.title = element_text(size=10, face="bold"),
  legend.title.align = 0.5,
  legend.margin = margin(c(2,2,2,2)),
  legend.key.height=unit(.4,"cm"),
  legend.text = element_text(size = 8),
  axis.text = element_text(size = 9, colour = "black"),
  axis.title = element_text(size = 10,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))

# ------------------------------------------------------------------------------------------

# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project")

# Just load the combined metadata and abundance tables for simplicity
# otu_data.df <- read.csv("Result_tables/other/OTU_counts_abundances_and_metadata.csv")
species_data.df <- read.csv("Result_tables/other/specie_counts_abundances_and_metadata.csv")
genus_data.df <- read.csv("Result_tables/other/genus_counts_abundances_and_metadata.csv")
family_data.df <- read.csv("Result_tables/other/family_counts_abundances_and_metadata.csv")
order_data.df <- read.csv("Result_tables/other/order_counts_abundances_and_metadata.csv")
class_data.df <- read.csv("Result_tables/other/class_counts_abundances_and_metadata.csv")
phylum_data.df <- read.csv("Result_tables/other/phylum_counts_abundances_and_metadata.csv")


# ------------------------------------------------------------------------------------------

generate_taxa_summary <- function(mydata, taxa_column, abundance_column= "Relative_abundance_rarified"){
  # Generate a summary table of each taxa for each group. Use this to filter to taxa of interest
  # For each taxa, count the :
  #   number of samples it is in
  #   max, min and summed abundance
  taxa_group_summary <-
    mydata %>% 
    dplyr::select(taxa_column, abundance_column) %>%
    group_by_(taxa_column) %>% # use the _ form of the function to allow variables to be used
    dplyr::mutate(In_N_samples = n()) %>%
    
    dplyr::summarise(In_N_samples = max(In_N_samples),
                     Max_abundance = max(get(abundance_column)),
                     Min_abundance = min(get(abundance_column)),
                     Mean_abundance = mean(get(abundance_column)),
                     Summed_abundance = sum(get(abundance_column))) %>%
    arrange(dplyr::desc(Max_abundance)) %>%
    group_by_(taxa_column) %>%
    as.data.frame()
  return(taxa_group_summary)
}

get_top_taxa <- function(taxa_summary, my_top_n = 10){
  taxa_column <- names(taxa_summary)[1]
  most_abundant_taxa <- taxa_summary %>% arrange(dplyr::desc(Mean_abundance)) %>% head(my_top_n) %>% select_(taxa_column) %>% .[[1]] %>% as.character()
  return(most_abundant_taxa)
}

# Function that takes the full abundance + metadata dataframe and the taxa summary dataframe, determines most abundant taxa and 
# relabels the low(er) abundance taxa to "Other"
relabel_low_abundance_taxa <- function(mydata, taxa_summary, my_top_n = 10){
  internal_data <- mydata
  taxa_column <- names(taxa_summary)[1]
  # Get the top_n taxa by the mean abundance
  most_abundant_taxa <- taxa_summary %>% arrange(dplyr::desc(Mean_abundance)) %>% head(my_top_n) %>% select_(taxa_column) %>% .[[1]] %>% as.character()
  internal_data[,taxa_column] <- as.character(internal_data[,taxa_column])
  internal_data[!internal_data[,taxa_column] %in%  most_abundant_taxa, taxa_column] <- "Other"
  return(internal_data)
}


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

class_data_processed.df <- class_data.df

# Count the number of samples for each sample site and filter (though probably not necessary as GNS/DNS have been filtered)
# samples_each_community <- class_data_processed.df %>%
#   group_by(Remote_Community) %>%
#   dplyr::summarise(Count = n_distinct(Sample)) %>%
#   as.data.frame()
# 
# class_data_processed.df <- class_data_processed.df[class_data_processed.df$Remote_Community %in% as.character(samples_each_type[samples_each_community$Count > 2,"Remote_Community"]),]

class_taxa_summary.df <- generate_taxa_summary(class_data_processed.df,
                                               taxa_column = "taxonomy_class",
                                               abundance_column = "Relative_abundance_rarified")

# Relabel lower abundance taxa to "Other"
class_data_processed.df <- relabel_low_abundance_taxa(class_data_processed.df, class_taxa_summary.df, my_top_n = 12)


# Define the taxa colours to ensure consistency in each sub-plot
top_taxa <- get_top_taxa(class_taxa_summary.df,12)

class_taxa_summary.df[class_taxa_summary.df$taxonomy_class %in% top_taxa,]



# Create a bar plot for each variable, annotated by another variable
create_stacked_bar_charts <- function(mydata, facet_variable, annotate_variable){
  # Assume "Other" is first
  taxa_colours.l <- setNames(c("grey", my_colour_pallete_12_soft), unique(mydata$taxonomy_class))
  barplots.l <- list()
  for (fv in unique(mydata[,facet_variable])){
    data_subset <- subset(mydata, get(facet_variable) == fv)
    # print(mydata)
    # print(fc)
    # print( unique(mydata[,facet_variable]))
    data_subset[,annotate_variable] <- factor(data_subset[,annotate_variable], levels = sort(unique(as.character(data_subset[,annotate_variable]))))
    data_subset <- data_subset[order(data_subset[, annotate_variable]),]
    annotation_data <- unique(data_subset[,c("Sample", annotate_variable, paste0(annotate_variable,"_colour"))])
    all_sample_colours <- annotation_data[,paste0(annotate_variable,"_colour")]
    
    variable_values <- unique(annotation_data[,annotate_variable])
    variable_shapes <- setNames(rep(c(25,24,23,22,21),length(variable_values))[1:length(variable_values)],variable_values)
    all_sample_shapes <- as.numeric(
      lapply(
        as.character(sort(annotation_data[,annotate_variable])),
        function(x) variable_shapes[x][[1]]
      )
    )
    
    myplot <- ggplot(data_subset, aes(x = Sample, y = Relative_abundance_rarified, fill = taxonomy_class)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = taxa_colours.l, name = "Taxonomy", guide = F) +
      
      ggtitle(fv) +
      ylab(NULL) + xlab(NULL) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(0,0,0,0), "null"),
        plot.title = element_text(hjust = 0.5, size =10),
        panel.background=element_blank(),
        panel.border=element_blank(),
        axis.line = element_blank()) +
      
      # Add points to the bottom of each bar as defined by the annotate variable
      annotate(geom='point', x = seq(1,length(unique(data_subset$Sample)), 1), y = -0.03, size=2,stroke= 0, shape = all_sample_shapes, color ="black", fill = all_sample_colours)
    barplots.l[[as.character(fv)]] <- myplot # Store the plot
  }
  
  # Create another plot to make the legend. Since this contains all samples, the legend will have all taxa. 
  all_sample_plot <- ggplot(mydata, aes(x = Sample, y = Relative_abundance_rarified, fill = taxonomy_class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(my_colour_pallete_12_soft, "grey"), name = "Taxonomy")
  
  # Create another plot to make the legend for the lesion point.
  colour_list <- unique(mydata[,c(annotate_variable,paste0(annotate_variable,"_colour"))])
  colour_list <- setNames(as.character(colour_list[,2]), as.character(colour_list[,1]))
  
  all_sample_plot_for_colour <- ggplot(mydata, aes(x = Sample, y = Relative_abundance_rarified, fill = get(annotate_variable))) +
    geom_point(aes(shape = get(annotate_variable)), color = "black", stroke = .1) +
    scale_shape_manual(values = c(25,24,23), name = annotate_variable) +
    scale_fill_manual(values = colour_list, name = annotate_variable)
  
  # Extract the legend
  my_legend_taxa <- cowplot::get_legend(all_sample_plot + theme(legend.position = "right", 
                                                                legend.text = element_text(size = 9),
                                                                legend.title = element_text(size =10, face = "bold"),
                                                                legend.justification = "center",
                                                                plot.margin = unit(c(0, 0, 0, 0), "cm")))
  my_legend_colour <- cowplot::get_legend(all_sample_plot_for_colour + theme(legend.position = "right", 
                                                                             legend.text = element_text(size = 9),
                                                                             legend.title = element_text(size =10, face = "bold"),
                                                                             legend.justification = "center",
                                                                             plot.margin = unit(c(0, 0, 0, 0), "cm")))
  my_legend <- plot_grid(my_legend_taxa, my_legend_colour, ncol = 1,nrow = 4) 
  
  
  # Make a grid of plots with the list of plots for both cohorts
  grid_plot <- cowplot::plot_grid(plotlist = barplots.l,ncol = 2)
  grid_plot <- plot_grid(grid_plot, my_legend, rel_widths = c(1,.4), ncol = 2)
  return(grid_plot)
}

my_grid_plot <- create_stacked_bar_charts(mydata = class_data_processed.df, facet_variable = "Remote_Community", annotate_variable = "Otitis_status")

my_grid_plot <- create_stacked_bar_charts(mydata = class_data_processed.df, facet_variable = "Nose", annotate_variable = "Otitis_status")
ggsave(my_grid_plot, filename = "Result_figures/abundance_analysis_plots/stacked_bar_charts.pdf", height = 30, width =80, units = "cm")

