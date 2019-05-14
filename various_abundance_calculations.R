# Script to calculate / format the abundances for different groupings
# Some basic plotting of abundance data

library(dplyr)
library(reshape2)
library(ggplot2)

# install.packages("VennDiagram")
library(VennDiagram)

############################################################
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
patient_pallete_270 <- c("#456c00","#77bbff","#75fb3f","#273300","#f5006f","#ac008b","#125700","#ffef77","#00278e","#3d005e","#d84100","#015686","#01289f","#ff8c4c","#0070b5","#8015cd","#feffd6","#02d8aa","#019cf4","#4f2f00","#bbffe9","#c52900","#1b002c","#a3ac00","#5d9200","#f29fff","#231500","#934cfc","#988a00","#002cbb","#ffeb36","#ffa758","#f1f9ff","#000045","#b4004b","#602900","#390048","#e6c400","#00ce3c","#ff7bd0","#8cff56","#e60051","#b89aff","#00474b","#d5fbff","#ff79c2","#1d0016","#00635d","#ff8e33","#992300","#ff6e91","#ffa081","#534a00","#61002d","#ffe1c1","#8c0072","#00405d","#89ffc6","#607500","#64ff6f","#002e52","#9b97ff","#b1ccff","#02c5cd","#5dffba","#beff45","#00112b","#b8ff74","#7f0027","#0074cd","#005c6f","#3f00af","#dd7900","#cced00","#77ffd6","#ffc5b5","#99ffb1","#01ea72","#f0ff88","#007f63","#abff9d","#391200","#003a74","#114b00","#0a0100","#ff5fff","#ffccfb","#00d6b7","#c7ff93","#1efa5a","#005437","#f6af00","#a60024","#ffb7e6","#ea0043","#c7ffbc","#72ab00","#789300","#585500","#c3ff14","#00f398","#ab4a00","#9b7600","#85e5ff","#006235","#130020","#006825","#ff735c","#007a7f","#02a3a4","#4856ff","#bf52ff","#00edbc","#a31bd6","#009642","#e93bee","#e400ae","#ffbdd2","#00cfc7","#f1ffaa","#009b7a","#dd00c9","#ff697d","#004a14","#ff72ac","#ff3d1f","#fffaa3","#5d0400","#027ba4","#01c774","#002655","#00941f","#0a32d7","#82acff","#ff8af3","#ff4165","#001104","#ffd6f2","#efebff","#aebc00","#3e0030","#c5abff","#00402e","#ff4bae","#0275f1","#be89ff","#ffd899","#00c765","#01b208","#97ffd4","#7e9fff","#00fde1","#0050c9","#ff8eb5","#c800cd","#005173","#ff2b95","#76ff7a","#ea0074","#001d70","#009856","#f100a8","#ba6b00","#0293df","#00462d","#ff6862","#f6ff65","#02bbda","#2c2200","#01a876","#e35a00","#e3000f","#ff819e","#5a0039","#a558ff","#e2ffb2","#784800","#016def","#b400a2","#00143c","#00212d","#403d00","#ff75fe","#975300","#166c00","#260008","#917fff","#ff8d89","#01bf7a","#ffa6bf","#800086","#90a100","#cce4ff","#dad800","#52c900","#46a700","#0c0039","#0b0052","#79009d","#003c85","#bb0034","#59e7ff","#af0064","#64001e","#c0007e","#000897","#bd8400","#2b007f","#318400","#31f1ff","#7c8600","#807300","#ffc072","#6f005f","#770040","#e62c00","#2e0019","#005599","#6535e1","#5b0099","#006bd5","#0142a1","#baaf00","#00ab2d","#ffcc40","#edffec","#ef0031","#153100","#abe9ff","#6bbd00","#e5ff4e","#ffdb43","#ffa5ef","#01c4f3","#ffbd8f","#84d100","#bbff84","#9fcdff","#7b0059","#ffe897","#ff8711","#ffa869","#febaff","#20003a","#94002b","#5387ff","#756dff","#fff494","#a5c1ff","#e0ffcf","#002417","#530076","#ff8459","#ffe4ec","#00b650","#0119b7","#c963ff","#a2ff64","#9c6800","#03b6f8","#00a0c2","#00240b","#6297ff","#bd0010","#fff7af","#7d2d00","#cf7aff","#af5600","#322c00","#500028")
my_colour_pallete_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
# my_colour_pallete_soft_8 <- c("#8b90bc","#76cc5d","#a85bd2","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
my_colour_pallete_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
############################################################

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project")

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Load unfiltered data
unfiltered_data.df <- read.csv(file = "Result_tables/other/project_otu_table_unfiltered.csv",header = T)

# Load filtered/processed abundance data with metadata
otu_data.df <- read.csv("Result_tables/other/OTU_counts_abundances_and_metadata.csv")
genus_data.df <- read.csv("Result_tables/other/genus_counts_abundances_and_metadata.csv")
family_data.df <- read.csv("Result_tables/other/family_counts_abundances_and_metadata.csv")
order_data.df <- read.csv("Result_tables/other/order_counts_abundances_and_metadata.csv")
class_data.df <- read.csv("Result_tables/other/class_counts_abundances_and_metadata.csv")
phylum_data.df <- read.csv("Result_tables/other/phylum_counts_abundances_and_metadata.csv")


# --------------------------------------------------------------------------------------------
### Calculate the (min, max, mean, median, stdev, #samples) abundances of each taxa at each taxa level
generate_taxa_summary <- function(mydata, taxa_column, group_by_columns){
  select_columns <- c(taxa_column, group_by_columns, "Sample", "Read_count", "Read_count_rarified", "Relative_abundance", "Relative_abundance_rarified")
  total_samples <- length(unique(mydata$Sample))
  taxa_group_summary <- 
    mydata %>%
    # dplyr::filter(retained = "yes") %>% # keep only those samples that were retained
    dplyr::select_(.dots = select_columns) %>%
    dplyr::group_by_(.dots = c(taxa_column, group_by_columns)) %>%
    dplyr::mutate(N_samples = n_distinct(Sample)) %>% # number of unique samples/index
    dplyr::group_by_(.dots = c(group_by_columns)) %>%
    dplyr::mutate(N_total_samples_in_group = n_distinct(Sample))  %>%
    dplyr::group_by_(.dots = c(group_by_columns, taxa_column)) %>%
    dplyr::select(-Sample) %>%
    dplyr::summarise(N_samples = max(N_samples),
                     N_total_samples_in_group = max(N_total_samples_in_group),
                     Percent_group_samples = round((max(N_samples) / max(N_total_samples_in_group))*100, 2),
                     Percent_total_samples = round((max(N_samples) / total_samples)*100, 2),
                     Mean_read_count = round(mean(Read_count), 2),
                     Median_read_count = median(Read_count),
                     Min_read_count = min(Read_count),
                     Max_read_count = max(Read_count),
                     Summed_read_count = sum(Read_count),
                     
                     Mean_read_count_rarified = round(mean(Read_count_rarified),2),
                     Median_read_count_rarified = median(Read_count_rarified),
                     Min_read_count_rarified = min(Read_count_rarified),
                     Max_read_count_rarified = max(Read_count_rarified),
                     Summed_read_count_rarified = sum(Read_count_rarified),
                     
                     Mean_relative_abundance = round(mean(Relative_abundance), 5),
                     Median_relative_abundance = round(median(Relative_abundance), 5),
                     Min_relative_abundance = round(min(Relative_abundance),5),
                     Max_relative_abundance = round(max(Relative_abundance),5),
                     Summed_relative_abundance = round(sum(Relative_abundance),5),
                     
                     Mean_relative_abundance_rarified = round(mean(Relative_abundance_rarified), 5),
                     Median_relative_abundance_rarified = round(median(Relative_abundance_rarified), 5),
                     Min_relative_abundance_rarified = round(min(Relative_abundance_rarified), 5),
                     Max_relative_abundance_rarified = round(max(Relative_abundance_rarified), 5),
                     Summed_relative_abundance_rarified = round(sum(Relative_abundance_rarified),5)
    ) %>%
    as.data.frame()
  return(taxa_group_summary)
}

# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# Sample type
phylum_community_summary <- generate_taxa_summary(phylum_data.df, "taxonomy_phylum", group_by_columns = c("Remote_Community"))
class_community_summary <- generate_taxa_summary(class_data.df, "taxonomy_class", group_by_columns = c("Remote_Community"))
order_community_summary <- generate_taxa_summary(order_data.df, "taxonomy_order", group_by_columns = c("Remote_Community"))
family_community_summary <- generate_taxa_summary(family_data.df, "taxonomy_family", group_by_columns = c("Remote_Community"))
genus_community_summary <- generate_taxa_summary(genus_data.df, "taxonomy_genus", group_by_columns = c("Remote_Community"))
otu_community_summary <- generate_taxa_summary(otu_data.df, "OTU.ID", group_by_columns = c("Remote_Community"))
otu_community_summary$taxonomy_species <- as.character(unlist(lapply(as.character(otu_community_summary$OTU.ID), function(x) otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species)))

# Write sample type summaries to file
write.csv(phylum_community_summary,file = "Result_tables/abundance_analysis_tables/phylum_summary.csv", row.names = F, quote = F)
write.csv(class_community_summary,file = "Result_tables/abundance_analysis_tables/class_summary.csv", row.names = F, quote = F)
write.csv(order_community_summary,file = "Result_tables/abundance_analysis_tables/order_summary.csv", row.names = F, quote = F)
write.csv(family_community_summary,file = "Result_tables/abundance_analysis_tables/family_summary.csv", row.names = F, quote = F)
write.csv(genus_community_summary,file = "Result_tables/abundance_analysis_tables/genus_summary.csv", row.names = F, quote = F)
write.csv(otu_community_summary,file = "Result_tables/abundance_analysis_tables/otu_summary.csv", row.names = F, quote = F)


# --------------------------
filter_summary_to_top_n <- function(taxa_summary, grouping_variables, abundance_column = "Summed_read_count_rarified", my_top_n = 10){
  # Get the top N taxa by rarified read count. Requires taxa summary table.
  out <- 
    taxa_summary %>%
    dplyr::group_by_(.dots = c(grouping_variables)) %>%
    dplyr::arrange(dplyr::desc(get(abundance_column))) %>%
    dplyr::top_n(my_top_n, get(abundance_column)) %>% 
    # dplyr::arrange_(.dots = c(grouping_variables),abundance_column) %>%
    dplyr::arrange_(.dots = c(grouping_variables)) %>%
    as.data.frame()
  return(out)
}

# Top 10 taxa
phylum_community_summary_top_10 <- filter_summary_to_top_n(taxa_summary = phylum_community_summary,
                                                           grouping_variables = c("Remote_Community"), 
                                                           abundance_column = "Mean_relative_abundance_rarified",
                                                           my_top_n = 10)
class_community_summary_top_10 <- filter_summary_to_top_n(class_community_summary,
                                                          grouping_variables = c("Remote_Community"), 
                                                          abundance_column = "Mean_relative_abundance_rarified",
                                                          my_top_n = 10)
order_community_summary_top_10 <- filter_summary_to_top_n(order_community_summary,
                                                          grouping_variables = c("Remote_Community"), 
                                                          abundance_column = "Mean_relative_abundance_rarified",
                                                          my_top_n = 10)
family_community_summary_top_10 <- filter_summary_to_top_n(family_community_summary,
                                                           grouping_variables = c("Remote_Community"), 
                                                           abundance_column = "Mean_relative_abundance_rarified",
                                                           my_top_n = 10)
genus_community_summary_top_10 <- filter_summary_to_top_n(genus_community_summary,
                                                          grouping_variables = c("Remote_Community"),
                                                          abundance_column = "Mean_relative_abundance_rarified",
                                                          my_top_n = 10)
otu_community_summary_top_10 <- filter_summary_to_top_n(otu_community_summary,
                                                        grouping_variables = c("Remote_Community"), 
                                                        abundance_column = "Mean_relative_abundance_rarified",
                                                        my_top_n = 10)

write.csv(phylum_community_summary_top_10,file = "Result_tables/abundance_analysis_tables/phylum_community_taxa_summary_top_10.csv", row.names = F, quote = F)
write.csv(class_community_summary_top_10,file = "Result_tables/abundance_analysis_tables/class_community_taxa_summary_top_10.csv", row.names = F, quote = F)
write.csv(order_community_summary_top_10,file = "Result_tables/abundance_analysis_tables/order_community_taxa_summary_top_10.csv", row.names = F, quote = F)
write.csv(family_community_summary_top_10,file = "Result_tables/abundance_analysis_tables/family_community_taxa_summary_top_10.csv", row.names = F, quote = F)
write.csv(genus_community_summary_top_10,file = "Result_tables/abundance_analysis_tables/genus_community_taxa_summary_top_10.csv", row.names = F, quote = F)
write.csv(otu_community_summary_top_10,file = "Result_tables/abundance_analysis_tables/otu_community_taxa_summary_top_10.csv", row.names = F, quote = F)

# ------------------------------------------------------------------------------------------------------------
# Function to generate a more simplified summary for any group(s) (ignore taxa).
# No need to calculate number/percent samples in group.
generate_summary <- function(mydata, group_by_columns){
  select_columns <- c(group_by_columns, "Sample", "Read_count", "Read_count_rarified", "Relative_abundance", "Relative_abundance_rarified")
  total_samples <- length(unique(mydata$Sample))
  group_summary <- 
    mydata %>%
    # dplyr::filter(retained = "yes") %>% # keep only those samples that were retained
    dplyr::select_(.dots = select_columns) %>%
    dplyr::group_by_(.dots = c(group_by_columns)) %>%
    dplyr::mutate(N_samples = n_distinct(Sample)) %>% # number of unique samples/index
    dplyr::select(-Sample) %>%
    dplyr::summarise(N_samples = max(N_samples),
                     Percent_total_samples = round((max(N_samples) / total_samples)*100, 2),
                     Mean_read_count = round(mean(Read_count), 2),
                     Median_read_count = median(Read_count),
                     Min_read_count = min(Read_count),
                     Max_read_count = max(Read_count),
                     Summed_read_count = sum(Read_count),
                     
                     Mean_read_count_rarified = round(mean(Read_count_rarified),2),
                     Median_read_count_rarified = median(Read_count_rarified),
                     Min_read_count_rarified = min(Read_count_rarified),
                     Max_read_count_rarified = max(Read_count_rarified),
                     Summed_read_count_rarified = sum(Read_count_rarified),
                     
                     Mean_relative_abundance = round(mean(Relative_abundance), 5),
                     Median_relative_abundance = round(median(Relative_abundance), 5),
                     Min_relative_abundance = round(min(Relative_abundance),5),
                     Max_relative_abundance = round(max(Relative_abundance),5),
                     Summed_relative_abundance = round(sum(Relative_abundance),5),
                     
                     Mean_relative_abundance_rarified = round(mean(Relative_abundance_rarified), 5),
                     Median_relative_abundance_rarified = round(median(Relative_abundance_rarified), 5),
                     Min_relative_abundance_rarified = round(min(Relative_abundance_rarified), 5),
                     Max_relative_abundance_rarified = round(max(Relative_abundance_rarified), 5),
                     Summed_relative_abundance_rarified = round(sum(Relative_abundance_rarified),5)
    ) %>%
    as.data.frame()
  return(group_summary)
}

# Calculate the read count total for each site for UC, CD and CONTROL samples
community_summary.df <- generate_summary(phylum_data.df, group_by_columns = c("Remote_Community"))

# Save 
write.csv(community_summary.df,file = "Result_tables/abundance_analysis_tables/remote_community_summary.csv", row.names = F, quote = F)

# Get the colour pallete for the sample sites
community_colours.df <- unique(phylum_data.df[,c("Remote_Community","Remote_Community_colour")])
community_colours.df <- setNames( as.character(community_colours.df$Remote_Community_colour), as.character(community_colours.df$Remote_Community))

# Plot the mean read count for each sample site as a bar graph
mean_read_count_bar_graph <- ggplot(community_summary.df, aes(x = factor(Remote_Community), y = Mean_read_count)) +
  geom_bar(position = "dodge", stat= "identity", aes(fill = factor(Remote_Community)), color = "black",size =.2) +
  geom_text(aes(label = paste0("n = ",N_samples)),vjust = -0.6,size = 3) +
  # geom_text(aes(label = count, x = type, y = count), position = position_dodge(width = 0.8), vjust = -0.6
  scale_fill_manual(values = community_colours.df, guide = F) +
  xlab("Remote Community") +
  ylab("Mean read count")
  # scale_y_continuous(breaks = seq(0,4000, 500), limits = c(0,3501))
mean_read_count_bar_graph
# ggsave(plot = mean_read_count_bar_graph, filename ="Result_figures/abundance_analysis_plots/paper_UC_CD_CONTROL/sample_type_mean_read_count_bar_graph_UC_CD_CONTROL.pdf", width = 5, height = 5)

# Plot the summed read count for each sample site as a bar graph
summed_read_count_bar_graph <- ggplot(community_summary.df, aes(x = factor(Remote_Community), y = Summed_read_count_rarified)) +
  geom_bar(position = "dodge", stat= "identity", aes(fill = factor(Remote_Community)), color = "black",size =.2) +
  geom_text(aes(label = paste0("n = ",N_samples)),vjust = -0.6,size = 3) +
  # geom_text(aes(label = count, x = type, y = count), position = position_dodge(width = 0.8), vjust = -0.6
  scale_fill_manual(values = community_colours.df, guide = F) +
  xlab("Remote Community") +
  ylab("Summed read count (rarefied)")
  # scale_y_continuous(breaks = seq(0,4000, 500), limits = c(0,3501))
summed_read_count_bar_graph
# ggsave(plot = summed_read_count_bar_graph, filename ="Result_figures/abundance_analysis_plots/sample_type_summed_read_count_bar_graph_UC_CD_CONTROL.pdf", width = 5, height = 5)

# Plot theMean_relative_abundance_rarified for each sample site as a bar graph
mean_relative_abundaunce_bar_graph <- ggplot(community_summary.df, aes(x = factor(Remote_Community), y = Mean_relative_abundance_rarified)) +
  geom_bar(position = "dodge", stat= "identity", aes(fill = factor(Remote_Community)), color = "black",size =.2) +
  geom_text(aes(label = paste0("n = ",N_samples)),vjust = -0.6,size = 3) +
  # geom_text(aes(label = count, x = type, y = count), position = position_dodge(width = 0.8), vjust = -0.6
  scale_fill_manual(values = community_colours.df, guide = F) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.1)) +
  xlab("Sample type") +
  ylab("Mean relative abundance (rarefied)")
  # common_theme
# scale_y_continuous(breaks = seq(0,4000, 500), limits = c(0,3501))
mean_relative_abundaunce_bar_graph
# ggsave(plot = mean_relative_abundaunce_bar_graph, filename ="Result_figures/abundance_analysis_plots/sample_type_mean_relative_abundance_bar_graph_UC_CD_CONTROL.pdf", width = 5, height = 5)


# ------------------------------------------------------------------------------------------------------------
# Calculate the prevalences etc. of each taxa, not grouped

# phylum_summary <- generate_summary(phylum_data.df,group_by_columns = c("taxonomy_phylum"))
# class_summary <- generate_summary(class_data.df,group_by_columns = c("taxonomy_class"))
# order_summary <- generate_summary(order_data.df,group_by_columns = c("taxonomy_order"))
# family_summary <- generate_summary(family_data.df,group_by_columns = c("taxonomy_family"))
# genus_summary <- generate_summary(genus_data.df,group_by_columns = c("taxonomy_genus"))

# --------------------------------------------------------------------------
### Identify the core taxa for each group and make a venn diagram

# Core taxa should be, for example, >=90% prevalent in the respective group

# ----------------------------------------------------
##  Core taxa within each Community

# Genus level, filtered to high prevalence
genus_community_summary_high_prev <- genus_community_summary[genus_community_summary$Percent_group_samples >= 80,]
core_genus_community_summary_high_prev <- as.character(unique(genus_community_summary_high_prev$taxonomy_genus))

# OTU level, filtered to high prevalence
otu_community_summary_high_prev <- otu_community_summary[otu_community_summary$Percent_group_samples >= 80,]
core_otu_community_summary_high_prev <- as.character(unique(otu_community_summary_high_prev$OTU.ID))

# The not filtered summaries are : otu_sample_type_CD_UC_control_summary and genus_sample_type_CD_UC_control_summary
# ---------------------------------------------------------------
# Make a venn diagram for each pair/triplet of groups

# First for each pair of sample types
colour_list <- unique(genus_data.df[,c("Remote_Community","Remote_Community_colour")])
colour_list <- setNames(as.character(colour_list[,2]), as.character(colour_list[,1]))

generate_sample_type_venn_plots <- function(summary_table,taxonomy_column, output_dir, prefix = ""){
  
  # The combined set of taxa from all groups
  all_taxa  <- as.character(unique(summary_table[,taxonomy_column]))
  
  # Get the combinations for all sample types
  combinations <- combn(as.character(unique(summary_table$Remote_Community)),2)
  shared_results.df <- data.frame()
  for (i in 1:ncol(combinations)){
    grid.newpage()
    # Set the two groups to compare
    group_1 <- combinations[1,i]
    group_2 <- combinations[2,i]
    # Get the taxa for both groups
    group_1_taxa <- as.character(summary_table[summary_table$Remote_Community == group_1,taxonomy_column])
    group_2_taxa <- as.character(summary_table[summary_table$Remote_Community == group_2,taxonomy_column])
    # Number of taxa from combined set in each group
    # same as length(group_1_taxa)
    N_taxa_group_1 <- sum(all_taxa %in% group_1_taxa, na.rm = T)
    N_taxa_group_2 <- sum(all_taxa %in% group_2_taxa, na.rm = T)
    # Number of taxa shared by both groups
    shared_taxa <- intersect(group_1_taxa, group_2_taxa)
    N_shared <- length(shared_taxa)
    for (shared_taxa_entry in shared_taxa){
      if (taxonomy_column == "OTU.ID"){
        tax_string <- as.character(summary_table[summary_table$OTU.ID == shared_taxa_entry,"taxonomy_species"][1])
        shared_results.df <- unique(rbind(shared_results.df, data.frame(group_1, group_2, shared_taxa_entry, tax_string)))
      } else{
        shared_results.df <- unique(rbind(shared_results.df, data.frame(group_1, group_2, shared_taxa_entry)))  
      }
      
    }
    pdf(file = paste0(output_dir, prefix, group_1, "_vs_", group_2, ".pdf"),
        height = 5,
        width = 5)
    # Venn diagram orders the sets by their size. So we need to flip the labels/colours to match
    if (N_taxa_group_2 < N_taxa_group_1){
      draw.pairwise.venn(
        area1 = N_taxa_group_1, 
        area2 = N_taxa_group_2,
        cross.area = N_shared,
        category = c(group_1, group_2),
        fill = as.character(c(colour_list[group_1], colour_list[group_2])),
        euler.d = F,
        scaled = F,
      )
    } else{
      draw.pairwise.venn(
        area1 = N_taxa_group_2, 
        area2 = N_taxa_group_1,
        cross.area = N_shared,
        category = c(group_2,group_1),
        fill = rev(as.character(c(colour_list[group_1], colour_list[group_2]))),
        euler.d = F,
        scaled = F,
      )
    }
    dev.off()
  }
  return(shared_results.df)
}

# Genus level, high prevalence
genus_community_high_prev_shared <- generate_sample_type_venn_plots(genus_community_summary_high_prev, 
                                "taxonomy_genus",
                                output_dir = "Result_figures/abundance_analysis_plots/venn_diagrams/genus/",
                                prefix = "genus_prevalence_80pc_remote_community__")
# OTU level, high prevalence
otu_community_high_prev_shared <- generate_sample_type_venn_plots(otu_community_summary_high_prev, 
                                "OTU.ID",
                                output_dir = "Result_figures/abundance_analysis_plots/venn_diagrams/otu/",
                                prefix = "otu_prevalence_80pc_remote_community__")

# Genus level, all  
genus_community_shared <- generate_sample_type_venn_plots(genus_community_summary, 
                                                          "taxonomy_genus",
                                                          output_dir = "Result_figures/abundance_analysis_plots/venn_diagrams/genus/",
                                                          prefix = "genus_all_remote_community__")

# OTU level, all
otu_community_shared <- generate_sample_type_venn_plots(otu_community_summary, 
                                "OTU.ID",
                                output_dir = "Result_figures/abundance_analysis_plots/venn_diagrams/otu/",
                                prefix = "otu_all_remote_community__")

# Save the tables describing the shared taxa
write.csv(genus_community_high_prev_shared, file = "Result_tables/abundance_analysis_tables/genus_community_shared_taxa_high_prev.csv",quote = F, row.names = F)
write.csv(genus_community_shared, file = "Result_tables/abundance_analysis_tables/genus_community_shared_taxa.csv",quote = F, row.names = F)
write.csv(otu_community_high_prev_shared, file = "Result_tables/abundance_analysis_tables/otu_community_shared_taxa_high_prev.csv",quote = F, row.names = F)
write.csv(otu_community_shared, file = "Result_tables/abundance_analysis_tables/otu_community_shared_taxa.csv",quote = F, row.names = F)





# ------------------------------------------------------------------------------------
# Now we want to do the same for disease groups. Since we only have three disease groups in the current analysis,
# specifically UC, CD and Control, we can make a triplet plot.

# Function assumes the input summary table has been reduced to a single sample type 
# (this could be done in the function)
colour_list <- unique(genus_data.df[,c("DX_Groups","DX_Groups_colour")])
colour_list <- setNames(as.character(colour_list[,2]), as.character(colour_list[,1]))

generate_sample_type_dx_groups_venn_plots <- function(summary_table, taxonomy_column, output_dir, prefix = ""){
  
  for (st in unique(summary_table$Sample_Type)){
    grid.newpage()
    print(st)
    temp <- subset(summary_table, Sample_Type == st)
    all_taxa  <- as.character(unique(temp[,taxonomy_column]))
    group_1 <- "CD"
    group_2 <- "CONTROL"
    group_3 <- "UC"
    
    pdf(file = paste0(output_dir, prefix,st, "__", group_1, "_vs_", group_2, "_vs_", group_3, ".pdf"),
        height = 5,
        width = 5)
    
    group_1_taxa <- as.character(temp[temp$DX_Groups == group_1,taxonomy_column])
    group_2_taxa <- as.character(temp[temp$DX_Groups == group_2,taxonomy_column])
    group_3_taxa <- as.character(temp[temp$DX_Groups == group_3,taxonomy_column])
    
    N_taxa_group_1 <- sum(all_taxa %in% group_1_taxa, na.rm = T) # 33, CD
    N_taxa_group_2 <- sum(all_taxa %in% group_2_taxa, na.rm = T) # 23, CONTROL
    N_taxa_group_3 <- sum(all_taxa %in% group_3_taxa, na.rm = T) # 22, UC
    
    N_taxa_group_1_in_group_2 <- sum(group_1_taxa %in% group_2_taxa)
    N_taxa_group_1_in_group_3 <- sum(group_1_taxa %in% group_3_taxa)
    N_taxa_group_2_in_group_3 <- sum(group_2_taxa %in% group_3_taxa)
    
    shared_taxa <- Reduce(intersect, list(group_1_taxa,group_2_taxa,group_3_taxa))
    N_shared <- length(shared_taxa)
    
    draw.triple.venn(
      area1 = N_taxa_group_1,
      area2 = N_taxa_group_2,
      area3 = N_taxa_group_3,
      n12 = N_taxa_group_1_in_group_2,
      n23 = N_taxa_group_2_in_group_3,
      n13 = N_taxa_group_1_in_group_3,
      n123 = N_shared,
      category = c(group_1,group_2, group_3),
      fill = as.character(c(colour_list[group_1], colour_list[group_2], colour_list[group_3]))
    )
    dev.off()
  }
}
generate_sample_type_dx_groups_venn_plots(genus_sample_type_dx_group_CD_UC_control_summary,
                                          "taxonomy_genus",
                                          output_dir = "Result_figures/abundance_analysis_plots/paper_UC_CD_CONTROL/venn_diagrams/genus_within_sample_type/",
                                          prefix = "genus_all_sample_type_")

generate_sample_type_dx_groups_venn_plots(genus_sample_type_dx_group_CD_UC_control_summary_high_prev,
                                          "taxonomy_genus",
                                          output_dir = "Result_figures/abundance_analysis_plots/paper_UC_CD_CONTROL/venn_diagrams/genus_within_sample_type/",
                                          prefix = "genus_prevalence_50pc_sample_type_")

generate_sample_type_dx_groups_venn_plots(otu_sample_type_dx_group_CD_UC_control_summary,
                                          "OTU.ID",
                                          output_dir = "Result_figures/abundance_analysis_plots/paper_UC_CD_CONTROL/venn_diagrams/otu_within_sample_type/",
                                          prefix = "otu_all_sample_type_")

generate_sample_type_dx_groups_venn_plots(otu_sample_type_dx_group_CD_UC_control_summary_high_prev,
                                          "OTU.ID",
                                          output_dir = "Result_figures/abundance_analysis_plots/paper_UC_CD_CONTROL/venn_diagrams/otu_within_sample_type/",
                                          prefix = "otu_prevalence_50pc_sample_type_")


########################
# testing

# temp <- subset(genus_sample_type_dx_group_CD_UC_control_summary_high_prev, Sample_Type == "DU")
# summary_table <- temp
# taxonomy_column <- "taxonomy_genus"
# 
# all_taxa  <- as.character(unique(summary_table[,taxonomy_column]))
# group_1 <- "CD"
# group_2 <- "CONTROL"
# group_3 <- "UC"
# 
# group_1_taxa <- as.character(summary_table[summary_table$DX_Groups == group_1,taxonomy_column])
# group_2_taxa <- as.character(summary_table[summary_table$DX_Groups == group_2,taxonomy_column])
# group_3_taxa <- as.character(summary_table[summary_table$DX_Groups == group_3,taxonomy_column])
# 
# N_taxa_group_1 <- sum(all_taxa %in% group_1_taxa, na.rm = T) # 33, CD
# N_taxa_group_2 <- sum(all_taxa %in% group_2_taxa, na.rm = T) # 23, CONTROL
# N_taxa_group_3 <- sum(all_taxa %in% group_3_taxa, na.rm = T) # 22, UC
# 
# N_taxa_group_1_in_group_2 <- sum(group_1_taxa %in% group_2_taxa)
# N_taxa_group_1_in_group_3 <- sum(group_1_taxa %in% group_3_taxa)
# N_taxa_group_2_in_group_3 <- sum(group_2_taxa %in% group_3_taxa)
# 
# shared_taxa <- Reduce(intersect, list(group_1_taxa,group_2_taxa,group_3_taxa))
# N_shared <- length(shared_taxa)
# # middle is highest
# # left second highest
# # right lowest
# grid.newpage()
# draw.triple.venn(
#   area1 = N_taxa_group_1,
#   area2 = N_taxa_group_2,
#   area3 = N_taxa_group_3,
#   n12 = N_taxa_group_1_in_group_2,
#   n23 = N_taxa_group_2_in_group_3,
#   n13 = N_taxa_group_1_in_group_3,
#   n123 = N_shared,
#   category = c(group_1,group_2, group_3),
#   fill = as.character(c(colour_list[group_1], colour_list[group_2], colour_list[group_3]))
# )
########################

# --------------------------------------------------------------------------------

# Core each sample type
# core_DNS <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "DNS","taxonomy_genus"])
# core_DU <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "DU","taxonomy_genus"])
# core_G <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "G","taxonomy_genus"])
# core_GNS <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "GNS","taxonomy_genus"])
# core_O <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "O","taxonomy_genus"])
# core_R <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "R","taxonomy_genus"])
# core_RC <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "RC","taxonomy_genus"])
# core_TI <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == "TI","taxonomy_genus"])

# N core taxa in each sample type
# DNS,DU,G,GNS,O,R,RC,TI
# N_core_DNS <- sum(core_taxa_sample_type_filtered_all %in% core_DNS, na.rm = T)
# N_core_DU <- sum(core_taxa_sample_type_filtered_all %in% core_DU, na.rm = T)
# N_core_G <- sum(core_taxa_sample_type_filtered_all %in% core_G, na.rm = T)
# N_core_GNS <- sum(core_taxa_sample_type_filtered_all %in% core_GNS, na.rm = T)
# N_core_O <- sum(core_taxa_sample_type_filtered_all %in% core_O, na.rm = T)
# N_core_R <- sum(core_taxa_sample_type_filtered_all %in% core_R, na.rm = T)
# N_core_RC <- sum(core_taxa_sample_type_filtered_all %in% core_RC, na.rm = T)
# N_core_TI <- sum(core_taxa_sample_type_filtered_all %in% core_TI, na.rm = T)



# shared_core_pairs.df <- data.frame()
# sample_type_combinations <- combn(as.character(unique(genus_abundances_sample_type_filtered$Sample_Type)),2)
# for (i in 1:ncol(sample_type_combinations)){
#   group_1 <- sample_type_combinations[1,i]
#   group_2 <- sample_type_combinations[2,i]
#   core_group_1 <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == group_1,"taxonomy_genus"])
#   core_group_2 <- as.character(genus_abundances_sample_type_filtered[genus_abundances_sample_type_filtered$Sample_Type == group_2,"taxonomy_genus"])
#   N_shared <- length(intersect(core_group_1, core_group_2))
#   shared_core_pairs.df <- rbind(shared_core_pairs.df, data.frame(group_1, group_2, N_shared)) 
# }
# temp <- data.frame(shared_core_pairs.df$group_2, shared_core_pairs.df$group_1, shared_core_pairs.df$N_shared)
# names(temp) <- c("group_1", "group_2", "N_shared")
# shared_core_pairs.df <- bind_rows(shared_core_pairs.df, temp)


# ---------------------------------------------------------------
# Make a 5-way venn diagram, just to test


# grid.newpage()
# draw.quintuple.venn(area1 = N_core_DU,
#                     area2 = N_core_G,
#                     area3 = N_core_O,
#                     area4 = N_core_R,
#                     area5 = N_core_TI,
#                     
#                     n12 = subset(shared_core_pairs.df, group_1 == "DU" & group_2 == "G")$N_shared,
#                     n13 = subset(shared_core_pairs.df, group_1 == "DU" & group_2 == "O")$N_shared,
#                     n14 = subset(shared_core_pairs.df, group_1 == "DU" & group_2 == "R")$N_shared,
#                     n15 = subset(shared_core_pairs.df, group_1 == "DU" & group_2 == "TI")$N_shared,
#                     
#                     n23 = subset(shared_core_pairs.df, group_1 == "G" & group_2 == "O")$N_shared,
#                     n24 = subset(shared_core_pairs.df, group_1 == "G" & group_2 == "R")$N_shared,
#                     n25 = subset(shared_core_pairs.df, group_1 == "G" & group_2 == "TI")$N_shared,
#                     
#                     n34 = subset(shared_core_pairs.df, group_1 == "O" & group_2 == "R")$N_shared,
#                     n35 = subset(shared_core_pairs.df, group_1 == "O" & group_2 == "TI")$N_shared,
#                     
#                     n45 = subset(shared_core_pairs.df, group_1 == "R" & group_2 == "TI")$N_shared,
#                     
#                     n123 = length(Reduce(intersect, list(core_DU,core_G,core_O))),
#                     n124 = length(Reduce(intersect, list(core_DU,core_G,core_R))),
#                     n125 = length(Reduce(intersect, list(core_DU,core_G,core_TI))),
#                     
#                     n134 = length(Reduce(intersect, list(core_DU,core_O,core_R))),
#                     n135 = length(Reduce(intersect, list(core_DU,core_O,core_TI))),
#                     
#                     n145 = length(Reduce(intersect, list(core_DU,core_R,core_TI))),
#                     
#                     n234 = length(Reduce(intersect, list(core_G,core_O,core_R))),
#                     n235 = length(Reduce(intersect, list(core_G,core_O,core_TI))),
#                     
#                     n245 = length(Reduce(intersect, list(core_G,core_R,core_TI))),
#                     
#                     n345 = length(Reduce(intersect, list(core_O,core_R,core_TI))),
#                     
#                     n1234 = length(Reduce(intersect, list(core_DU,core_G,core_O,core_R))),
#                     n1235 = length(Reduce(intersect, list(core_DU,core_G,core_O,core_TI))),
#                     n1245 = length(Reduce(intersect, list(core_DU,core_G,core_R,core_TI))),
#                     n1345 = length(Reduce(intersect, list(core_DU,core_O,core_R,core_TI))),
#                     n2345 = length(Reduce(intersect, list(core_G,core_O,core_R,core_TI))),
#                     n12345 = length(Reduce(intersect, list(core_DU,core_G,core_O,core_R,core_TI))),
#                     
#                     category = c("DU", "G", "O", "R", "TI"),
#                     fill = as.character(c(colour_list["DU"], colour_list["G"], colour_list["O"], colour_list["R"], colour_list["TI"]))
#                     # fill = unique(as.character(genus_data.df$Sample_Type_colour))[1:5]
# )
