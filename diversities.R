# Diversity calculations (Shannon, Chao1, Simpson) for each sample.
# Significance tests of each discrete group comparing diversity indices.
# Generate boxplots for discrete 
# Comparison of continuous variables vs diversity indices

library(vegan)
library(reshape2)
library(dplyr)
library(ggplot2)
library(FSA)
#install.packages("FSA")

# source('http://bioconductor.org/biocLite.R')
# biocLite('phyloseq')
library(phyloseq)
library(nlme)


common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
  panel.background = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white", size = 1),
  strip.text = element_text(size = 6),
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

# ******************************************************************************************
# Taken from : https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# move an OTU table from vegan to phyloseq  
# otu_table(PhyloseqObject) <- otu_table(veganOTUobject, taxa_are_rows=TRUE)  
# move sample data from vegan to phyloseq
# sample_data(PhyloseqObject) <- as.data.frame(veganSampleDataObject)


# Assumes metric is a column, e.g Mean_relative_abundance
generate_diversity_boxplot <- function(mymetadata, variable, metric, fill_pallete = my_colour_pallete_206_distinct, variable_colours_available = F){
  internal_metadata <- mymetadata[!is.na(mymetadata[variable]),]
  variable_values <- factor(as.character(unique(internal_metadata[[variable]])))
  if (variable_colours_available == T){
    color_col_name <- paste0(variable, "_colour")
    variable_colours <- setNames(as.character(unique(internal_metadata[[color_col_name]])), as.character(unique(internal_metadata[[variable]])))
  } else{
    variable_colours <- setNames(fill_pallete[1:length(variable_values)], variable_values)  
  }
  
  myplot <- ggplot(internal_metadata, aes(x = get(variable),
                                          y = get(metric), 
                                          fill = factor(get(variable)))) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.75, preserve = "single"), 
                 size = .2, 
                 # width = .5, 
                 linetype= "dashed") + 
    
    geom_boxplot(outlier.shape = 1,
                 outlier.size = .5,
                 position = position_dodge(width = 0.75, preserve = "single"), 
                 aes(ymin=..lower.., 
                     ymax=..upper..,
                     colour = factor(get(variable))), 
                 coef = 0, 
                 size = .2) + 
    
    geom_boxplot(outlier.shape = NA,
                 # outlier.size = .5,
                 position = position_dodge(width = 0.75, preserve = "single"), 
                 aes(ymin=..lower.., 
                     ymax=..upper..), 
                 coef = 0, 
                 size = .2) +
    # geom_point()+
    scale_color_manual(values = variable_colours, name = variable) +
    scale_fill_manual(values = variable_colours, name = variable) +
    # xlab(gsub("_", " ", variable)) +
    xlab(variable) +
    ylab(metric) +
    coord_flip() +
    common_theme +
    theme(plot.title = element_text(hjust = 0.5))
  # theme(
  #   axis.text.x = element_text(angle = 90),
  #   axis.text.y = element_text(size = axis_text_y),
  # )
  return(myplot)
}

# ******************************************************************************************

# --------------------------------------------------------------------------------
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project/")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Load the count matrices
otu_rare.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts_rarefied.csv", row.names =  1))

# Since we likely removed samples from the count matrix
# in the main script, remove them from the metadata.df here
metadata.df <- metadata.df[rownames(metadata.df) %in%colnames(otu_rare.m),]

# Remove samples that are not in the metadata.
otu_rare.m <- otu_rare.m[,colnames(otu_rare.m) %in% rownames(metadata.df)]

# Define the discrete variables
discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season","Nose")

# create phyloseq object
otu_rare_phyloseq <- otu_table(otu_rare.m, taxa_are_rows=TRUE)

# Estimate Chao1 richness
otu_rare_chao1.df <- estimate_richness(otu_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_rare_chao1.df <- otu_rare_chao1.df[rownames(UC_CONTROL_CD_samples_metadata.df),]

# Combine the metadata and the diversity metrics into a single dataframe
full=cbind(metadata.df, otu_rare_chao1.df)

# Retain only complete cases (remove entries with missing data)
# full <- full[complete.cases(full),]

# ------------------------------------------
# Generate plots
# For each Discrete variable

for (myvar in discrete_variables){
  myplot <- generate_diversity_boxplot(full, variable = myvar,fill_pallete = my_colour_pallete_10_distinct,metric = "Chao1",variable_colours_available = F) + 
    guides(fill = F, color = F) + 
    ggtitle("Chao1") +
    scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
    
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
  
  myplot <- generate_diversity_boxplot(full, variable = myvar,fill_pallete = my_colour_pallete_10_distinct,metric = "Shannon",variable_colours_available = F) + 
    guides(fill = F, color = F) + 
    ggtitle("Shannon") +
    scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")

  myplot <- generate_diversity_boxplot(full, variable = myvar,fill_pallete = my_colour_pallete_10_distinct,metric = "Simpson",variable_colours_available = F) +
    guides(fill = F, color = F) +
    ggtitle("Simpson") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
}

full_subset <- subset(full, Remote_Community == 0)

# For each variable in each community
for (community in unique(metadata.df$Remote_Community)){
  for (myvar in discrete_variables){
    if (myvar == "Remote_Community") {next}
    full_subset <- subset(full, Remote_Community == community)
    full_subset[,myvar] <- factor(full_subset[,myvar])
    myplot <- generate_diversity_boxplot(full_subset, variable = myvar,fill_pallete = my_colour_pallete_10_distinct,metric = "Chao1",variable_colours_available = F) + 
      guides(fill = F, color = F) + 
      ggtitle("Chao1") +
      scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
    
    ggsave(filename = paste0("Result_figures/diversity_analysis/community_",community,"__",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
    
    myplot <- generate_diversity_boxplot(full_subset, variable = myvar,fill_pallete = my_colour_pallete_10_distinct,metric = "Shannon",variable_colours_available = F) + 
      guides(fill = F, color = F) + 
      ggtitle("Shannon") +
      scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
    ggsave(filename = paste0("Result_figures/diversity_analysis/community_",community,"__",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
    
    myplot <- generate_diversity_boxplot(full_subset, variable = myvar,fill_pallete = my_colour_pallete_10_distinct,metric = "Simpson",variable_colours_available = F) +
      guides(fill = F, color = F) +
      ggtitle("Simpson") +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
    ggsave(filename = paste0("Result_figures/diversity_analysis/community_",community,"__",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
  }  
}


# ------------------------------------------

# Generate the diversity summary tables for each variable
# For each discrete variable, calculate the diversity index mean, max, min, median and stdev
# and write to file

for (var in discrete_variables) {
  diversity_summary <- full %>% 
    dplyr::group_by_(var) %>%
    dplyr::summarise(
      Shannon_Mean=mean(Shannon), 
      Shannon_Max=max(Shannon), 
      Shannon_Min=min(Shannon), 
      Shannon_Median=median(Shannon), 
      Shannon_Std=sd(Shannon),
      
      Simpson_Mean=mean(Simpson), 
      Simpson_Max=max(Simpson), 
      Simpson_Min=min(Simpson), 
      Simpson_Median=median(Simpson), 
      Simpson_Std=sd(Simpson),
      
      Chao1_Mean=mean(Chao1), 
      Chao1_Max=max(Chao1), 
      Chao1_Min=min(Chao1), 
      Chao1_Median=median(Chao1), 
      Chao1_Std=sd(Chao1),
      
      N_samples=n()
      ) %>% 
    as.data.frame()
  outfilename <- paste0("Result_tables/diversity_analysis/variable_summaries/", var, "_diversity_summary.csv")
  write.csv(x = diversity_summary, outfilename, row.names = F,quote = F)
}

# Repeat but within each community
for (var in discrete_variables) {
  if (var  == "Remote_Community") {next}
  diversity_summary <- full %>% 
    dplyr::group_by_("Remote_Community", var) %>%
    dplyr::summarise(
      Shannon_Mean=mean(Shannon), 
      Shannon_Max=max(Shannon), 
      Shannon_Min=min(Shannon), 
      Shannon_Median=median(Shannon), 
      Shannon_Std=sd(Shannon),
      
      Simpson_Mean=mean(Simpson), 
      Simpson_Max=max(Simpson), 
      Simpson_Min=min(Simpson), 
      Simpson_Median=median(Simpson), 
      Simpson_Std=sd(Simpson),
      
      Chao1_Mean=mean(Chao1), 
      Chao1_Max=max(Chao1), 
      Chao1_Min=min(Chao1), 
      Chao1_Median=median(Chao1), 
      Chao1_Std=sd(Chao1),
      
      N_samples=n()
    ) %>% 
    as.data.frame()
  outfilename <- paste0("Result_tables/diversity_analysis/variable_summaries_within_community/", var, "_diversity_summary_within_community.csv")
  write.csv(x = diversity_summary, outfilename, row.names = F,quote = F)
}


# Now that we have calculated the diversities for each sample, we can test if diversity distributions are significantly different between groups:

# Tests that can be used to compare multiple discrete groups include:
# kruskal-wallis (non-parametric, data does not need to be normal, typically used for more than two groups)
# anova
# t-tests (for two-groups)
# Wilcoxon rank sum test (non-parametric)
# general linear models
# general mix linear models (includes a random effect)

# ------------------------------
# First the Kruskal-Wallis test

# Each community+group (not pairwise, all groups together)
all_variables_by_community_kruskal_comparison <- data.frame("Remote_Community" = character(),
                                                              "Variable" = character(),
                                                              "Shannon_chi-squared" = character(), 
                                                              "Shannon_p-value" = character(),
                                                              "Simpson_chi-squared" = character(), 
                                                              "Simpson_p-value" = character(),
                                                              "Chao1_chi-squared" = character(), 
                                                              "Chao1_p-value" = character()
)

# Each sample_type+group vs sample_type+group
all_variables_group_pairs_by_community_kruskal_comparison <- data.frame("Remote_Community" = character(),
                                                                          "Grouping_variable" = character(),
                                                                          "Group1" = character(), 
                                                                          "Group2" = character(), 
                                                                          "Shannon_chi-squared" = character(), 
                                                                          "Shannon_p-value" = character(),
                                                                          "Simpson_chi-squared" = character(), 
                                                                          "Simpson_p-value" = character(),
                                                                          "Chao1_chi-squared" = character(), 
                                                                          "Chao1_p-value" = character()
                                                                          )

# For each sample type, limit to pairs of groups in each variable
for (rc in unique(full$Remote_Community)){
  # For each discrete variable
  for (variable in discrete_variables) {
    if (variable  == "Remote_Community") {next}
    metadata_subset.df <- subset(full, Remote_Community == rc) # Subset the metadata to only entries in the sample type
    metadata_subset.df <- metadata_subset.df[!is.na(metadata_subset.df[[variable]]),] # remove NA entries
    if (length(unique(full[[variable]])) == 1) {next} # If only one unique group left, skip
    # ....................................
    # All groups together
    # Perform the Kruskal wallace test on the Shannon diversity
    kw_shannon_test <- kruskal.test(Shannon~get(variable), data = metadata_subset.df)
    # Perform the Kruskal wallace test on the Simpson diversity
    kw_simpson_test <- kruskal.test(Simpson~get(variable), data = metadata_subset.df)
    # Perform the Kruskal wallace test on the Chao1 diversity
    kw_chao1_test <- kruskal.test(Chao1~get(variable), data = metadata_subset.df)
    
    # Store the Sample_Type, variable, groups, chi-squared value and p-value
    all_variables_by_community_kruskal_comparison <- rbind(all_variables_by_community_kruskal_comparison,
                                                             data.frame("Remote_Community" = rc,
                                                                        "Variable" = variable,
                                                                        "Shannon_chi-squared" = round(kw_shannon_test$statistic,5),
                                                                        "Shannon_p-value" = round(kw_shannon_test$p.value,5),
                                                                        "Simpson_chi-squared" = round(kw_simpson_test$statistic,5),
                                                                        "Simpson_p-value" = round(kw_simpson_test$p.value,5),
                                                                        "Chao1_chi-squared" = round(kw_chao1_test$statistic,5),
                                                                        "Chao1_p-value" = round(kw_chao1_test$p.value,5)
                                                             ))
    # ....................................
    
    
    group_combinations <- combn(unique(full[[variable]]), 2) # Determine the group combinations
    for (i in 1:ncol(group_combinations)){ 
      group_1 <- as.character(group_combinations[1,i])
      group_2 <- as.character(group_combinations[2,i])
      if (any(is.na(c(group_1, group_2)))){ # If either group is NA, skip
        next
      }

      # Subset the metadata to entries for the two groups
      metadata_subset_by_group.df <- subset(metadata_subset.df, 
                                            (get(variable) %in% c(group_1, group_2)))
      if (length(unique(metadata_subset_by_group.df[[variable]])) <2){
        next
      }
      # Perform the Kruskal wallace test on the Shannon diversity
      kw_shannon_test <- kruskal.test(Shannon~get(variable), data = metadata_subset_by_group.df)
      # Perform the Kruskal wallace test on the Simpson diversity
      kw_simpson_test <- kruskal.test(Simpson~get(variable), data = metadata_subset_by_group.df)
      # Perform the Kruskal wallace test on the Chao1 diversity
      kw_chao1_test <- kruskal.test(Chao1~get(variable), data = metadata_subset_by_group.df)
      
      # Store the Sample_Type, variable, groups, chi-squared value and p-value
      all_variables_group_pairs_by_community_kruskal_comparison <- rbind(all_variables_group_pairs_by_community_kruskal_comparison,
                                                                     data.frame("Remote_Community" = rc,
                                                                                "Grouping_variable" = variable,
                                                                                "Group1" = group_1,
                                                                                "Group2" = group_2,
                                                                                "Shannon_chi-squared" = round(kw_shannon_test$statistic,5),
                                                                                "Shannon_p-value" = round(kw_shannon_test$p.value,5),
                                                                                "Simpson_chi-squared" = round(kw_simpson_test$statistic,5),
                                                                                "Simpson_p-value" = round(kw_simpson_test$p.value,5),
                                                                                "Chao1_chi-squared" = round(kw_chao1_test$statistic,5),
                                                                                "Chao1_p-value" = round(kw_chao1_test$p.value,5)
                                                                     ))
    }
  }
}

# Order by the smallest p-value for each metric
all_variables_by_community_kruskal_comparison <- all_variables_by_community_kruskal_comparison[order(apply(all_variables_by_community_kruskal_comparison[,c("Shannon_p.value","Simpson_p.value", "Chao1_p.value")], 1, min)),]
all_variables_group_pairs_by_community_kruskal_comparison <- all_variables_group_pairs_by_community_kruskal_comparison[order(apply(all_variables_group_pairs_by_community_kruskal_comparison[,c("Shannon_p.value","Simpson_p.value", "Chao1_p.value")], 1, min)),]

# Order by the variable / grouping variable and sample type
all_variables_by_community_kruskal_comparison <- all_variables_by_community_kruskal_comparison[order(all_variables_by_community_kruskal_comparison$Variable, all_variables_by_community_kruskal_comparison$Remote_Community),]
all_variables_group_pairs_by_community_kruskal_comparison <- all_variables_group_pairs_by_community_kruskal_comparison[order(all_variables_group_pairs_by_community_kruskal_comparison$Grouping_variable, all_variables_group_pairs_by_community_kruskal_comparison$Remote_Community),]

# Remove rownames
rownames(all_variables_by_community_kruskal_comparison) <- NULL
rownames(all_variables_group_pairs_by_community_kruskal_comparison) <- NULL

# Save
write.csv(all_variables_by_community_kruskal_comparison, file = "Result_tables/diversity_analysis/Groups_by_community_kruskal.csv", row.names = F, quote = F)
write.csv(all_variables_group_pairs_by_community_kruskal_comparison, file = "Result_tables/diversity_analysis/Groups_pairs_by_community_kruskal.csv", row.names = F, quote = F)

# ----------------------------

# Now Mann–Whitney U / Wilcox test


# Community vs Community
community_comparison <- data.frame("Group_1" = character(),
                                     "Group_2" = character(),
                                     "Shannon_p-value" = character(),
                                     "Simpson_p-value" = character(),
                                     "Chao1_p-value" = character())

community_combinations <- combn(as.character(unique(full$Remote_Community)), 2)
for (i in 1:ncol(community_combinations)) {
  group_1 <- community_combinations[1,i]
  group_2 <- community_combinations[2,i]
  group_1_meta <- subset(full, Remote_Community == group_1)
  group_2_meta <- subset(full, Remote_Community == group_2)
  
  # Test on the Shannon diversity
  wilcox_shannon_test <- wilcox.test(group_1_meta$Shannon, group_2_meta$Shannon)
  # Test on the Simpson diversity
  wilcox_simpson_test <- wilcox.test(group_1_meta$Simpson, group_2_meta$Simpson)
  # Test on the Chao1 diversity
  wilcox_chao1_test <- wilcox.test(group_1_meta$Chao1, group_2_meta$Chao1)
  
  community_comparison <- rbind(community_comparison, data.frame("Group_1" = group_1, 
                                                                     "Group_2" = group_2, 
                                                                     "Shannon_p-value" = wilcox_shannon_test$p.value,
                                                                     "Simpson_p-value" = wilcox_simpson_test$p.value,
                                                                     "Chao1_p-value" = wilcox_chao1_test$p.value
                                                                     ))
}
write.csv(community_comparison, file = "Result_tables/diversity_analysis/community_wilcox.csv", row.names = F, quote = F)



# disease state vs disease state
dx_groups_comparison <- data.frame("Group_1" = character(),
                                     "Group_2" = character(),
                                     "Shannon_p-value" = character(),
                                     "Simpson_p-value" = character(),
                                     "Chao1_p-value" = character())

dx_groups_combinations <- combn(as.character(unique(full$DX_Groups)), 2)
for (i in 1:ncol(dx_groups_combinations)) {
  group_1 <- dx_groups_combinations[1,i]
  group_2 <- dx_groups_combinations[2,i]
  group_1_meta <- subset(full, DX_Groups == group_1)
  group_2_meta <- subset(full, DX_Groups == group_2)
  
  # Test on the Shannon diversity
  wilcox_shannon_test <- wilcox.test(group_1_meta$Shannon, group_2_meta$Shannon)
  # Test on the Simpson diversity
  wilcox_simpson_test <- wilcox.test(group_1_meta$Simpson, group_2_meta$Simpson)
  # Test on the Chao1 diversity
  wilcox_chao1_test <- wilcox.test(group_1_meta$Chao1, group_2_meta$Chao1)
  
  dx_groups_comparison <- rbind(dx_groups_comparison, data.frame("Group_1" = group_1, 
                                                                     "Group_2" = group_2, 
                                                                     "Shannon_p-value" = wilcox_shannon_test$p.value,
                                                                     "Simpson_p-value" = wilcox_simpson_test$p.value,
                                                                     "Chao1_p-value" = wilcox_chao1_test$p.value
  ))
}
write.csv(dx_groups_comparison, file = "Result_tables/diversity_analysis/dx_groups_wilcox.csv", row.names = F, quote = F)

# DX groups vs dx groups within each sample type
sample_type_dx_groups_comparison <- data.frame("Sample_Type" = character(),
                                   "Group_1" = character(),
                                   "Group_2" = character(),
                                   "Shannon_p-value" = character(),
                                   "Simpson_p-value" = character(),
                                   "Chao1_p-value" = character())
for (st in unique(full$Sample_Type)){
  data_subset <- subset(full, Sample_Type == st)
  dx_groups_combinations <- combn(as.character(unique(data_subset$DX_Groups)), 2)
  for (i in 1:ncol(dx_groups_combinations)) {
    group_1 <- dx_groups_combinations[1,i]
    group_2 <- dx_groups_combinations[2,i]
    group_1_meta <- subset(data_subset, DX_Groups == group_1)
    group_2_meta <- subset(data_subset, DX_Groups == group_2)
    
    # Test on the Shannon diversity
    wilcox_shannon_test <- wilcox.test(group_1_meta$Shannon, group_2_meta$Shannon)
    # Test on the Simpson diversity
    wilcox_simpson_test <- wilcox.test(group_1_meta$Simpson, group_2_meta$Simpson)
    # Test on the Chao1 diversity
    wilcox_chao1_test <- wilcox.test(group_1_meta$Chao1, group_2_meta$Chao1)
    
    sample_type_dx_groups_comparison <- rbind(sample_type_dx_groups_comparison, data.frame("Sample_Type" = st,
                                                                   "Group_1" = group_1, 
                                                                   "Group_2" = group_2, 
                                                                   "Shannon_p-value" = wilcox_shannon_test$p.value,
                                                                   "Simpson_p-value" = wilcox_simpson_test$p.value,
                                                                   "Chao1_p-value" = wilcox_chao1_test$p.value
    ))
  }
}
write.csv(sample_type_dx_groups_comparison, file = "Result_tables/diversity_analysis/sample_type_dx_groups_wilcox.csv", row.names = F, quote = F)


# Variable

a <- subset(full, Sample_Type == "DU")
b <- subset(full, Sample_Type == "TI")
a$Chao1

wilcox.test(a$Chao1, b$Chao1)


a <- subset(full, Sample_Type == "DU" & DX_Groups == "UC")
b <- subset(full, Sample_Type == "DU" & DX_Groups == "CONTROL")

temp <- wilcox.test(a$Chao1, b$Chao1)

kruskal.test(Chao1~Sample_Type, data = full)
kruskal.test(Shannon~Sample_Type, data = full)
kruskal.test(Simpson~Sample_Type, data = full)



summary(aov(formula = Chao1~Sample_Type, data = full))
summary(aov(formula = Chao1~DX_Groups, data = full))
# ANOVA on 
summary(aov(formula = Chao1~Sample_Type + DX_Groups + AGE + BMI + Gender, data = full))

# lme requires a random effect variable to be specified. Simply make this the sample ID (Index or Sample_No, the latter is the patient ID)
full$Index <- rownames(full)


# Run linear mixed model for Chao1 index
chao1full_model <- lme(fixed = Chao1 ~ DX_Groups + Sample_Type + Gender + AGE + BMI,
                       data = full, random = ~ 1 | Index) 
summary(chao1full_model)

# Run linear model to test significance of variables across all variables
# lm will apply an anova or t-test where appropriate
lm_sample_type_chao <- lm(Chao1 ~ Sample_Type, data = full)
lm_dx_groups_chao <- lm(Chao1 ~ DX_Groups, data = full)

summary(aov(formula = Chao1~Sample_Type, data = full))
summary(lm_sample_type_chao)
summary(aov(formula = Chao1~DX_Groups, data = full))
summary(lm_dx_groups_chao)

lm_sample_type_simpson <- lm(Simpson ~ Sample_Type, data = full)
lm_dx_groups_simpson <- lm(Simpson ~ DX_Groups, data = full)



summary(lm_sample_type_simpson)
summary(lm_dx_groups_simpson)


# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------











# 
# # For each discrete variable, calculate the Shannon diversity index mean, max, min, median and stdev
# # and write to file
# for (var in discrete_variables) {
#   # metadata.df[c(var, "Shannon_diversity")]
#   diversity_summary <- metadata.df %>% 
#     group_by_(var) %>%
#     summarise(Shannon_Mean=mean(Shannon_diversity), 
#               Shannon_Max=max(Shannon_diversity), 
#               Shannon_Min=min(Shannon_diversity), 
#               Shannon_Median=median(Shannon_diversity), 
#               Shannon_Std=sd(Shannon_diversity),
#               N_samples=n()) %>% 
#     as.data.frame()
#   outfilename <- paste0("Result_tables/diversity_analysis/variable_summaries/", var, "_diversity_summary.csv")
#   write.csv(x = diversity_summary, outfilename, row.names = F,quote = F)
# }
# 
# # For each discrete variable, within each sample type (location), 
# # calculate the diversity index mean, max,min,median and stdev and write to file
# for (var in discrete_variables) {
#   if (var  == "Sample_Type") {next}
#   diversity_summary <- metadata.df %>%
#     group_by_("Sample_Type", var) %>%
#     summarise(Shannon_Mean=mean(Shannon_diversity),
#               Shannon_Max=max(Shannon_diversity),
#               Shannon_Min=min(Shannon_diversity),
#               Shannon_Median=median(Shannon_diversity),
#               Shannon_Std=sd(Shannon_diversity),
#               N_samples=n()) %>% 
#     as.data.frame()
#   outfilename <- paste0("Result_tables/diversity_analysis/variable_summaries_within_sample_type/", var, "_diversity_summary_within_sample_type.csv")
#   write.csv(x = diversity_summary, outfilename, row.names = F,quote = F)
# }
# 
# # --------------------------------------------------------------------------------
# ## Calculate the significance values with kruskal wallis
# 
# # Each group vs group
# all_variable_groups_kruskal_comparison <- NULL
# for (variable in discrete_variables){
#   group_combinations <- combn(unique(metadata.df[[variable]]), 2)
#   for (i in 1:ncol(group_combinations)){
#     group_1 <- as.character(group_combinations[1,i])
#     group_2 <- as.character(group_combinations[2,i])
#     if (any(is.na(c(group_1, group_2)))){
#       next
#     }
#     metadata_subset.df <- subset(metadata.df, get(variable) %in% c(group_1, group_2))
#     kw_shannon_test <- kruskal.test(Shannon_diversity~get(variable), data = metadata_subset.df)
#     all_variable_groups_kruskal_comparison <- rbind(all_variable_groups_kruskal_comparison, 
#                                                     data.frame(variable, 
#                                                                group_1, 
#                                                                group_2, 
#                                                                "shannon", 
#                                                                round(kw_shannon_test$statistic,5),
#                                                                round(kw_shannon_test$p.value,5)
#                                                     ))
#   }
# }
# names(all_variable_groups_kruskal_comparison) <- c("Grouping_variable","Group1", "Group2", "Diversity_index", "chi-squared", "p-value")
# all_variable_groups_kruskal_comparison <- all_variable_groups_kruskal_comparison[order(all_variable_groups_kruskal_comparison$`p-value`),]
# rownames(all_variable_groups_kruskal_comparison) <- NULL
# write.csv(all_variable_groups_kruskal_comparison, file = "Result_tables/diversity_analysis/all_variable_groups_kruskal.csv", row.names = F, quote = F)
# 
# # --------------------------------------------------------------------------------
# 
# # Each sample_type+group vs sample_type+group
# all_variable_groups_by_sample_type_kruskal_comparison <- NULL
# 
# # For each Sample_Type
# for (st in unique(metadata.df$Sample_Type)){
#   # For each discrete variable
#   for (variable in discrete_variables) {
#     if (variable  == "Sample_Type") {next}
#     metadata_subset.df <- subset(metadata.df, Sample_Type == st) # Subset the metadata to only entries in the sample type
#     if (length(unique(metadata_subset.df[[variable]])) == 1) {next} # If only one unique group left, skip
#     group_combinations <- combn(unique(metadata_subset.df[[variable]]), 2) # Determine the group combinations
#     for (i in 1:ncol(group_combinations)){ 
#       group_1 <- as.character(group_combinations[1,i])
#       group_2 <- as.character(group_combinations[2,i])
#       if (any(is.na(c(group_1, group_2)))){ # If either group is NA, skip
#         next
#       }
#       # Subset the metadata to entries for the two groups
#       metadata_subset_by_group.df <- subset(metadata_subset.df, 
#                                    (get(variable) %in% c(group_1, group_2)))
#       # Perform the Kruskal wallace test on the Shannon diversity
#       kw_shannon_test <- kruskal.test(Shannon_diversity~get(variable), data = metadata_subset_by_group.df)
#       # Store the Sample_Type, variable, groups, chi-squared value and p-value
#       all_variable_groups_by_sample_type_kruskal_comparison <- rbind(all_variable_groups_by_sample_type_kruskal_comparison, 
#                                                       data.frame(st,
#                                                                  variable,
#                                                                  group_1, 
#                                                                  group_2, 
#                                                                  "shannon", 
#                                                                  round(kw_shannon_test$statistic,5),
#                                                                  round(kw_shannon_test$p.value,5)
#                                                       ))
#     }
#   }
# }
# # all_variable_groups_by_sample_type_kruskal_comparison
# names(all_variable_groups_by_sample_type_kruskal_comparison) <- c("Sample_Type","Grouping_variable","Group1", "Group2", "Diversity_index", "chi-squared", "p-value")
# all_variable_groups_by_sample_type_kruskal_comparison <- all_variable_groups_by_sample_type_kruskal_comparison[order(all_variable_groups_by_sample_type_kruskal_comparison$`p-value`),]
# rownames(all_variable_groups_by_sample_type_kruskal_comparison) <- NULL
# write.csv(all_variable_groups_by_sample_type_kruskal_comparison, file = "Result_tables/diversity_analysis/all_variable_groups_by_sample_type_kruskal.csv", row.names = F, quote = F)
# 
# 
# # temp <- metadata.df
# # temp <- subset(temp, Sample_Type == "G")
# # temp$Wheat_related_symptoms_and_group
# # group_1 <- "FGID_No"
# # group_2 <- "Control_no"
# # temp$Wheat_related_symptoms_and_group %in% c(group_1, group_2)
# # temp <- subset(temp, Wheat_related_symptoms_and_group %in% c(group_1, group_2))
# # kruskal.test(Shannon_diversity~Wheat_related_symptoms_and_group, data = temp)
# 
# # --------------------------------------------------------------------------------
# 
# # Dunn test group vs group
# multiple_groups_dunn_test <- NULL
# for (variable in discrete_variables){
#   # print(paste0(variable, " : ", n_groups))
#   metadata_subset.df <- metadata.df[!is.na(metadata.df[[variable]]),]
#   n_groups = length(unique(metadata_subset.df[,variable]))
#   # if (any(is.na(metadata_subset.df[,variable]))){
#   #   next
#   # }
#   if (n_groups > 2){
#     
#     results <- dunnTest(Shannon_diversity~get(variable), data = metadata_subset.df, method = "bh",alpha = 0.05)
#     results_groups_sep <- separate(results$res,Comparison, into = c("Group1", "Group2"), sep = " - ")
#     multiple_groups_dunn_test <- rbind(multiple_groups_dunn_test, data.frame("Variable" = variable, results_groups_sep, "Index" = "Shannon", "N_groups" = n_groups))
#   }
# }
# multiple_groups_dunn_test <- multiple_groups_dunn_test[order(multiple_groups_dunn_test$Index, multiple_groups_dunn_test$P.adj),]
# # multiple_groups_dunn_test
# write.csv(multiple_groups_dunn_test, file = "Result_tables/diversity_analysis/all_variable_groups_dunn.csv", row.names = F, quote = F)
# 
# 
# # --------------------------------------------------------------------------------
# 
# # Each sample_type+group vs sample_type_group
# multiple_groups_by_sample_type_dunn_test <- NULL
# for (st in unique(metadata_subset.df$Sample_Type)){
#     for (variable in discrete_variables){
#       if (variable  == "Sample_Type") {next}
#       metadata_subset.df <- metadata.df[!is.na(metadata.df[[variable]]),]
#       metadata_subset.df <- subset(metadata_subset.df, Sample_Type == st)
#       n_groups = length(unique(metadata_subset.df[,variable]))
#       # if (any(is.na(metadata_subset.df[,variable]))){
#       #   next
#       # }
#       if (n_groups > 2){
#         results <- dunnTest(Shannon_diversity~get(variable), data = metadata_subset.df, method = "bh",alpha = 0.05)
#         results_groups_sep <- separate(results$res,Comparison, into = c("Group1", "Group2"), sep = " - ")
#         multiple_groups_by_sample_type_dunn_test <- rbind(multiple_groups_by_sample_type_dunn_test, data.frame("Sample_Type" = st, "Variable" = variable, results_groups_sep, "Index" = "Shannon", "N_groups" = n_groups))
#         }
#       }
#     }
# multiple_groups_by_sample_type_dunn_test <- multiple_groups_by_sample_type_dunn_test[order(multiple_groups_by_sample_type_dunn_test$Index, multiple_groups_by_sample_type_dunn_test$P.adj),]
# # multiple_groups_by_sample_type_dunn_test
# write.csv(multiple_groups_by_sample_type_dunn_test, file = "Result_tables/diversity_analysis/all_variable_groups_by_sample_type_dunn.csv", row.names = F, quote = F)
# # --------------------------------------------------------------------------------
# 
# # Now generate box plots
# # Want to show the diversity for each sample within each group
# 
# # Ensure each variable is factorised
# metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)
# 
# # generate_diversity_boxplot(mymetadata = metadata.df,
#                            # variable = "DX_Groups",
#                            # fill_pallete = my_colour_pallete_15)
# 
# for (var in discrete_variables){
#   myplot <- generate_diversity_boxplot(mymetadata = metadata.df,
#                              variable = var,
#                              fill_pallete = my_colour_pallete_15)
#   out_filename <- paste0("Result_figures/diversity_analysis/", var,"_diversity_boxplot.pdf")
#   ggsave(filename = out_filename, plot = myplot, width = 15, height = 10, units = "cm")
# }
# 
# for (var in discrete_variables){
#   myplot <- generate_diversity_boxplot(mymetadata = metadata.df,
#                                        variable = var,
#                                        fill_pallete = my_colour_pallete_15) + facet_wrap(~Sample_Type)
#   out_filename <- paste0("Result_figures/diversity_analysis/by_sample_type/", var,"_diversity_boxplot_by_sampletype.pdf")
#   ggsave(filename = out_filename, plot = myplot, width = 25, height = 15, units = "cm")
# }
# 
# 
# # generate_diversity_boxplot(mymetadata = subset(metadata.df, Sample_Type == "DNS"),
# #                           variable = "Smoking",
# #                           fill_pallete = my_colour_pallete_15) + facet_wrap(~Sample_Type)
# # 
# # 
# # metadata.df$Sample_Type_wheat <- factor(with(metadata.df, paste0(Sample_Type, "_", Wheat_related_symptoms_and_group)))
# # generate_diversity_boxplot(mymetadata = metadata.df,
# #                            variable = "Wheat_related_symptoms_and_group",
# #                            fill_pallete = my_colour_pallete_10_distinct) + facet_wrap(~Sample_Type)
# # 
# # generate_diversity_boxplot(mymetadata = metadata.df,
# #                            variable = "Positive_for_Methane",
# #                            fill_pallete = my_colour_pallete_10_distinct) + facet_wrap(~Sample_Type)
# 
# # -------------------------------------------------------------------------------
# 
# discrete_variables <- c("DX_Groups","Sample_Type", "Wheat_related_symptoms_and_group","Gender","Smoking", "Acid_Eructation", "Dysphagia", "Fullness", "Early_Satiety",
#                         "Postprandial_Pain", "Epigastric_Pain", "Retrosternal_Discomfort", "Pain_defecation", 
#                         "Difficulty_Defecation", "Constipation", "Loose_stool", "Incontinence", "Urgency_defecation",
#                         "Diarrhea", "Loss_apetite", "Abd_cramps", "Sickeness", "Nausea", "Vomiting", "Bloating", 
#                         "Gas_flatulence", "Belching", "Headache", "Fatigue", "Back_Pain", "Depression", "Sleep", "Anxiety", "H_pylori_gastritis",
#                         "Positive_for_H2","Positive_for_Methane","Final_Diagnosis_H2_or_CH4","PPI_Medications")
# 
# continuous_variables <- c("AGE","DU_bacterial_load","TI_bacterial_load","DNS_bacterial_load","Eos_count",
#                           "Intraepithelial_Lymphocytes_counts","Hydrogen_.H2._Baseline","H2_Peak","Peak_Time_H2",
#                           "Methane_.CH4._Baeline","Methane_Peak","Peak_Time_Methane","NCV","NC_total_score",
#                           "NC_ab_pain","NC_fullness","NC_nausea","weight","height","BMI","SAGIS__Diarr_Tot",
#                           "SAGIS_Consti_Tot","SAGIS_Total")
# 
# # Calculate the correlation scores per variable
# correlation_scores_per_variable <- data.frame()
# for (var in continuous_variables){
#   metadata_subset.df <- metadata.df[!is.na(metadata.df[,var]),]
#   values <- as.numeric(metadata_subset.df[,var])
#   n_samples <- length(values)
#   if (n_samples <= 10){next}
#   pearson_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "pearson")
#   spearman_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "spearman")
#   correlation_scores_per_variable <- rbind(correlation_scores_per_variable, data.frame(Variable = var, 
#                                                                                        pearson = pearson_score,
#                                                                                        spearman = spearman_score,
#                                                                                        Number_of_samples = n_samples))
# }
# correlation_scores_per_variable <- correlation_scores_per_variable[rev(order(abs(correlation_scores_per_variable$pearson))),]
# write.csv(correlation_scores_per_variable, file = "Result_tables/diversity_analysis/Correlation_scores/correlation_scores_continuous_variable.csv", row.names = F, quote = F)
# 
# # Calculate the correlation scores per variable for samples in the same sampletype
# correlation_scores_per_variable_per_sample_type <- data.frame()
# for (st in unique(metadata.df$Sample_Type)){
#   for (var in continuous_variables){
#     metadata_subset.df <- metadata.df[!is.na(metadata.df[,var]),]
#     metadata_subset.df <- subset(metadata_subset.df, Sample_Type == st)
#     values <- as.numeric(metadata_subset.df[,var])
#     n_samples <- length(values)
#     if (n_samples <= 10){next}
#     pearson_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "pearson")
#     spearman_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "spearman")
#     correlation_scores_per_variable_per_sample_type <- rbind(correlation_scores_per_variable_per_sample_type, data.frame(Variable = var, 
#                                                                                          pearson = pearson_score,
#                                                                                          spearman = spearman_score,
#                                                                                          Number_of_samples = n_samples,
#                                                                                          Sample_Type = st))
#   }
# }
# correlation_scores_per_variable_per_sample_type <- correlation_scores_per_variable_per_sample_type[rev(order(abs(correlation_scores_per_variable_per_sample_type$pearson))),]
# write.csv(correlation_scores_per_variable_per_sample_type, file = "Result_tables/diversity_analysis/Correlation_scores/correlation_scores_continuous_variable_per_sampletype.csv", row.names = F, quote = F)
# 
# 
# mydataframe <- data.frame()
# for (dis_var in discrete_variables){ # For each discrete variable
#   for (dis_group in unique(metadata.df[, dis_var])){ # For each group in the discrete variable
#     if (is.na(dis_group)) {next} # If the group is NA, next
#     # metadata_dis_subset.df <- subset(metadata.df, get(dis_var) == dis_group) # Get only entries for this group
#       for (var in continuous_variables){
#         metadata_subset.df <- metadata.df[!is.na(metadata.df[,var]),] # get rows where the continuous value is not NA
#         metadata_subset.df <- subset(metadata_subset.df, get(dis_var) == dis_group) # Get rows that match the sample type and discrete state
#         values <- as.numeric(metadata_subset.df[,var])
#         n_samples <- length(values)
#         if (n_samples <= 10){next}
#         pearson_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "pearson")
#         spearman_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "spearman")
#         if (is.na(pearson_score) &  is.na(spearman_score)) {next}
#         mydataframe <- rbind(mydataframe, data.frame(Variable = var, 
#                                                      pearson = pearson_score,
#                                                      spearman = spearman_score,
#                                                      Number_of_samples = n_samples,
#                                                      Discrete_variable = dis_var,
#                                                      Discrete_group = dis_group))
# }
#     }
#   }
# mydataframe <- mydataframe[rev(order(abs(mydataframe$pearson))),]
# write.csv(mydataframe, file = "Result_tables/diversity_analysis/Correlation_scores/correlation_scores_continuous_variable_per_discrete_variable_group.csv", row.names = F, quote = F)
# 
# 
# 
# # Calculate the correlation scores per continuous variable for samples in the same sampletype and discrete state group
# correlation_scores_per_variable_per_sample_type_per_discrete_var <- data.frame()
# for (dis_var in discrete_variables){ # For each discrete variable
#   for (dis_group in unique(metadata.df[, dis_var])){ # For each group in the discrete variable
#     if (is.na(dis_group)) {next} # If the group is NA, next
#     # metadata_dis_subset.df <- subset(metadata.df, get(dis_var) == dis_group) # Get only entries for this group
#     for (st in unique(metadata.df$Sample_Type)){
#       for (var in continuous_variables){
#         metadata_subset.df <- metadata.df[!is.na(metadata.df[,var]),] # get rows where the continuous value is not NA
#         metadata_subset.df <- subset(metadata_subset.df, (Sample_Type == st) & get(dis_var) == dis_group) # Get rows that match the sample type and discrete state
#         values <- as.numeric(metadata_subset.df[,var])
#         n_samples <- length(values)
#         if (n_samples <= 10){next}
#         pearson_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "pearson")
#         spearman_score <- cor(values, metadata_subset.df$Shannon_diversity, method = "spearman")
#         if (is.na(pearson_score) &  is.na(spearman_score)) {next}
#         correlation_scores_per_variable_per_sample_type_per_discrete_var <- rbind(correlation_scores_per_variable_per_sample_type_per_discrete_var, data.frame(Variable = var, 
#                                                                                                                              pearson = pearson_score,
#                                                                                                                              spearman = spearman_score,
#                                                                                                                              Number_of_samples = n_samples,
#                                                                                                                              Sample_Type = st,
#                                                                                                                              Discrete_variable = dis_var,
#                                                                                                                              Discrete_group = dis_group))
#       }
#     }
#   }
# }
# correlation_scores_per_variable_per_sample_type_per_discrete_var <- correlation_scores_per_variable_per_sample_type_per_discrete_var[rev(order(abs(correlation_scores_per_variable_per_sample_type_per_discrete_var$pearson))),]
# write.csv(correlation_scores_per_variable_per_sample_type_per_discrete_var, file = "Result_tables/diversity_analysis/Correlation_scores/correlation_scores_continuous_variable_per_sampletype_per_discrete_variable_group.csv", row.names = F, quote = F)
# 
# 
# 
# 
# 
# metadata.df[!is.na(metadata.df[,"SAGIS_Consti_Tot"]),][,"SAGIS_Consti_Tot"]
# 
# ggplot(metadata.df, aes(x = Methane_.CH4._Baeline, y = Shannon_diversity, color = DX_Groups)) +
#   geom_point() +
#   facet_wrap(~DX_Groups) +
#   common_theme
# 
# 
# 
# generate_diversity_scatter_plot <- function(mymetadata, 
#                                             continuous_variable, #x axis
#                                             facet_variable = NULL,
#                                             fill_pallete = my_colour_pallete_206_distinct,
#                                             # axis_text_y = 9,
#                                             point_alpha = 0.6, point_stroke = .3,
#                                             fill_variable = NULL # variable to colour by
#                                             ){
#   internal_metadata <- mymetadata
#   internal_metadata <- internal_metadata[!is.na(internal_metadata[[variable]]),]
#   internal_metadata[variable] <- as.numeric(internal_metadata[[variable]])
#   
#   xlimits <- c(min(internal_metadata[continuous_variable],na.rm = T), max(internal_metadata[continuous_variable],na.rm = T)*1.05)
#   ylimits <- c(min(internal_metadata["Shannon_diversity"],na.rm = T), max(internal_metadata["Shannon_diversity"],na.rm = T)*1.05)
#   myplot <- ggplot(internal_metadata, aes(x = as.numeric(get(continuous_variable)), y= Shannon_diversity))
#   
#   if (!is.null(fill_variable)){
#     myplot <- myplot + 
#       geom_point(shape = 1, color = "black", alpha = point_alpha, stroke = point_stroke) +
#       geom_point(shape = 16, aes(colour = get(fill_variable)), alpha = point_alpha, stroke = 0) +
#       geom_smooth(se = F, method = "lm", na.rm = T, aes(color = get(fill_variable)), linetype = "dashed", size = .3) +
#       scale_colour_manual(values = fill_pallete, name = fill_variable)
#       
#   } else{
#     myplot <- myplot + 
#       geom_point(shape = 21, fill = "grey30", color = "black", stroke = point_stroke, alpha = point_alpha) +
#       geom_smooth(se = F, method = "lm", na.rm = T, color = "grey20", linetype = "dashed", size = .3)
#   }
#   myplot <- myplot +
#     scale_y_continuous(limits = ylimits) +
#     scale_x_continuous(limits = xlimits) +
#     ylab("Shannon_diversity") +
#     xlab(continuous_variable) +
#     common_theme 
#   if (!is.null(facet_variable)){
#     myplot <- myplot + facet_wrap(~get(facet_variable), scales = "free_x")
#   }
#   return(myplot)
#   
# }
# 
# generate_diversity_scatter_plot(metadata.df, 
#                                 continuous_variable = "Methane_.CH4._Baeline", 
#                                 fill_variable = "DX_Groups",
#                                 facet_variable = "Sample_Type",
#                                 fill_pallete = my_colour_pallete_20)
# 
# 
# generate_diversity_scatter_plot(metadata.df, 
#                                 continuous_variable = "DU_bacterial_load", 
#                                 fill_variable = NULL,#"DX_Groups",
#                                 facet_variable = NULL,#"Sample_Type",
#                                 fill_pallete = my_colour_pallete_10_distinct)

