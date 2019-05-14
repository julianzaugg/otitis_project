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



# Function that calculates the geometric mean with some error-protection bits. 
# DESeq2 does not appear to work (will throw an error) if every OTU (or genus or genome etc.) 
# contains at least one count of zero in every row of the count data.
# Specifically, the function "dds<-DESeq(dds, betaPrior = FALSE)" will fail
# One way to address this is to use the function below as input to DESeq2 to transform the data.
# Calculate the geometric means prior to estimating the size factors
gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

# Center log ratio transform
clr = function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
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
  
  myplot <- ggplot(internal_metadata, aes(x = factor(get(variable)),
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
discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season","Nose", 
"Otitis_status_OM_6mo","Remote_Community_Otitis_status","OM_6mo_Type_OM","Remote_Community_Season")

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
  full_subset[,myvar] <- factor(full_subset[,myvar])
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

# full_subset_Otitis_status_OM_6mo <- full[complete.cases(full[,c("Otitis_status", "OM_6mo")]),]
# full_subset_Otitis_status_OM_6mo$Otitis_status_OM_6mo <- with(full_subset, paste0(Otitis_status, "_", OM_6mo))
# 
# myplot <- 
#   generate_diversity_boxplot(full_subset_Otitis_status_OM_6mo, 
#                              variable = "Otitis_status_OM_6mo",
#                              fill_pallete = my_colour_pallete_10_distinct,
#                              metric = "Chao1",
#                              variable_colours_available = F) +
#   guides(fill = F, color = F) +
#   ggtitle("Chao1") +
#   scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
# ggsave(filename = paste0("Result_figures/diversity_analysis/Otitis_status_OM_6mo__Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
# 

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

# ------------------------------------------------------------------------------------------------------------------------------
# Calculate the beta-diversity for each variable 
# Centre-log transform the counts first and use a euclidean distance. This should be equivalent or superior to 
# a bray curtis transform/distance used on counts. 
# Normally we would use un-rarefied data, however since our data is ntot truly rarified, this should be ok.
# As far as I can tell in the literature, the euclidean distance between CLR values is an appropriate beta diversity measure
beta_diversity_significances <- data.frame("Variable" = character(),
                                           "P_value" = numeric(),
                                           "R_value" = numeric()
)
for (myvar in discrete_variables){
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  otu_rare_subset.m <- otu_rare.m[,rownames(metadata_subset.df)]
  # print(myvar)
  temp <- with(metadata_subset.df, anosim(t(clr(otu_rare_subset.m)),get(myvar), distance = "euclidean",permutations = 999))
  beta_diversity_significances <- rbind(beta_diversity_significances, data.frame("Variable" = myvar,
                                                                                 "P_value" = temp$signif,
                                                                                 "R_value" = temp$statistic))
}

write.csv(beta_diversity_significances, file = "Result_tables/diversity_analysis/variable_beta_diversity_significance.csv", row.names = F, quote = F)

# Calculate the beta-diversity for each variable within each community 
beta_diversity_significances <- data.frame("Remote_Community" = character(),
                                           "Variable" = character(),
                                           "P_value" = numeric(),
                                           "R_value" = numeric()
)
for (community in unique(metadata.df$Remote_Community)){
  for (myvar in discrete_variables){
    if (myvar == "Remote_Community") {next}
    metadata_subset.df <- subset(metadata.df, Remote_Community == community)
    metadata_subset.df <- metadata_subset.df[!is.na(metadata_subset.df[,myvar]),]
    otu_rare_subset.m <- otu_rare.m[,rownames(metadata_subset.df)]
    # print(all(colnames(otu_rare_subset.m) == rownames(metadata_subset.df)))
    # print(myvar)
    temp <- with(metadata_subset.df, anosim(t(clr(otu_rare_subset.m)),get(myvar), distance = "euclidean",permutations = 999))
    beta_diversity_significances <- rbind(beta_diversity_significances, data.frame("Remote_Community" = community,
                                                                                   "Variable" = myvar,
                                                                                   "P_value" = temp$signif,
                                                                                   "R_value" = temp$statistic))
  }
}
write.csv(beta_diversity_significances, file = "Result_tables/diversity_analysis/within_community_variable_beta_diversity_significance.csv", row.names = F, quote = F)

# with(metadata.df, anosim(t(otu_rare.m),Remote_Community,distance = "bray",permutations = 10000))
# with(metadata.df, anosim(bray.m,Remote_Community,permutations = 10000))
# p.adjust(p, method = p.adjust.methods, n = length(p))
#betad <- betadiver(t(otu_rare.m), "z")
# ------------------------------------------------------------------------------------------------------------------------------


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

# Compare each group in within each variable 
variable_comparison <- data.frame("Variable" = character(),
                                  "Group_1" = character(),
                                  "Group_2" = character(),
                                  "N_group_1" = character(),
                                  "N_group_2" = character(),
                                  "Shannon_p-value" = character(),
                                  "Simpson_p-value" = character(),
                                  "Chao1_p-value" = character())


for (myvar in discrete_variables) {
  group_combinations <- combn(as.character(unique(full[,myvar])), 2)
  for (i in 1:ncol(group_combinations)) {
    group_1 <- group_combinations[1,i]
    group_2 <- group_combinations[2,i]
    if (any(is.na(c(group_1, group_2)))) {next} 
    group_1_meta <- subset(full, get(myvar) == group_1)
    group_2_meta <- subset(full, get(myvar) == group_2)
    
    n_group_1 <- dim(group_1_meta)[1]
    n_group_2 <- dim(group_2_meta)[1]
    
    # Test on the Shannon diversity
    wilcox_shannon_test <- wilcox.test(group_1_meta$Shannon, group_2_meta$Shannon)
    # Test on the Simpson diversity
    wilcox_simpson_test <- wilcox.test(group_1_meta$Simpson, group_2_meta$Simpson)
    # Test on the Chao1 diversity
    wilcox_chao1_test <- wilcox.test(group_1_meta$Chao1, group_2_meta$Chao1)
    
    
    variable_comparison <- rbind(variable_comparison, data.frame(  "Variable" = myvar,
                                                                   "Group_1" = group_1, 
                                                                   "Group_2" = group_2,
                                                                   "N_Group_1" = n_group_1,
                                                                   "N_Group_2" = n_group_2,
                                                                   "Shannon_p-value" = wilcox_shannon_test$p.value,
                                                                   "Simpson_p-value" = wilcox_simpson_test$p.value,
                                                                   "Chao1_p-value" = wilcox_chao1_test$p.value
    ))
  }
}
write.csv(variable_comparison, file = "Result_tables/diversity_analysis/variable_group_comparison_wilcox.csv", row.names = F, quote = F)


# Compare each group in for each variable within each community
variable_comparison <- data.frame("Community" = character(),
                                  "Variable" = character(),
                                  "Group_1" = character(),
                                  "Group_2" = character(),
                                  "N_group_1" = character(),
                                  "N_group_2" = character(),
                                  "Shannon_p-value" = character(),
                                  "Simpson_p-value" = character(),
                                  "Chao1_p-value" = character())

for (community in unique(full$Remote_Community)){
  for (myvar in discrete_variables) {
    if (myvar == "Remote_Community") {next}
    full_community <- subset(full, Remote_Community == community)
    group_combinations <- combn(as.character(unique(full_community[,myvar])), 2)
    for (i in 1:ncol(group_combinations)) {
      group_1 <- group_combinations[1,i]
      group_2 <- group_combinations[2,i]
      if (any(is.na(c(group_1, group_2)))) {next} 
      group_1_meta <- subset(full_community, get(myvar) == group_1)
      group_2_meta <- subset(full_community, get(myvar) == group_2)
      n_group_1 <- dim(group_1_meta)[1]
      n_group_2 <- dim(group_2_meta)[1]
      
      # Test on the Shannon diversity
      wilcox_shannon_test <- wilcox.test(group_1_meta$Shannon, group_2_meta$Shannon)
      # Test on the Simpson diversity
      wilcox_simpson_test <- wilcox.test(group_1_meta$Simpson, group_2_meta$Simpson)
      # Test on the Chao1 diversity
      wilcox_chao1_test <- wilcox.test(group_1_meta$Chao1, group_2_meta$Chao1)
      
      
      variable_comparison <- rbind(variable_comparison, data.frame( "Community" = community,
                                                                    "Variable" = myvar,
                                                                    "Group_1" = group_1, 
                                                                    "Group_2" = group_2, 
                                                                    "N_Group_1" = n_group_1,
                                                                    "N_Group_2" = n_group_2,
                                                                    "Shannon_p-value" = wilcox_shannon_test$p.value,
                                                                    "Simpson_p-value" = wilcox_simpson_test$p.value,
                                                                    "Chao1_p-value" = wilcox_chao1_test$p.value
      ))
    }
  }
}

write.csv(variable_comparison, file = "Result_tables/diversity_analysis/within_community_variable_group_comparison_wilcox.csv", row.names = F, quote = F)



summary(aov(formula = Chao1~Remote_Community +  Nose, data = full))

# lme requires a random effect variable to be specified. Simply make this the sample ID (Index or Sample_No, the latter is the patient ID)
full$Index <- rownames(full)

# Run linear mixed model for Chao1 index
chao1full_model <- lme(fixed = Chao1 ~ Remote_Community,
                       data = full, random = ~ 1 | Index)
summary(chao1full_model)

# Run linear model to test significance of variables across all variables
# lm will apply an anova or t-test where appropriate
summary(lm(Chao1 ~ Remote_Community, data = full))
summary(lm(Chao1 ~ Season, data = full))
summary(lm(Chao1 ~ Nose, data = full))

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------