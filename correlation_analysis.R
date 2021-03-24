# Calculate the correlations in abundances between taxa

library(devtools)
# install_github("zdk123/SpiecEasi")
# library(SpiecEasi)
library(Matrix)
library(igraph)
library(psych)
library(ggraph)
library(tidygraph)

otu_relabeller_function <- function(my_labels){
  taxonomy_strings <- unlist(lapply(my_labels, function(x) {
    as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_genus)
  }
  ))
  unlist(lapply(taxonomy_strings, function(x) {
    phylostring <- unlist(strsplit(x, split = ";"))
    paste(phylostring[3], phylostring[6], sep = ";")
  }))
}

genus_relabeller_function <- function(my_labels){
  unlist(lapply(my_labels, 
                function(x) {
                  phylostring <- unlist(strsplit(x, split = ";"))
                  # paste(phylostring[2],phylostring[3], phylostring[6], sep = ";")
                  # paste(phylostring[3], phylostring[6], sep = ";")
                  paste(phylostring[3], phylostring[5], phylostring[6], sep = ";")
                }))
}

combined_otu_labeller <- function(x){
  # print(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
  first_resolved_taxonomy(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
}

genus_relabeller_network <- function(x){
  # print(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
  gsub(".*__(.*)","\\1", first_resolved_taxonomy(x))
}

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T, row.names = "Sequence_file_ID_clean")
metadata.df <- metadata.df[!is.na(metadata.df$Otitis_Status),]

# Load feature taxonomy map
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)
rownames(otu_taxonomy_map.df) <- otu_taxonomy_map.df$OTU.ID

# Define descrete variables
discrete_variables <- c("Nose","Tympanic_membrane", "Otitis_Status",
                        "Season","Community","Gold_Star",
                        "H.influenzae_culture","M.catarrhalis_culture","S.pneumoniae_culture",
                        "Otitis_Status__Gold_Star", "Tympanic_membrane__Gold_Star",
                        "Community__Season","Community__Gold_Star","Community__Otitis_Status",
                        "H.Influenzae_qPCR", "M.catarrhalis_qPCR", "S.pneumoniae_qPCR",
                        "Corynebacterium_pseudodiphtheriticum","Dolosigranulum_pigrum","N_HRV")
discrete_variables_to_add_with_counts <- c("Community","Gold_Star","Season","Nose")


# Load count matrices
otu.df <- read.csv("Result_tables/count_tables/OTU_counts.csv", header =T)
genus.df <- read.csv("Result_tables/count_tables/Genus_counts.csv", header =T)
otu.m <- df2matrix(otu.df)
genus.m <- df2matrix(genus.df)




# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                           Generate fastspar inputs, only required for group subsets

prepare_input_variables <- c("Community", "Nose", "Otitis_Status", "Community__Gold_Star", "Community__Otitis_Status")

<<<<<<< HEAD
for (variable in prepare_input_variables){
  for (group in as.character(unique(metadata.df[,variable]))){
    sample_list <- as.character(subset(metadata.df, get(variable) == group)$Index)
    if (length(sample_list) < 4){
      print(paste0("Variable ", variable, ", group ", group, " has less than 2 samples, skipping"))
      next()
    }
    otu_subset_data.df <- otu.df[,c("OTU.ID", sample_list)]
    genus_subset_data.df <- genus.df[,c("taxonomy_genus", sample_list)]
    otu_subset_data.df <- otu_subset_data.df[otu_subset_data.df$OTU.ID %in% rownames(df2matrix(otu_subset_data.df)[which(apply(df2matrix(otu_subset_data.df), 1, sum) >= 50),]),]
    genus_subset_data.df <- genus_subset_data.df[genus_subset_data.df$taxonomy_genus %in% rownames(df2matrix(genus_subset_data.df)[which(apply(df2matrix(genus_subset_data.df), 1, sum) >= 50),]),]
    names(otu_subset_data.df)[1] <- "#OTU ID"
    names(genus_subset_data.df)[1] <- "#OTU ID"

    write.table(x = otu_subset_data.df, file = paste0("Result_tables/fastspar_inputs/", variable, "/", variable, "___",gsub("/|\\s", "_",group), "___otu_counts_fastspar.tsv"), sep = "\t", quote = F, row.names = F)
    write.table(x = genus_subset_data.df, file = paste0("Result_tables/fastspar_inputs/", variable, "/", variable, "___",gsub("/|\\s", "_",group), "___genus_counts_fastspar.tsv"), sep = "\t", quote = F, row.names = F)
  }
}
=======
# dim(genus.m[,rownames(metadata.df[metadata.df$Community__Otitis_Status == "Remote__Never OM",])])

# for (variable in prepare_input_variables){
#   for (group in as.character(unique(metadata.df[,variable]))){
#     sample_list <- as.character(subset(metadata.df, get(variable) == group)$Index)
#     if (length(sample_list) < 4){
#       print(paste0("Variable ", variable, ", group ", group, " has less than 2 samples, skipping"))
#       next()
#     }
#     otu_subset_data.df <- otu.df[,c("OTU.ID", sample_list)]
#     genus_subset_data.df <- genus.df[,c("taxonomy_genus", sample_list)]
#     # print(dim(otu_subset_data.df))
#     # print(dim(genus_subset_data.df))
#     otu_subset_data.df <- otu_subset_data.df[otu_subset_data.df$OTU.ID %in% rownames(df2matrix(otu_subset_data.df)[which(apply(df2matrix(otu_subset_data.df), 1, sum) >= 50),]),]
#     genus_subset_data.df <- genus_subset_data.df[genus_subset_data.df$taxonomy_genus %in% rownames(df2matrix(genus_subset_data.df)[which(apply(df2matrix(genus_subset_data.df), 1, sum) >= 50),]),]
#     # print(dim(otu_subset_data.df))
#     # print(dim(genus_subset_data.df))
#     names(otu_subset_data.df)[1] <- "#OTU ID"
#     names(genus_subset_data.df)[1] <- "#OTU ID"
# 
#     write.table(x = otu_subset_data.df, file = paste0("Result_tables/fastspar_inputs/", variable, "/", variable, "___",gsub("/|\\s", "_",group), "___otu_counts_fastspar.tsv"), sep = "\t", quote = F, row.names = F)
#     write.table(x = genus_subset_data.df, file = paste0("Result_tables/fastspar_inputs/", variable, "/", variable, "___",gsub("/|\\s", "_",group), "___genus_counts_fastspar.tsv"), sep = "\t", quote = F, row.names = F)
#   }
# }
>>>>>>> ac8d1af8b0bcc3ccad540a152432b91af0083eb7

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                                                             NETWORKS
<<<<<<< HEAD
# 
# NOTE - FASTSPAR SHOULD HAVE BEEN RUN SEPARATELY USING INPUTS GENERATED ABOVE
=======
>>>>>>> ac8d1af8b0bcc3ccad540a152432b91af0083eb7
otu_cor_files <- list.files("Additional_results/fastspar/")[grepl("___otu___correlation.tsv", list.files("Additional_results/fastspar/"))]

source("code/helper_functions.R")
<<<<<<< HEAD
for (cor_file in otu_cor_files){

  otu_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
  otu_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
  file_name_split <- strsplit(cor_file, split = "___")[[1]]
  variable <- file_name_split[1]
  group <- file_name_split[2]
  print(cor_file)
  print(variable)
  print(dim(otu_fastspar_cor.m))
  otu_correlation_network.l <- generate_correlation_network(cor_matrix = otu_fastspar_cor.m,
                                                            p_matrix = otu_fastspar_pval.m,
                                                            relabeller_function = combined_otu_labeller,
                                                            p_value_threshold = 0.01,
                                                            cor_threshold = 0.6,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            node_label_segment_colour = "purple",
                                                            label_colour = "black",
                                                            label_size = 3,
                                                            plot_height = 10,
                                                            plot_width = 10,
                                                            edge_width_min = .5,
                                                            edge_width_max = 2.5,
                                                            edge_alpha = 1,
                                                            # network_layout = "stress",
                                                            network_layout = "fr",
                                                            # network_layout = "kk",
                                                            # exclude_to_from_df = edges_to_remove.df,
                                                            plot_title = paste0(variable, ": ", group, "; ASV correlation"),
                                                            filename= paste0("Result_figures/correlation_analysis/networks/otu/",variable,"___",group,"___feature_correlation_network.pdf"),
                                                            myseed = 1,
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = "pdf")
}
=======
# for (cor_file in otu_cor_files){
#   
#   otu_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
#                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
#   otu_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
#                                               sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
#   file_name_split <- strsplit(cor_file, split = "___")[[1]]
#   variable <- file_name_split[1]
#   group <- file_name_split[2]
#   print(cor_file)
#   print(variable)
#   print(dim(otu_fastspar_cor.m))
#   otu_correlation_network.l <- generate_correlation_network(cor_matrix = otu_fastspar_cor.m,
#                                                             p_matrix = otu_fastspar_pval.m,
#                                                             relabeller_function = combined_otu_labeller,
#                                                             p_value_threshold = 0.01,
#                                                             cor_threshold = 0.6,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             node_label_segment_colour = "purple",
#                                                             label_colour = "black",
#                                                             label_size = 3,
#                                                             plot_height = 10,
#                                                             plot_width = 10,
#                                                             edge_width_min = .5,
#                                                             edge_width_max = 2.5,
#                                                             edge_alpha = 1,
#                                                             # network_layout = "stress",
#                                                             network_layout = "fr",
#                                                             # network_layout = "kk",
#                                                             # exclude_to_from_df = edges_to_remove.df,
#                                                             plot_title = paste0(variable, ": ", group, "; ASV correlation"),
#                                                             filename= paste0("Result_figures/correlation_analysis/networks/otu/",variable,"___",group,"___feature_correlation_network.pdf"),
#                                                             myseed = 1,
#                                                             edgetype = "link",
#                                                             show_p_label = F,
#                                                             file_type = "pdf")
# }
>>>>>>> ac8d1af8b0bcc3ccad540a152432b91af0083eb7

genus_cor_files <- list.files("Additional_results/fastspar/")[grepl("___genus___correlation.tsv", list.files("Additional_results/fastspar/"))]
genus_cor_files <- grep("Nose|Otitis", genus_cor_files,value =T)
file_type <- "svg"
for (cor_file in genus_cor_files){
  genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
                                               sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
  genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
                                                sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
  file_name_split <- strsplit(cor_file, split = "___")[[1]]
  variable <- file_name_split[1]
  group <- file_name_split[2]
  print(variable)
  genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                              p_matrix = genus_fastspar_pval.m,
                                                              # relabeller_function = genus_relabeller_function,
                                                              relabeller_function = first_resolved_taxonomy,
                                                              # relabeller_function = genus_relabeller_network,
                                                              p_value_threshold = 0.05,
                                                              cor_threshold = 0.5,
                                                              node_size = 4,
                                                              node_colour = "grey20",
                                                              node_fill = "grey20",
                                                              node_label_segment_colour = "purple",
                                                              label_colour = "black",
                                                              label_size = 4,
                                                              plot_height = 10,
                                                              plot_width = 10,
                                                              edge_width_min = .5,
                                                              edge_width_max = 2.5,
                                                              edge_alpha = 1,
                                                              # network_layout = "stress",
                                                              network_layout = "fr",
                                                              # network_layout = "kk",
                                                              # exclude_to_from_df = edges_to_remove.df,
                                                              # plot_title = paste0(variable, ": ", group, "; Genus correlation"),
                                                              filename= paste0("Result_figures/correlation_analysis/networks/genus/",variable,"___",group,"___genus_correlation_network.",file_type),
                                                              myseed = 1,
                                                              edgetype = "link",
                                                              show_p_label = F,
                                                              file_type = file_type)
}

genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Genus_correlation.tsv",
<<<<<<< HEAD
=======
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Genus_pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
file_type <- "pdf"
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                            p_matrix = genus_fastspar_pval.m,
                                                            relabeller_function = first_resolved_taxonomy,
                                                            # relabeller_function = genus_relabeller_network,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.3,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            node_label_segment_colour = "purple",
                                                            label_colour = "black",
                                                            label_size = 4,
                                                            plot_height = 10,
                                                            plot_width = 10,
                                                            edge_width_min = .5,
                                                            edge_width_max = 2.5,
                                                            edge_alpha = 1,
                                                            # network_layout = "stress",
                                                            network_layout = "fr",
                                                            # network_layout = "kk",
                                                            # exclude_to_from_df = edges_to_remove.df,
                                                            # plot_title = paste0(variable, ": ", group, "; Genus correlation"),
                                                            filename= paste0("Result_figures/correlation_analysis/networks/genus/Genus_correlation_network.",file_type),
                                                            myseed = 1,
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = file_type
                                                            )

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                                           vs g__Dolosigranulum
#
# Break down Remote + Otitis status, 
# genus level, 
# abundances for "g__Moraxella", "g__Haemophilus", "g__Corynebacterium","g__Streptococcus" vs "g__Dolosigranulum"
# Include in insets plot_feature_correlations_external() results for g__Dolosigranulum

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

taxa_of_interest <- c("g__Dolosigranulum","g__Moraxella", "g__Haemophilus", "g__Corynebacterium","g__Streptococcus")

# ------------------------------------------------------------------------
# Loop over each result and save to file
genus_cor_files <- list.files("Additional_results/fastspar/")[grepl("___genus___correlation.tsv", list.files("Additional_results/fastspar/"))]
# genus_cor_files <- grep("Community__Otitis_Status.*genus", genus_cor_files,value =T)
genus_cor_files <- grep("^Otitis_Status.*genus", genus_cor_files,value =T)
# genus_cor_files <- grep("^Community___.*genus", genus_cor_files,value =T)

file_type = "svg"
for (cor_file in genus_cor_files){
  genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
                                               sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
  genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
                                                sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
  file_name_split <- strsplit(cor_file, split = "___")[[1]]
  variable <- file_name_split[1]
  group <- file_name_split[2]
  print(variable)
  print(group)
  
  rownames(genus_fastspar_cor.m) <- unlist(lapply(rownames(genus_fastspar_cor.m), first_resolved_taxonomy))
  colnames(genus_fastspar_cor.m) <- unlist(lapply(colnames(genus_fastspar_cor.m), first_resolved_taxonomy))
  rownames(genus_fastspar_pval.m) <- unlist(lapply(rownames(genus_fastspar_pval.m), first_resolved_taxonomy))
  colnames(genus_fastspar_pval.m) <- unlist(lapply(colnames(genus_fastspar_pval.m), first_resolved_taxonomy))
  
  print(dim(genus_fastspar_pval.m))
  plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
                                     feature = "g__Dolosigranulum",
                                     p_value_matrix = genus_fastspar_pval.m,
                                     top_n = 25,
                                     plot_width = 7, plot_height = 7,
                                     include_self = F,
                                     include_title = F,
                                     # filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
                                     filename= paste0("Result_figures/correlation_analysis/by_feature/",variable,"___",group,"___genus_dolosigranulum_correlations_top_25.",file_type),
                                     format = file_type)
  
  genus_fastspar_cor.m <- genus_fastspar_cor.m[taxa_of_interest,taxa_of_interest]
  print(genus_fastspar_cor.m)
  genus_fastspar_pval.m <- genus_fastspar_pval.m[taxa_of_interest,taxa_of_interest]
  
  plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
                                     feature = "g__Dolosigranulum",
                                     p_value_matrix = genus_fastspar_pval.m,
                                     top_n = 10,
                                     plot_width = 4, plot_height = 3,
                                     include_self = F,
                                     include_title = F,
                                     # filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
                                     filename= paste0("Result_figures/correlation_analysis/by_feature/",variable,"___",group,"___genus_dolosigranulum_correlations.",file_type),
                                     format = file_type)
}


genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/Genus_correlation.tsv"),
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/Genus_pvalues.tsv"),
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
rownames(genus_fastspar_cor.m) <- unlist(lapply(rownames(genus_fastspar_cor.m), first_resolved_taxonomy))
colnames(genus_fastspar_cor.m) <- unlist(lapply(colnames(genus_fastspar_cor.m), first_resolved_taxonomy))
rownames(genus_fastspar_pval.m) <- unlist(lapply(rownames(genus_fastspar_pval.m), first_resolved_taxonomy))
colnames(genus_fastspar_pval.m) <- unlist(lapply(colnames(genus_fastspar_pval.m), first_resolved_taxonomy))
# genus_fastspar_cor.m[taxa_of_interest,]
# genus_fastspar_cor.m <- genus_fastspar_cor.m[taxa_of_interest,taxa_of_interest]
# genus_fastspar_pval.m <- genus_fastspar_pval.m[taxa_of_interest,taxa_of_interest]
#
# source("code/helper_functions.R")
plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
                                   feature = "g__Dolosigranulum",
                                   p_value_matrix = genus_fastspar_pval.m,
                                   top_n = 100,
                                   plot_width = 30, plot_height = 25,
                                   include_self = F,
                                   include_title = F,
                                   filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
                                   format = "pdf")

# ------------------------------------------------------------------------

genus_palette <- setNames(c("#a84b54","#a2b432", "#c057c5","#6f64cf","#5aae36"),
                          c("g__Dolosigranulum", "g__Moraxella", "g__Haemophilus","g__Corynebacterium", "g__Streptococcus"))

genus.df <- read.csv("Result_tables/count_tables/Genus_counts.csv", header =T)
# genus.df$taxonomy_genus %in% rownames(genus_fastspar_cor.m)
genus.df <- m2df(clr(df2matrix(genus.df)), "taxonomy_genus")
genus_melt.df <- melt(genus.df,variable.name = "Sample", value.name = "Counts")
genus_melt.df <- genus_melt.df[genus_melt.df$Sample %in% metadata.df$Index,]
genus_melt.df <- left_join(genus_melt.df, metadata.df, by = c("Sample" = "Index"))
genus_melt.df$Label <- unlist(lapply(genus_melt.df$taxonomy_genus, first_resolved_taxonomy))
genus_melt.df <- subset(genus_melt.df, Label %in% taxa_of_interest)

# sample delosi
temp <- genus_melt.df[genus_melt.df$Label == "g__Dolosigranulum",]
temp <- dcast(temp, Sample~Label, value.var = "Counts")

temp2 <- genus_melt.df[genus_melt.df$Label != "g__Dolosigranulum",]
genus_melt.df <- left_join(temp2, temp, by = "Sample")

genus_melt.df$g__Dolosigranulum <- log10(genus_melt.df$g__Dolosigranulum)
genus_melt.df$Counts <- log10(genus_melt.df$Counts)
genus_melt.df$Counts[is.infinite(genus_melt.df$Counts)] <- 0


remote_effusion.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__Effusion")
remote_hxom.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__HxOM")
remote_neverom.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__Never OM")
remote_perforation.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__Perforation")
rural_effusion.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__Effusion")
rural_hxom.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__HxOM")
rural_neverom.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__Never OM")
rural_perforation.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__Perforation")

remote_effusion_plot <- 
  ggplot(remote_effusion.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

remote_hxom_plot <- 
  ggplot(remote_hxom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

remote_neverom_plot <- 
  ggplot(remote_neverom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

remote_perforation_plot <- 
  ggplot(remote_perforation.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

rural_effusion_plot <- 
  ggplot(rural_effusion.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

rural_hxom_plot <- 
  ggplot(rural_hxom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

rural_neverom_plot <- 
  ggplot(rural_neverom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

rural_perforation_plot <- 
  ggplot(rural_perforation.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
  geom_point() +
  xlab("log10(Dolosigranulum CLR transformed counts)") +
  ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
  scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
  facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
  scale_y_continuous(limits = c(-0.05,2.5)) +
  theme_bw() + theme(legend.title.align = 0.5)

base <- "Result_figures/correlation_analysis/"

filetype = "svg"
ggsave(plot = remote_effusion_plot,filename = paste0(base, "remote_effusion_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
ggsave(plot = remote_hxom_plot,filename = paste0(base, "remote_hxom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
ggsave(plot = remote_neverom_plot,filename = paste0(base, "remote_neverom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
ggsave(plot = remote_perforation_plot,filename = paste0(base, "remote_perforation_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
ggsave(plot = rural_effusion_plot,filename = paste0(base, "rural_effusion_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
ggsave(plot = rural_hxom_plot,filename = paste0(base, "rural_hxom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
ggsave(plot = rural_neverom_plot,filename = paste0(base, "rural_neverom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
ggsave(plot = rural_perforation_plot,filename = paste0(base, "rural_perforation_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)



# 
# 
# myplot <- 
#   ggplot(genus_melt.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
#   # ggplot(genus_melt.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
#   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
#   geom_point() +
#   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
#   guides(fill=guide_legend(title="Genus"),
#          colour=guide_legend(title="Genus"),
#          shape=guide_legend(title="Genus")) +
#   xlab("log2(Dolosigranulum CLR transformed counts)") +
#   ylab(expression(paste(log[2]~"(Genus CLR transformed counts)"))) + 
#   scale_shape_manual(values = c(25,24,23,22,21)) +
#   scale_fill_manual(values = genus_palette) +
#   scale_colour_manual(values = genus_palette) +
#   # scale_x_continuous(limits = c(2.5,20), breaks = seq(0,20,2.5)) +
#   # scale_y_continuous(limits = c(0,20), breaks = seq(0,20,2.5)) +
#   scale_x_continuous(limits = c(0.3,1.2)) +
#   scale_y_continuous(limits = c(-0.05,2.5)) +
#   theme_bw() +
#   theme(legend.title.align = 0.5)
# myplot
# ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.pdf",
#        plot = myplot,
#        device = "pdf",
#        height = 15,
#        width = 10)
# 
# ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.svg",
#        plot = myplot,
#        device = "svg",
#        height = 15,
#        width = 10)




# ------------------------------------------------------------------------------





# Create all abundance plots
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header =T)
genus_rel_melt.df <- melt(genus_rel.df,variable.name = "Sample", value.name = "Relative_abundance")
genus_rel_melt.df <- genus_rel_melt.df[genus_rel_melt.df$Sample %in% metadata.df$Index,]
genus_rel_melt.df <- left_join(genus_rel_melt.df, metadata.df, by = c("Sample" = "Index"))
genus_rel_melt.df$Label <- unlist(lapply(genus_rel_melt.df$taxonomy_genus, first_resolved_taxonomy))
genus_rel_melt.df <- subset(genus_rel_melt.df, Label %in% taxa_of_interest)

# sample delosi
temp <- genus_rel_melt.df[genus_rel_melt.df$Label == "g__Dolosigranulum",]
temp <- dcast(temp, Sample~Label, value.var = "Relative_abundance")

temp2 <- genus_rel_melt.df[genus_rel_melt.df$Label != "g__Dolosigranulum",]
genus_rel_melt.df <- left_join(temp2, temp, by = "Sample")


# ------------------------------------------------------------------------------


myplot <- 
  ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
                              fill = Label,
                              shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
  geom_point() +
  facet_wrap(Otitis_Status~Community, scales = "free", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),
         colour=guide_legend(title="Genus"),
         shape=guide_legend(title="Genus")) +
  xlab("Dolosigranulum relative abundance") +
  ylab("Genus relative abundance") + 
  scale_shape_manual(values = c(25,24,23,22,21)) +
  scale_fill_manual(values = genus_palette) +
  scale_colour_manual(values = genus_palette) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  theme_bw() +
  theme(legend.title.align = 0.5)
myplot
ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.pdf",
       plot = myplot,
       device = "pdf",
       height = 15,
       width = 10)

ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.svg",
       plot = myplot,
       device = "svg",
       height = 15,
       width = 10)


myplot <- 
  ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
                                fill = Label,
                                shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
  geom_point() +
  facet_wrap(~Otitis_Status, scales = "free", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),
         colour=guide_legend(title="Genus"),
         shape=guide_legend(title="Genus")) +
  xlab("Dolosigranulum relative abundance") +
  ylab("Genus relative abundance") + 
  scale_shape_manual(values = c(25,24,23,22,21)) +
  scale_fill_manual(values = genus_palette) +
  scale_colour_manual(values = genus_palette) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  theme_bw() +
  theme(legend.title.align = 0.5)
myplot
ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_otitis_status_correlations.pdf",
       plot = myplot,
       device = "pdf",
       height = 8,
       width = 10)

ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_otitis_status_correlations.svg",
       plot = myplot,
       device = "svg",
       height = 8,
       width = 10)


myplot <- 
  ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
                                fill = Label,
                                shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
  geom_point() +
  facet_wrap(~Community, scales = "free", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),
         colour=guide_legend(title="Genus"),
         shape=guide_legend(title="Genus")) +
  xlab("Dolosigranulum relative abundance") +
  ylab("Genus relative abundance") + 
  scale_shape_manual(values = c(25,24,23,22,21)) +
  scale_fill_manual(values = genus_palette) +
  scale_colour_manual(values = genus_palette) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  theme_bw() +
  theme(legend.title.align = 0.5)
myplot
ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_correlations.pdf",
       plot = myplot,
       device = "pdf",
       height = 5,
       width = 10)

ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_correlations.svg",
       plot = myplot,
       device = "svg",
       height = 5,
       width = 10)



myplot <- 
  ggplot(genus_melt.df, aes(x = log10(g__Dolosigranulum), y = log10(Counts), 
                                fill = Label,
                                shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
  geom_point() +
  facet_wrap(~Community, scales = "free", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),
         colour=guide_legend(title="Genus"),
         shape=guide_legend(title="Genus")) +
  xlab("Dolosigranulum relative abundance") +
  ylab("Genus relative abundance") + 
  scale_shape_manual(values = c(25,24,23,22,21)) +
  scale_fill_manual(values = genus_palette) +
  scale_colour_manual(values = genus_palette) +
  # scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  # scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  theme_bw() +
  theme(legend.title.align = 0.5)
myplot


myplot <- 
  ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
                            fill = Label,
                            shape = Label)) +
  geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
  geom_point() +
  # facet_wrap(~Community, scales = "free", ncol = 2) +
  guides(fill=guide_legend(title="Genus"),
         colour=guide_legend(title="Genus"),
         shape=guide_legend(title="Genus")) +
  xlab("Dolosigranulum relative abundance") +
  ylab("Genus relative abundance") + 
  scale_shape_manual(values = c(25,24,23,22,21)) +
  scale_fill_manual(values = genus_palette) +
  scale_colour_manual(values = genus_palette) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  theme_bw() +
  theme(legend.title.align = 0.5)
myplot

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# GENUS, Community___Remote
genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community___Remote___genus___correlation.tsv",
                                           sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community___Remote___genus___pvalues.tsv",
                                            sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))

genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                          p_matrix = genus_fastspar_pval.m,
                                                          relabeller_function = first_resolved_taxonomy,
                                                          
                                                          p_value_threshold = 0.05,
                                                          cor_threshold = 0.4,
                                                          node_size = 4,
                                                          node_colour = "grey20",
                                                          node_fill = "grey20",
                                                          label_colour = "black",
                                                          label_size = 3,
                                                          plot_height = 10,
                                                          plot_width = 10,
                                                          edge_width_min = .5,
                                                          edge_width_max = 2.5,
                                                          network_layout = "fr",
                                                          # exclude_to_from_df = edges_to_remove.df,
                                                          plot_title = "Community: Remote; Genus correlation",
                                                          filename="Result_figures/correlation_analysis/networks/Community___Remote___genus_correlation_network.pdf",
                                                          myseed = 1, 
                                                          edgetype = "link",
                                                          show_p_label = F,
                                                          file_type = "pdf")
genus_correlation_network.l$network_plot

# plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
#                                    feature = "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Weeksellaceae;g__Ornithobacterium",
#                                    p_value_matrix = genus_fastspar_pval.m,
#                                    top_n = 10)

source("code/helper_functions.R")
heatmap(genus_fastspar_cor.m)
plot_corrplot(correlation_matrix = genus_fastspar_cor.m,
              p_value_matrix = genus_fastspar_pval.m,
              p_value_threshold = .05,
              relabeller_function = first_resolved_taxonomy,
              label_size = .5,
              make_insig_na = F,
              order = "original",
              file_type = "pdf",
              # filename = "out.pdf",
              plot_height = 10,
              plot_width = 10)

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# GENUS, Community___Rural
genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community___Rural___genus___correlation.tsv",
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community___Rural___genus___pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))

# dim(subset(metadata.df, Community == "Remote"))
# dim(subset(metadata.df, Community == "Rural"))
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                            p_matrix = genus_fastspar_pval.m,
                                                            relabeller_function = first_resolved_taxonomy,
                                                            
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.4,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            label_colour = "black",
                                                            label_size = 3,
                                                            plot_height = 10,
                                                            plot_width = 10,
                                                            edge_width_min = .5,
                                                            edge_width_max = 2.5,
                                                            network_layout = "fr",
                                                            # exclude_to_from_df = edges_to_remove.df,
                                                            plot_title = "Community: Rural; Genus correlation",
                                                            filename="Result_figures/correlation_analysis/networks/Community___Rural___genus_correlation_network.pdf",
                                                            myseed = 1, 
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = "pdf")

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# GENUS, Community__Gold_Star___Remote__Healthy
genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Healthy___genus___correlation.tsv",
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Healthy___genus___pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))

# dim(subset(metadata.df, Community == "Remote"))
# dim(subset(metadata.df, Community == "Rural"))
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                            p_matrix = genus_fastspar_pval.m,
                                                            relabeller_function = first_resolved_taxonomy,
                                                            
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.4,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            label_colour = "black",
                                                            label_size = 3,
                                                            plot_height = 10,
                                                            plot_width = 10,
                                                            edge_width_min = .5,
                                                            edge_width_max = 2.5,
                                                            network_layout = "fr",
                                                            # exclude_to_from_df = edges_to_remove.df,
                                                            plot_title = "Community__Gold_Star: Remote__Healthy; Genus correlation",
                                                            filename="Result_figures/correlation_analysis/networks/Community_Gold_Star___Remote__Healthy___genus_correlation_network.pdf",
                                                            myseed = 1, 
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = "pdf")

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# GENUS, Community__Gold_Star___Remote__Hx_Current_OM_URTI
genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Hx_Current_OM_URTI___genus___correlation.tsv",
>>>>>>> ac8d1af8b0bcc3ccad540a152432b91af0083eb7
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Genus_pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
file_type <- "pdf"
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                            p_matrix = genus_fastspar_pval.m,
                                                            relabeller_function = first_resolved_taxonomy,
                                                            # relabeller_function = genus_relabeller_network,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.3,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            node_label_segment_colour = "purple",
                                                            label_colour = "black",
                                                            label_size = 4,
                                                            plot_height = 10,
                                                            plot_width = 10,
                                                            edge_width_min = .5,
                                                            edge_width_max = 2.5,
                                                            edge_alpha = 1,
                                                            # network_layout = "stress",
                                                            network_layout = "fr",
                                                            # network_layout = "kk",
                                                            # exclude_to_from_df = edges_to_remove.df,
                                                            # plot_title = paste0(variable, ": ", group, "; Genus correlation"),
                                                            filename= paste0("Result_figures/correlation_analysis/networks/genus/Genus_correlation_network.",file_type),
                                                            myseed = 1,
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = file_type
                                                            )

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                                           vs g__Dolosigranulum
#
# Break down Remote + Otitis status, 
# genus level, 
# abundances for "g__Moraxella", "g__Haemophilus", "g__Corynebacterium","g__Streptococcus" vs "g__Dolosigranulum"
# Include in insets plot_feature_correlations_external() results for g__Dolosigranulum

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

taxa_of_interest <- c("g__Dolosigranulum","g__Moraxella", "g__Haemophilus", "g__Corynebacterium","g__Streptococcus")

# ------------------------------------------------------------------------
# Loop over each result and save to file
genus_cor_files <- list.files("Additional_results/fastspar/")[grepl("___genus___correlation.tsv", list.files("Additional_results/fastspar/"))]
# genus_cor_files <- grep("Community__Otitis_Status.*genus", genus_cor_files,value =T)
genus_cor_files <- grep("^Otitis_Status.*genus", genus_cor_files,value =T)
# genus_cor_files <- grep("^Community___.*genus", genus_cor_files,value =T)

file_type = "svg"
for (cor_file in genus_cor_files){
  genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
                                               sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
  genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
                                                sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
  file_name_split <- strsplit(cor_file, split = "___")[[1]]
  variable <- file_name_split[1]
  group <- file_name_split[2]
  print(variable)
  print(group)
  
  rownames(genus_fastspar_cor.m) <- unlist(lapply(rownames(genus_fastspar_cor.m), first_resolved_taxonomy))
  colnames(genus_fastspar_cor.m) <- unlist(lapply(colnames(genus_fastspar_cor.m), first_resolved_taxonomy))
  rownames(genus_fastspar_pval.m) <- unlist(lapply(rownames(genus_fastspar_pval.m), first_resolved_taxonomy))
  colnames(genus_fastspar_pval.m) <- unlist(lapply(colnames(genus_fastspar_pval.m), first_resolved_taxonomy))
  
  print(dim(genus_fastspar_pval.m))
  plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
                                     feature = "g__Dolosigranulum",
                                     p_value_matrix = genus_fastspar_pval.m,
                                     top_n = 25,
                                     plot_width = 7, plot_height = 7,
                                     include_self = F,
                                     include_title = F,
                                     # filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
                                     filename= paste0("Result_figures/correlation_analysis/by_feature/",variable,"___",group,"___genus_dolosigranulum_correlations_top_25.",file_type),
                                     format = file_type)
  
  genus_fastspar_cor.m <- genus_fastspar_cor.m[taxa_of_interest,taxa_of_interest]
  print(genus_fastspar_cor.m)
  genus_fastspar_pval.m <- genus_fastspar_pval.m[taxa_of_interest,taxa_of_interest]
  
  plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
                                     feature = "g__Dolosigranulum",
                                     p_value_matrix = genus_fastspar_pval.m,
                                     top_n = 10,
                                     plot_width = 4, plot_height = 3,
                                     include_self = F,
                                     include_title = F,
                                     # filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
                                     filename= paste0("Result_figures/correlation_analysis/by_feature/",variable,"___",group,"___genus_dolosigranulum_correlations.",file_type),
                                     format = file_type)
}


genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/Genus_correlation.tsv"),
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/Genus_pvalues.tsv"),
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
rownames(genus_fastspar_cor.m) <- unlist(lapply(rownames(genus_fastspar_cor.m), first_resolved_taxonomy))
colnames(genus_fastspar_cor.m) <- unlist(lapply(colnames(genus_fastspar_cor.m), first_resolved_taxonomy))
rownames(genus_fastspar_pval.m) <- unlist(lapply(rownames(genus_fastspar_pval.m), first_resolved_taxonomy))
colnames(genus_fastspar_pval.m) <- unlist(lapply(colnames(genus_fastspar_pval.m), first_resolved_taxonomy))
# genus_fastspar_cor.m[taxa_of_interest,]
# genus_fastspar_cor.m <- genus_fastspar_cor.m[taxa_of_interest,taxa_of_interest]
# genus_fastspar_pval.m <- genus_fastspar_pval.m[taxa_of_interest,taxa_of_interest]
#
# source("code/helper_functions.R")
# plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
#                                    feature = "g__Dolosigranulum",
#                                    p_value_matrix = genus_fastspar_pval.m,
#                                    top_n = 100,
#                                    plot_width = 30, plot_height = 25,
#                                    include_self = F,
#                                    include_title = F,
#                                    filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
#                                    format = "pdf")

# ------------------------------------------------------------------------
# BELOW FOR TESTING, NOT FOR PUBLICATION

# # genus_palette <- setNames(c("#a84b54","#a2b432", "#c057c5","#6f64cf","#5aae36"),
# #                           c("g__Dolosigranulum", "g__Moraxella", "g__Haemophilus","g__Corynebacterium", "g__Streptococcus"))
# # 
# # genus.df <- read.csv("Result_tables/count_tables/Genus_counts.csv", header =T)
# # # genus.df$taxonomy_genus %in% rownames(genus_fastspar_cor.m)
# # genus.df <- m2df(clr(df2matrix(genus.df)), "taxonomy_genus")
# # genus_melt.df <- melt(genus.df,variable.name = "Sample", value.name = "Counts")
# # genus_melt.df <- genus_melt.df[genus_melt.df$Sample %in% metadata.df$Index,]
# # genus_melt.df <- left_join(genus_melt.df, metadata.df, by = c("Sample" = "Index"))
# # genus_melt.df$Label <- unlist(lapply(genus_melt.df$taxonomy_genus, first_resolved_taxonomy))
# # genus_melt.df <- subset(genus_melt.df, Label %in% taxa_of_interest)
# # 
# # # sample delosi
# # temp <- genus_melt.df[genus_melt.df$Label == "g__Dolosigranulum",]
# # temp <- dcast(temp, Sample~Label, value.var = "Counts")
# # 
# # temp2 <- genus_melt.df[genus_melt.df$Label != "g__Dolosigranulum",]
# # genus_melt.df <- left_join(temp2, temp, by = "Sample")
# # 
# # genus_melt.df$g__Dolosigranulum <- log10(genus_melt.df$g__Dolosigranulum)
# # genus_melt.df$Counts <- log10(genus_melt.df$Counts)
# # genus_melt.df$Counts[is.infinite(genus_melt.df$Counts)] <- 0
# # 
# # 
# # remote_effusion.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__Effusion")
# # remote_hxom.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__HxOM")
# # remote_neverom.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__Never OM")
# # remote_perforation.df <- subset(genus_melt.df, Community__Otitis_Status == "Remote__Perforation")
# # rural_effusion.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__Effusion")
# # rural_hxom.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__HxOM")
# # rural_neverom.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__Never OM")
# # rural_perforation.df <- subset(genus_melt.df, Community__Otitis_Status == "Rural__Perforation")
# # 
# # remote_effusion_plot <- 
# #   ggplot(remote_effusion.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # remote_hxom_plot <- 
# #   ggplot(remote_hxom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # remote_neverom_plot <- 
# #   ggplot(remote_neverom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # remote_perforation_plot <- 
# #   ggplot(remote_perforation.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # rural_effusion_plot <- 
# #   ggplot(rural_effusion.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # rural_hxom_plot <- 
# #   ggplot(rural_hxom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # rural_neverom_plot <- 
# #   ggplot(rural_neverom.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # rural_perforation_plot <- 
# #   ggplot(rural_perforation.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) +
# #   geom_point() +
# #   xlab("log10(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[10]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) + scale_fill_manual(values = genus_palette) + scale_colour_manual(values = genus_palette) +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),colour=guide_legend(title="Genus"), shape=guide_legend(title="Genus")) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() + theme(legend.title.align = 0.5)
# # 
# # base <- "Result_figures/correlation_analysis/"
# # 
# # filetype = "svg"
# # ggsave(plot = remote_effusion_plot,filename = paste0(base, "remote_effusion_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# # ggsave(plot = remote_hxom_plot,filename = paste0(base, "remote_hxom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# # ggsave(plot = remote_neverom_plot,filename = paste0(base, "remote_neverom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# # ggsave(plot = remote_perforation_plot,filename = paste0(base, "remote_perforation_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# # ggsave(plot = rural_effusion_plot,filename = paste0(base, "rural_effusion_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# # ggsave(plot = rural_hxom_plot,filename = paste0(base, "rural_hxom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# # ggsave(plot = rural_neverom_plot,filename = paste0(base, "rural_neverom_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# # ggsave(plot = rural_perforation_plot,filename = paste0(base, "rural_perforation_delosi_correlations.",filetype),width = 6, height = 4, device = filetype)
# 
# 
# 
# # 
# # 
# # myplot <- 
# #   ggplot(genus_melt.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   # ggplot(genus_melt.df, aes(x = g__Dolosigranulum, y = Counts, fill = Label,shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
# #   geom_point() +
# #   facet_wrap(~Community__Otitis_Status, scales = "free_x", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),
# #          colour=guide_legend(title="Genus"),
# #          shape=guide_legend(title="Genus")) +
# #   xlab("log2(Dolosigranulum CLR transformed counts)") +
# #   ylab(expression(paste(log[2]~"(Genus CLR transformed counts)"))) + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) +
# #   scale_fill_manual(values = genus_palette) +
# #   scale_colour_manual(values = genus_palette) +
# #   # scale_x_continuous(limits = c(2.5,20), breaks = seq(0,20,2.5)) +
# #   # scale_y_continuous(limits = c(0,20), breaks = seq(0,20,2.5)) +
# #   scale_x_continuous(limits = c(0.3,1.2)) +
# #   scale_y_continuous(limits = c(-0.05,2.5)) +
# #   theme_bw() +
# #   theme(legend.title.align = 0.5)
# # myplot
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.pdf",
# #        plot = myplot,
# #        device = "pdf",
# #        height = 15,
# #        width = 10)
# # 
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.svg",
# #        plot = myplot,
# #        device = "svg",
# #        height = 15,
# #        width = 10)
# 
# 
# 
# 
# 
# # Create all abundance plots
# # genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header =T)
# # genus_rel_melt.df <- melt(genus_rel.df,variable.name = "Sample", value.name = "Relative_abundance")
# # genus_rel_melt.df <- genus_rel_melt.df[genus_rel_melt.df$Sample %in% metadata.df$Index,]
# # genus_rel_melt.df <- left_join(genus_rel_melt.df, metadata.df, by = c("Sample" = "Index"))
# # genus_rel_melt.df$Label <- unlist(lapply(genus_rel_melt.df$taxonomy_genus, first_resolved_taxonomy))
# # genus_rel_melt.df <- subset(genus_rel_melt.df, Label %in% taxa_of_interest)
# # 
# # # sample delosi
# # temp <- genus_rel_melt.df[genus_rel_melt.df$Label == "g__Dolosigranulum",]
# # temp <- dcast(temp, Sample~Label, value.var = "Relative_abundance")
# # 
# # temp2 <- genus_rel_melt.df[genus_rel_melt.df$Label != "g__Dolosigranulum",]
# # genus_rel_melt.df <- left_join(temp2, temp, by = "Sample")
# # 
# # 
# # # ------------------------------------------------------------------------------
# # 
# # 
# # myplot <- 
# #   ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
# #                               fill = Label,
# #                               shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
# #   geom_point() +
# #   facet_wrap(Otitis_Status~Community, scales = "free", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),
# #          colour=guide_legend(title="Genus"),
# #          shape=guide_legend(title="Genus")) +
# #   xlab("Dolosigranulum relative abundance") +
# #   ylab("Genus relative abundance") + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) +
# #   scale_fill_manual(values = genus_palette) +
# #   scale_colour_manual(values = genus_palette) +
# #   scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   theme_bw() +
# #   theme(legend.title.align = 0.5)
# # myplot
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.pdf",
# #        plot = myplot,
# #        device = "pdf",
# #        height = 15,
# #        width = 10)
# # 
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_otitis_status_correlations.svg",
# #        plot = myplot,
# #        device = "svg",
# #        height = 15,
# #        width = 10)
# # 
# # 
# # myplot <- 
# #   ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
# #                                 fill = Label,
# #                                 shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
# #   geom_point() +
# #   facet_wrap(~Otitis_Status, scales = "free", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),
# #          colour=guide_legend(title="Genus"),
# #          shape=guide_legend(title="Genus")) +
# #   xlab("Dolosigranulum relative abundance") +
# #   ylab("Genus relative abundance") + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) +
# #   scale_fill_manual(values = genus_palette) +
# #   scale_colour_manual(values = genus_palette) +
# #   scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   theme_bw() +
# #   theme(legend.title.align = 0.5)
# # myplot
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_otitis_status_correlations.pdf",
# #        plot = myplot,
# #        device = "pdf",
# #        height = 8,
# #        width = 10)
# # 
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_otitis_status_correlations.svg",
# #        plot = myplot,
# #        device = "svg",
# #        height = 8,
# #        width = 10)
# # 
# # 
# # myplot <- 
# #   ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
# #                                 fill = Label,
# #                                 shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
# #   geom_point() +
# #   facet_wrap(~Community, scales = "free", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),
# #          colour=guide_legend(title="Genus"),
# #          shape=guide_legend(title="Genus")) +
# #   xlab("Dolosigranulum relative abundance") +
# #   ylab("Genus relative abundance") + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) +
# #   scale_fill_manual(values = genus_palette) +
# #   scale_colour_manual(values = genus_palette) +
# #   scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   theme_bw() +
# #   theme(legend.title.align = 0.5)
# # myplot
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_correlations.pdf",
# #        plot = myplot,
# #        device = "pdf",
# #        height = 5,
# #        width = 10)
# # 
# # ggsave(filename = "Result_figures/correlation_analysis/dolosigranulum_vs_genus_community_correlations.svg",
# #        plot = myplot,
# #        device = "svg",
# #        height = 5,
# #        width = 10)
# # 
# # 
# # 
# # myplot <- 
# #   ggplot(genus_melt.df, aes(x = log10(g__Dolosigranulum), y = log10(Counts), 
# #                                 fill = Label,
# #                                 shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
# #   geom_point() +
# #   facet_wrap(~Community, scales = "free", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),
# #          colour=guide_legend(title="Genus"),
# #          shape=guide_legend(title="Genus")) +
# #   xlab("Dolosigranulum relative abundance") +
# #   ylab("Genus relative abundance") + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) +
# #   scale_fill_manual(values = genus_palette) +
# #   scale_colour_manual(values = genus_palette) +
# #   # scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   # scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   theme_bw() +
# #   theme(legend.title.align = 0.5)
# # myplot
# # 
# # 
# # myplot <- 
# #   ggplot(genus_rel_melt.df, aes(x = g__Dolosigranulum*100, y = Relative_abundance*100, 
# #                             fill = Label,
# #                             shape = Label)) +
# #   geom_smooth( method= "lm",se =F, aes(colour = Label), lwd = .5) + 
# #   geom_point() +
# #   # facet_wrap(~Community, scales = "free", ncol = 2) +
# #   guides(fill=guide_legend(title="Genus"),
# #          colour=guide_legend(title="Genus"),
# #          shape=guide_legend(title="Genus")) +
# #   xlab("Dolosigranulum relative abundance") +
# #   ylab("Genus relative abundance") + 
# #   scale_shape_manual(values = c(25,24,23,22,21)) +
# #   scale_fill_manual(values = genus_palette) +
# #   scale_colour_manual(values = genus_palette) +
# #   scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
# #   theme_bw() +
# #   theme(legend.title.align = 0.5)
# # myplot
# 
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# # GENUS, Community___Remote
# genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community___Remote___genus___correlation.tsv",
#                                            sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
# genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community___Remote___genus___pvalues.tsv",
#                                             sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
# library(tidygraph)
# library(ggraph)
# genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
#                                                           p_matrix = genus_fastspar_pval.m,
#                                                           relabeller_function = first_resolved_taxonomy,
#                                                           
#                                                           p_value_threshold = 0.05,
#                                                           cor_threshold = 0.4,
#                                                           node_size = 4,
#                                                           node_colour = "grey20",
#                                                           node_fill = "grey20",
#                                                           label_colour = "black",
#                                                           label_size = 3,
#                                                           plot_height = 10,
#                                                           plot_width = 10,
#                                                           edge_width_min = .5,
#                                                           edge_width_max = 2.5,
#                                                           network_layout = "fr",
#                                                           # exclude_to_from_df = edges_to_remove.df,
#                                                           plot_title = "Community: Remote; Genus correlation",
#                                                           filename="Result_figures/correlation_analysis/networks/Community___Remote___genus_correlation_network.pdf",
#                                                           myseed = 1, 
#                                                           edgetype = "link",
#                                                           show_p_label = F,
#                                                           file_type = "pdf")
# genus_correlation_network.l$network_plot
# 
# source("code/helper_functions.R")
# # heatmap(genus_fastspar_cor.m)
# # plot_corrplot(correlation_matrix = genus_fastspar_cor.m,
# #               p_value_matrix = genus_fastspar_pval.m,
# #               p_value_threshold = .05,
# #               relabeller_function = first_resolved_taxonomy,
# #               label_size = .5,
# #               make_insig_na = F,
# #               order = "original",
# #               file_type = "pdf",
# #               # filename = "out.pdf",
# #               plot_height = 10,
# #               plot_width = 10)
# 
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# # GENUS, Community___Rural
# genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community___Rural___genus___correlation.tsv",
#                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
# genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community___Rural___genus___pvalues.tsv",
#                                               sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
# 
# # dim(subset(metadata.df, Community == "Remote"))
# # dim(subset(metadata.df, Community == "Rural"))
# genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
#                                                             p_matrix = genus_fastspar_pval.m,
#                                                             relabeller_function = first_resolved_taxonomy,
#                                                             
#                                                             p_value_threshold = 0.05,
#                                                             cor_threshold = 0.4,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             label_colour = "black",
#                                                             label_size = 3,
#                                                             plot_height = 10,
#                                                             plot_width = 10,
#                                                             edge_width_min = .5,
#                                                             edge_width_max = 2.5,
#                                                             network_layout = "fr",
#                                                             # exclude_to_from_df = edges_to_remove.df,
#                                                             plot_title = "Community: Rural; Genus correlation",
#                                                             filename="Result_figures/correlation_analysis/networks/Community___Rural___genus_correlation_network.pdf",
#                                                             myseed = 1, 
#                                                             edgetype = "link",
#                                                             show_p_label = F,
#                                                             file_type = "pdf")
# 
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# # GENUS, Community__Gold_Star___Remote__Healthy
# genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Healthy___genus___correlation.tsv",
#                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
# genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Healthy___genus___pvalues.tsv",
#                                               sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
# 
# # dim(subset(metadata.df, Community == "Remote"))
# # dim(subset(metadata.df, Community == "Rural"))
# genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
#                                                             p_matrix = genus_fastspar_pval.m,
#                                                             relabeller_function = first_resolved_taxonomy,
#                                                             
#                                                             p_value_threshold = 0.05,
#                                                             cor_threshold = 0.4,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             label_colour = "black",
#                                                             label_size = 3,
#                                                             plot_height = 10,
#                                                             plot_width = 10,
#                                                             edge_width_min = .5,
#                                                             edge_width_max = 2.5,
#                                                             network_layout = "fr",
#                                                             # exclude_to_from_df = edges_to_remove.df,
#                                                             plot_title = "Community__Gold_Star: Remote__Healthy; Genus correlation",
#                                                             filename="Result_figures/correlation_analysis/networks/Community_Gold_Star___Remote__Healthy___genus_correlation_network.pdf",
#                                                             myseed = 1, 
#                                                             edgetype = "link",
#                                                             show_p_label = F,
#                                                             file_type = "pdf")
# 
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# # GENUS, Community__Gold_Star___Remote__Hx_Current_OM_URTI
# genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Hx_Current_OM_URTI___genus___correlation.tsv",
#                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
# genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Hx_Current_OM_URTI___genus___pvalues.tsv",
#                                               sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
# 
# # dim(subset(metadata.df, Community == "Remote"))
# # dim(subset(metadata.df, Community == "Rural"))
# genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
#                                                             p_matrix = genus_fastspar_pval.m,
#                                                             relabeller_function = first_resolved_taxonomy,
#                                                             
#                                                             p_value_threshold = 0.05,
#                                                             cor_threshold = 0.4,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             label_colour = "black",
#                                                             label_size = 3,
#                                                             plot_height = 10,
#                                                             plot_width = 10,
#                                                             edge_width_min = .5,
#                                                             edge_width_max = 2.5,
#                                                             network_layout = "fr",
#                                                             # exclude_to_from_df = edges_to_remove.df,
#                                                             plot_title = "Community__Gold_Star: Remote__Hx_Current_OM_URTI; Genus correlation",
#                                                             filename="Result_figures/correlation_analysis/networks/Community_Gold_Star___Remote__Hx_Current_OM_URTI_genus_correlation_network.pdf",
#                                                             myseed = 1, 
#                                                             edgetype = "link",
#                                                             show_p_label = F,
#                                                             file_type = "pdf")
# 
# 
# 
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------------
# otu_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/OTU_correlation.tsv",
#                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
# otu_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/OTU_pvalues.tsv",
#                                               sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
# 
# 
# # combined_otu_labeller <- function(x){
# #   print(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
# #   first_resolved_taxonomy(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
# # }
# # rownames(otu_fastspar_pval.m) %in% otu_taxonomy_map.df$OTU.ID
# # plot_corrplot(correlation_matrix = otu_fastspar_cor.m,
# #               p_value_matrix = otu_fastspar_pval.m,
# #               p_value_threshold = 0.01,
# #               
# #               relabeller_function = otu_relabeller_function)
# source("code/helper_functions.R")
# otu_correlation_network.l <- generate_correlation_network(cor_matrix = otu_fastspar_cor.m,
#                                                             p_matrix = otu_fastspar_pval.m,
#                                                             relabeller_function = otu_relabeller_function,
#                                                             p_value_threshold = 0.05,
#                                                             cor_threshold = 0.5,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             label_colour = "black",
#                                                             label_size = 3,
#                                                             plot_height = 10,
#                                                             plot_width = 10,
#                                                             edge_width_min = .5,
#                                                             edge_width_max = 2.5,
#                                                             network_layout = "fr",
#                                                             # exclude_to_from_df = edges_to_remove.df,
#                                                             # filename="Result_figures/correlation_analysis/networks/genus_filtered_fastspar_cor_network.pdf",
#                                                             myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")
# 
# otu_correlation_network.l$network_plot
# 
# 
# genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Genus_correlation.tsv",
#                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
# genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Genus_pvalues.tsv",
#                                               sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
# source("code/helper_functions.R")
# 
# 
# temp_cor <- genus_fastspar_cor.m
# diag(temp_cor) <- NA
# dim(temp_cor)
# dim(genus_fastspar_cor.m)
# cor_retain <- rownames(temp_cor[apply(temp_cor, 1,function(x) max(x,na.rm = T)) >= 0.5,])
# temp_cor <- temp_cor[cor_retain,cor_retain]
# 
# genus_fastspar_pval.m <- genus_fastspar_pval.m[cor_retain,cor_retain]
# p_retain <- rownames(genus_fastspar_pval.m[apply(genus_fastspar_pval.m,1,function(x)min(x,na.rm=T)) <= 0.05,])
# genus_fastspar_pval.m <- genus_fastspar_pval.m[p_retain,p_retain]
# genus_fastspar_cor.m <- genus_fastspar_cor.m[rownames(genus_fastspar_pval.m),rownames(genus_fastspar_pval.m)]
# 
# # pdf("out.pdf",width = 10, height = 10)
# # corrplot(genus_fastspar_cor.m,p.mat = genus_fastspar_pval.m,sig.level = 0.05,tl.cex = .3,mar=c(0,0,0,0),type = "lower",tl.col = "black")
# # dev.off()
# dim(genus_fastspar_cor.m)
# dim(genus_fastspar_pval.m)
# source("code/helper_functions.R")
# plot_corrplot(correlation_matrix = genus_fastspar_cor.m,
#               p_value_matrix = genus_fastspar_pval.m,
#               p_value_threshold = .01,
#               relabeller_function = genus_relabeller_function,
#               label_size = .8,
#               make_insig_na = T,
#               order = "original",
#               file_type = "pdf",
#               filename = "out.pdf",
#               plot_height = 10,
#               plot_width = 10)
# dev.off()
# genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
#                                                             p_matrix = genus_fastspar_pval.m,
#                                                             relabeller_function = first_resolved_taxonomy,
#                                                             p_value_threshold = 0.05,
#                                                             cor_threshold = 0.5,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             label_colour = "black",
#                                                             label_size = 3,
#                                                             plot_height = 10,
#                                                             plot_width = 10,
#                                                             edge_width_min = .5,
#                                                             edge_width_max = 2.5,
#                                                             network_layout = "fr",
#                                                             # exclude_to_from_df = edges_to_remove.df,
#                                                             # filename="Result_figures/correlation_analysis/networks/genus_filtered_fastspar_cor_network.pdf",
#                                                             myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")
# 
# genus_correlation_network.l$network_plot
# 
