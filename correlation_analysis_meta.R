# Correlation for Otitis culture paper. Separate to 16S study.
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

library(devtools)
# install_github("zdk123/SpiecEasi")
# library(SpiecEasi)
library(Matrix)
# library(igraph)
# install.packages("psych")
library(psych)
# install.packages("ggraph")
library(ggraph)
library(tidygraph)

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")

for (i in c("Gold_0","Gold_1","Nose_0","Nose_1","Nose_2","All_samples", "Nose_0_1","Nose_0_2","Nose_1_2","Nose_12_combined", "Nose_0_12_combined")){
  dir.create(file.path(paste0("culture_paper_analysis/results/",i, "/by_feature")),recursive = T)
  dir.create(file.path(paste0("culture_paper_analysis/results/",i, "/networks")),recursive = T)
}

# Load the data
data.m <- as.matrix(read.table("data/Culture_paper/Culture_Viral_qPCR_binary_results.tsv", sep ="\t", header = T,row.names = 1))
metadata.df <- read.table("data/Culture_paper/Culture_Viral_qPCR_metadata.tsv", sep ="\t", header = T)
data_cleaned_names.df <- read.table("data/Culture_paper/Culture_Viral_qPCR_Taxonomy_clean_publication.tsv", sep = "\t", header = T)

rownames(data_cleaned_names.df) <- data_cleaned_names.df$Original_no_space
rownames(data.m) <- data_cleaned_names.df[rownames(data.m),]$Cleaned.up.version

# metadata.df <- read.table("data/metadata.tsv", sep ="\t", header = T)
metadata.df$NAME_2 <- gsub("-", ".", metadata.df$NAME)
rownames(metadata.df) <- metadata.df$NAME_2
data.m <- rbind(data.m, Gold_Star = metadata.df[colnames(data.m),]$Gold_Star)
data.m <- rbind(data.m, Nose = metadata.df[colnames(data.m),]$Nose)
rownames(data.m)[rownames(data.m) == "Gold_Star"] <- "Gold Star"

# --------------------------------------------------------------------------------
# For correlation networks, we want to remove edges between IQR variables for the same species
# Mcat <- as.data.frame(t(combn(grep("Mcat", rownames(data.m),value =T),2)))
# H_inf <- as.data.frame(t(combn(grep("H.inf", rownames(data.m),value =T),2)))
# Spn <- as.data.frame(t(combn(grep("Spn", rownames(data.m),value =T),2)))
Mcat <- as.data.frame(t(combn(grep("M. cat", rownames(data.m),value =T),2)))
H_inf <- as.data.frame(t(combn(grep("H. inf", rownames(data.m),value =T),2)))
Spn <- as.data.frame(t(combn(grep("S. pn", rownames(data.m),value =T),2)))
# Nose <- as.data.frame(t(combn(grep("Nose", rownames(data.m),value =T),2)))
# Gold_star <- as.data.frame(t(combn(grep("Gold", rownames(data.m),value =T),2)))

Nose <- as.data.frame(melt(unique(rownames(data.m), colnames(data.frame))))
Nose$Nose <- "Nose"
names(Nose) <- c("V1", "V2")

Gold_star <- as.data.frame(melt(unique(rownames(data.m), colnames(data.frame))))
Gold_star$Nose <- "Gold Star"
names(Gold_star) <- c("V1", "V2")

edges_to_remove.df <- rbind(Mcat, H_inf,Spn,Nose,Gold_star)
edges_to_remove.df <- rbind(Mcat, H_inf,Spn,Nose,Gold_star)

# ------------------------------------------------------------
feature1 <- "Staphylococcus haemolyticus"
feature2 <- "Staphylococcus hominis"
temp <- calculate_feature_correlations(data.m,
                               feature = feature1,
                               method = "pearson",adjust = "BH")

temp[feature2,]
# isSymmetric.matrix(cor(t(data.m)))
# temp <- calculate_correlation_matrix(data.m,method = "pearson", adjust = "BH")
# temp$cor_matrix[feature2,feature1]
# isSymmetric.matrix(round(temp$cor_padj_matrix,10))
# temp$cor_pval_matrix[feature2,feature1]
# temp$cor_padj_matrix[feature2,feature1]
# temp$cor_padj_matrix[feature1,feature2]

# ------------------------------------------------------------

# data.m["N_RSV_B",]
# data.m["N_FLU_B",]
# data.m["N_FLU_A",]

# zv <- apply(data.m, 1, function(x) length(unique(x)) == 1)
# data.m <- data.m[!zv, ]

# ------------------------------------


## Nose
# Network plot; 1 for each value – 1,2, and 3 (Nose) DONE
# Network plot; 1 for each value of the data transformed into a binary feature… 0=asymptomatic & 1=symptomatic (combination of the 2 & 3 values) DONE
  # 0 plot not required as done already 
# Correlation bar plot; 1 vs 2, and 1 vs 3 (by feature for each)
# Correlation bar plot; 0 (asymptomatic) vs 1 (symptomatic) (by feature for 01)
# ------------------------------------
calculate_prevalences <- function(mymatrix){
  apply(mymatrix, 1, function(x) {length(which(x > 0))}) /length(colnames(mymatrix))
}
prevalence_threshold = 0.2

# Create additional datasets

# Gold star
data_g0.m <- data.m[,data.m["Gold Star",] == 0]
data_g1.m <- data.m[,data.m["Gold Star",] == 1]

# Nose
data_n0.m <- data.m[,data.m["Nose",] == 0]
data_n1.m <- data.m[,data.m["Nose",] == 1]
data_n2.m <- data.m[,data.m["Nose",] == 2]

# Combinations of nose groups
data_n01.m <- data.m[,data.m["Nose",] %in% c(0,1)]
data_n02.m <- data.m[,data.m["Nose",] %in% c(0,2)]
data_n12.m <- data.m[,data.m["Nose",] %in% c(1,2)]

# Compare nose 0 against 1+2
data_n0_12_combined.m <- data.m
data_n0_12_combined.m["Nose",][data_n0_12_combined.m["Nose",] %in% c(1,2)] <- 1

# Nose 1 and 2 samples, but with 1 and 2 combined
data_n12_combined.m <- data_n0_12_combined.m[,data_n0_12_combined.m["Nose",] == 1]

# round(melt(calculate_prevalences(data_g1.m)) *100, 2)
# dim(data_n2.m) * .3
# melt(apply(data.m, 1, sum))
# melt(apply(data_n2.m, 1, sum))
# dim(data.m)

# Now that the individual datasets have been created, remove variables that are not prevalent enough
data.m <- data.m[calculate_prevalences(data.m) > prevalence_threshold,,drop =F]
data_g0.m <- data_g0.m[calculate_prevalences(data_g0.m) > prevalence_threshold,,drop =F]
data_g1.m <- data_g1.m[calculate_prevalences(data_g1.m) > prevalence_threshold,,drop =F]
data_n0.m <- data_n0.m[calculate_prevalences(data_n0.m) > prevalence_threshold,,drop =F]
data_n1.m <- data_n1.m[calculate_prevalences(data_n1.m) > prevalence_threshold,,drop =F]
data_n2.m <- data_n2.m[calculate_prevalences(data_n2.m) > prevalence_threshold,,drop =F]
data_n01.m <- data_n01.m[calculate_prevalences(data_n01.m) > prevalence_threshold,,drop =F]
data_n02.m <- data_n02.m[calculate_prevalences(data_n02.m) > prevalence_threshold,,drop =F]
data_n12.m <- data_n12.m[calculate_prevalences(data_n12.m) > prevalence_threshold,,drop =F]
data_n0_12_combined.m <- data_n0_12_combined.m[calculate_prevalences(data_n0_12_combined.m) > prevalence_threshold,,drop =F]
data_n12_combined.m <- data_n12_combined.m[calculate_prevalences(data_n12_combined.m) > prevalence_threshold,,drop =F]


# Just for correlation bar plots
# Correlation bar plot; 1 vs 2, and 1 vs 3 (by feature for each)
# Correlation bar plot; 0 (asymptomatic) vs 1 (symptomatic) (by feature for 01)
# data_just_for_feature_analysis_n01.m <- data_n01.m[!rownames(data_n01.m) == "Nose", ]
# data_just_for_feature_analysis_n02.m <- data_n02.m[!rownames(data_n02.m) == "Nose", ]
# data_just_for_feature_analysis_n0_12_combined.m <- data_n0_12_combined.m[!rownames(data_n0_12_combined.m) == "Nose", ]
# 
# data_just_for_feature_analysis_n01.m <- rbind(data_just_for_feature_analysis_n01.m, Nose = metadata.df[colnames(data_just_for_feature_analysis_n01.m),]$Nose)
# data_just_for_feature_analysis_n02.m <- rbind(data_just_for_feature_analysis_n02.m, Nose = metadata.df[colnames(data_just_for_feature_analysis_n02.m),]$Nose)
# data_just_for_feature_analysis_n0_12_combined.m <- rbind(data_just_for_feature_analysis_n0_12_combined.m, Nose = metadata.df[colnames(data_just_for_feature_analysis_n0_12_combined.m),]$Nose)
# data_just_for_feature_analysis_n0_12_combined.m["Nose",][data_just_for_feature_analysis_n0_12_combined.m["Nose",] %in% c(1,2)] <- 1

# ------------------------------------------------------------------------------------
# Calculate correlation results and generate correlation plots for each variable
generate_feature_results <- function(mydata.m, prefix =""){
  for (feature in rownames(mydata.m)){
    feature_filename <- gsub("/_", "_", feature)
    feature_filename <- gsub("/", "_", feature_filename)
    feature_filename <- gsub(" ", "_", feature_filename)
    calculate_feature_correlations(mydata.m,
                                   feature = feature,
                                   filename = paste0("culture_paper_analysis/results/",
                                                     prefix, feature_filename, "_feature_correlations.csv"),
                                   method = "pearson")
    

    feature_plot_filemame <- paste0("culture_paper_analysis/results/",prefix, feature_filename, "_feature_correlations.pdf")
    plot_feature_correlations(mydata.m, feature = feature,top_n = 25, method = "pearson",filename = feature_plot_filemame,
                              plot_width = 25, plot_height = 25)
  }
}
# For all samples
generate_feature_results(data.m, prefix = "All_samples/by_feature/All_")

# For sets of samples
# datasets <- list(data_g0.m,data_g1.m,data_n0.m,data_n1.m,data_n2.m)
generate_feature_results(data_g0.m, prefix = "Gold_0/by_feature/Gold_0_")
generate_feature_results(data_g1.m, prefix = "Gold_1/by_feature/Gold_1_")
generate_feature_results(data_n0.m, prefix = "Nose_0/by_feature/Nose_0_")
generate_feature_results(data_n1.m, prefix = "Nose_1/by_feature/Nose_1_")
generate_feature_results(data_n2.m, prefix = "Nose_2/by_feature/Nose_2_")

generate_feature_results(data_n01.m, prefix = "Nose_0_1/by_feature/Nose_0_1_")
generate_feature_results(data_n02.m, prefix = "Nose_0_2/by_feature/Nose_0_2_")
generate_feature_results(data_n12.m, prefix = "Nose_1_2/by_feature/Nose_1_2_")

generate_feature_results(data_n12_combined.m, prefix = "Nose_12_combined/by_feature/Nose_12_combined_")
generate_feature_results(data_n0_12_combined.m, prefix = "Nose_0_12_combined/by_feature/Nose_0_12_combined_")


# Calculate correlation matrix and p-values
correlation_results <- calculate_correlation_matrix(data.m, method = "pearson", adjust = "BH")

correlation_results_g0 <- calculate_correlation_matrix(data_g0.m, method = "pearson", adjust = "BH")
correlation_results_g1 <- calculate_correlation_matrix(data_g1.m, method = "pearson", adjust = "BH")

correlation_results_n0 <- calculate_correlation_matrix(data_n0.m, method = "pearson", adjust = "BH")
correlation_results_n1 <- calculate_correlation_matrix(data_n1.m, method = "pearson", adjust = "BH")
correlation_results_n2 <- calculate_correlation_matrix(data_n2.m, method = "pearson", adjust = "BH")

correlation_results_n01 <- calculate_correlation_matrix(data_n01.m, method = "pearson", adjust = "BH")
correlation_results_n02 <- calculate_correlation_matrix(data_n02.m, method = "pearson", adjust = "BH")
correlation_results_n12 <- calculate_correlation_matrix(data_n12.m, method = "pearson", adjust = "BH")
correlation_results_n12_combined <- calculate_correlation_matrix(data_n12_combined.m, method = "pearson", adjust = "BH")
correlation_results_n0_12_combined <- calculate_correlation_matrix(data_n0_12_combined.m, method = "pearson", adjust = "BH")

cor.m <- correlation_results$cor_matrix
cor_pval.m <- correlation_results$cor_padj_matrix

cor_g0.m <- correlation_results_g0$cor_matrix
cor_pval_g0.m <- correlation_results_g0$cor_padj_matrix

cor_g1.m <- correlation_results_g1$cor_matrix
cor_pval_g1.m <- correlation_results_g1$cor_padj_matrix

# cor_pval2.m <- correlation_results$cor_pval_matrix
# cor_pval2.m[feature2, feature1]
# cor_pval.m[feature2, feature1]
# calculate_feature_correlations(data.m,
#                                feature = feature1,
#                                method = "pearson",adjust = "BH")[feature2,]

cor_n0.m <- correlation_results_n0$cor_matrix
cor_pval_n0.m <- correlation_results_n0$cor_padj_matrix

cor_n1.m <- correlation_results_n1$cor_matrix
cor_pval_n1.m <- correlation_results_n1$cor_padj_matrix

cor_n2.m <- correlation_results_n2$cor_matrix
cor_pval_n2.m <- correlation_results_n2$cor_padj_matrix

cor_n01.m <- correlation_results_n01$cor_matrix
cor_pval_n01.m <- correlation_results_n01$cor_padj_matrix

cor_n02.m <- correlation_results_n02$cor_matrix
cor_pval_n02.m <- correlation_results_n02$cor_padj_matrix

cor_n12.m <- correlation_results_n12$cor_matrix
cor_pval_n12.m <- correlation_results_n12$cor_padj_matrix

cor_n12_combined.m <- correlation_results_n12_combined$cor_matrix
cor_pval_n12_combined.m <- correlation_results_n12_combined$cor_padj_matrix

cor_n0_12_combined.m <- correlation_results_n0_12_combined$cor_matrix
cor_pval_n0_12_combined.m <- correlation_results_n0_12_combined$cor_padj_matrix

# ------------------------------------------------------------------------------------
# Create network plots

# Network all samples
correlation_network.l <- generate_correlation_network(cor_matrix = cor.m,
                                                      p_matrix = cor_pval.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/All_samples/networks/All_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Gold 0
correlation_network.l <- generate_correlation_network(cor_matrix = cor_g0.m,
                                                      p_matrix = cor_pval_g0.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Gold_0/networks/Gold_0_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Gold 1
source("code/helper_functions.R")
correlation_network.l <- generate_correlation_network(cor_matrix = cor_g1.m,
                                                      p_matrix = cor_pval_g1.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Gold_1/networks/Gold_1_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Nose 0
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n0.m,
                                                      p_matrix = cor_pval_n0.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_0/networks/Nose_0_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Nose 1
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n1.m,
                                                      p_matrix = cor_pval_n1.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_1/networks/Nose_1_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Nose 2
source("code/helper_functions.R")
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n2.m,
                                                      p_matrix = cor_pval_n2.m,
                                                      p_value_threshold = .05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_2/networks/Nose_2_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")
# plot_corrplot(cor_n2.m,label_size = .4)

# Nose 01 samples
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n01.m,
                                                      p_matrix = cor_pval_n01.m,
                                                      p_value_threshold = .05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_0_1/networks/Nose_0_1_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Nose 02 samples
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n02.m,
                                                      p_matrix = cor_pval_n02.m,
                                                      p_value_threshold = .05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_0_2/networks/Nose_0_2_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Nose 12 samples
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n12.m,
                                                      p_matrix = cor_pval_n12.m,
                                                      p_value_threshold = .05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_1_2/networks/Nose_1_2_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Nose 12 combined samples
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n12_combined.m,
                                                      p_matrix = cor_pval_n12_combined.m,
                                                      p_value_threshold = .05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_12_combined/networks/Nose_12_combined_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

# Nose 0 and 12 combined samples
correlation_network.l <- generate_correlation_network(cor_matrix = cor_n0_12_combined.m,
                                                      p_matrix = cor_pval_n0_12_combined.m,
                                                      p_value_threshold = .05,
                                                      cor_threshold = 0.3,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      edge_width_min = .5,
                                                      edge_width_max = 2.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/Nose_0_12_combined/networks/Nose_0_12_combined_correlation_graph.pdf",
                                                      myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")


# TODO make corrplots for each datasets
# TODO get numbers per dataset
plot_corrplot(correlation_matrix = cor_g1.m,label_size = .5,
              p_value_matrix = cor_pval_g1.m,p_value_threshold = 0.05)

