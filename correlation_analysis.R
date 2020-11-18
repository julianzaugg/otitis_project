# Calculate the correlations in abundances between taxa

library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
library(igraph)
# install.packages("psych")
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

# Load combined data (counts, abundances and metadata)
# otu_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata.csv", header = T)
# genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv", header = T)

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# Generate fastspar inputs, only required for group subsets

prepare_input_variables <- c("Community", "Nose", "Otitis_Status", "Community__Gold_Star")

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

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------

otu_cor_files <- list.files("Additional_results/fastspar/")[grepl("___otu___correlation.tsv", list.files("Additional_results/fastspar/"))]
# lapply(rownames(otu_fastspar_cor.m),combined_otu_labeller)
source("code/helper_functions.R")
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
                                                            cor_threshold = 0.5,
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

genus_cor_files <- list.files("Additional_results/fastspar/")[grepl("___genus___correlation.tsv", list.files("Additional_results/fastspar/"))]
genus_cor_files <- grep("Nose|Otitis", genus_cor_files,value =T)
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
                                                              relabeller_function = first_resolved_taxonomy,
                                                              # relabeller_function = genus_relabeller_network,
                                                              p_value_threshold = 0.05,
                                                              cor_threshold = 0.4,
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
                                                              filename= paste0("Result_figures/correlation_analysis/networks/genus/",variable,"___",group,"___genus_correlation_network.pdf"),
                                                              myseed = 1,
                                                              edgetype = "link",
                                                              show_p_label = F,
                                                              file_type = "pdf")
}

genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Rural__Healthy___genus___covariance.tsv",
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Rural__Healthy___genus___pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                            p_matrix = genus_fastspar_pval.m,
                                                            relabeller_function = first_resolved_taxonomy,
                                                            
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.5)
# my_layout <- create_layout(genus_correlation_network.l$network_data, layout = 'igraph', algorithm = 'nicely')
cols_f <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))

ggraph(genus_correlation_network.l$network_data,layout = "stress") +
  geom_edge_link(aes(colour = Correlation, width = abs(Correlation))) +
  geom_node_point(aes(fill = name), shape = 21,size = 3)


colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
edge_width_min = .5
edge_width_max = 2.5
break_length <- length(seq(-1,1,.2))
correlation_graph_plot <- ggraph(genus_correlation_network.l$network_data, layout = "stress") +
  geom_edge_link(aes(colour = Correlation, width = abs(Correlation))) +
  geom_node_point(aes(fill = name), shape = 21,size = 3) +
  scale_edge_width_continuous(name="Correlation", range = c(edge_width_min,edge_width_max),
                              breaks = seq(-1,1,.2)) +
  scale_edge_colour_gradientn(colours = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                                               "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                               "#4393C3", "#2166AC", "#053061")))(11),
                              limits = c(-1,1), # limit colours to full corrrelation range
                              breaks = seq(-1,1,.2), # Break colours from -1 to 1 in steps of 0.2
                              guide = guide_edge_colourbar(barwidth = 0.5, barheight = 10)) +
  scale_fill_manual(values = colour_palette_206_distinct) +
  theme_graph(background = "white")
  edge_widths <- abs(c(rev(rev(seq(-edge_width_max, -edge_width_min, length.out = break_length/2))[-1]),
                       rev(seq(edge_width_max, edge_width_min, length.out = break_length/2))))
  correlation_graph_plot <- correlation_graph_plot + guides(edge_color = guide_legend(override.aes = list(edge_width = edge_widths)),
         edge_width = F)
  correlation_graph_plot








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
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Community__Gold_Star___Remote__Hx_Current_OM_URTI___genus___pvalues.tsv",
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
                                                            plot_title = "Community__Gold_Star: Remote__Hx_Current_OM_URTI; Genus correlation",
                                                            filename="Result_figures/correlation_analysis/networks/Community_Gold_Star___Remote__Hx_Current_OM_URTI_genus_correlation_network.pdf",
                                                            myseed = 1, 
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = "pdf")



# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
otu_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/OTU_correlation.tsv",
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
otu_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/OTU_pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))


# combined_otu_labeller <- function(x){
#   print(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
#   first_resolved_taxonomy(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
# }
# rownames(otu_fastspar_pval.m) %in% otu_taxonomy_map.df$OTU.ID
# plot_corrplot(correlation_matrix = otu_fastspar_cor.m,
#               p_value_matrix = otu_fastspar_pval.m,
#               p_value_threshold = 0.01,
#               
#               relabeller_function = otu_relabeller_function)
source("code/helper_functions.R")
otu_correlation_network.l <- generate_correlation_network(cor_matrix = otu_fastspar_cor.m,
                                                            p_matrix = otu_fastspar_pval.m,
                                                            relabeller_function = otu_relabeller_function,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.5,
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
                                                            # filename="Result_figures/correlation_analysis/networks/genus_filtered_fastspar_cor_network.pdf",
                                                            myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

otu_correlation_network.l$network_plot


genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Genus_correlation.tsv",
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Genus_pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
source("code/helper_functions.R")


temp_cor <- genus_fastspar_cor.m
diag(temp_cor) <- NA
dim(temp_cor)
dim(genus_fastspar_cor.m)
cor_retain <- rownames(temp_cor[apply(temp_cor, 1,function(x) max(x,na.rm = T)) >= 0.5,])
temp_cor <- temp_cor[cor_retain,cor_retain]

genus_fastspar_pval.m <- genus_fastspar_pval.m[cor_retain,cor_retain]
p_retain <- rownames(genus_fastspar_pval.m[apply(genus_fastspar_pval.m,1,function(x)min(x,na.rm=T)) <= 0.05,])
genus_fastspar_pval.m <- genus_fastspar_pval.m[p_retain,p_retain]
genus_fastspar_cor.m <- genus_fastspar_cor.m[rownames(genus_fastspar_pval.m),rownames(genus_fastspar_pval.m)]

# pdf("out.pdf",width = 10, height = 10)
# corrplot(genus_fastspar_cor.m,p.mat = genus_fastspar_pval.m,sig.level = 0.05,tl.cex = .3,mar=c(0,0,0,0),type = "lower",tl.col = "black")
# dev.off()
dim(genus_fastspar_cor.m)
dim(genus_fastspar_pval.m)
source("code/helper_functions.R")
plot_corrplot(correlation_matrix = genus_fastspar_cor.m,
              p_value_matrix = genus_fastspar_pval.m,
              p_value_threshold = .01,
              relabeller_function = genus_relabeller_function,
              label_size = .8,
              make_insig_na = T,
              order = "original",
              file_type = "pdf",
              filename = "out.pdf",
              plot_height = 10,
              plot_width = 10)
dev.off()
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                            p_matrix = genus_fastspar_pval.m,
                                                            relabeller_function = first_resolved_taxonomy,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.5,
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
                                                            # filename="Result_figures/correlation_analysis/networks/genus_filtered_fastspar_cor_network.pdf",
                                                            myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

genus_correlation_network.l$network_plot














# Generate taxonomy summaries
otu_taxa_summary.df <- generate_taxa_summary(mydata = otu_data.df, taxa_column = "OTU.ID")
genus_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df, taxa_column = "taxonomy_genus")
genus_taxa_summary_nose.df <- generate_taxa_summary(mydata = genus_data.df, taxa_column = "taxonomy_genus", group_by_columns = c("Nose"))
genus_taxa_summary_gold_star.df <- generate_taxa_summary(mydata = genus_data.df, taxa_column = "taxonomy_genus", group_by_columns = c("Gold_Star"))
# nose_taxa_summary.df <- generate_taxa_summary(mydata = genus_data.df, taxa_column = "taxonomy_genus", group_by_columns = c("Sample","Nose"))
# write.csv(genus_taxa_summary.df,file = "genus_nose_summary.csv", quote = F, row.names = F)
# write.csv(nose_taxa_summary.df,file = "genus_nose_summary_per_sample.csv", quote = F, row.names = F)

# Filter taxa summaries by prevalence
otu_taxa_summary_filtered.df <- otu_taxa_summary.df %>% filter(Percent_group_samples > 20)
genus_taxa_summary_filtered.df <- genus_taxa_summary.df %>% filter(Percent_group_samples > 20)
genus_taxa_summary_nose_filtered.df <- genus_taxa_summary_nose.df %>% filter(Percent_group_samples > 20)
genus_taxa_summary_gold_star_filtered.df <- genus_taxa_summary_gold_star.df %>% filter(Percent_group_samples > 20)

# Transform read counts
otu_clr.m <-  clr(otu.m)
genus_clr.m <-  clr(genus.m)

# Add metadata to count matrix (optional)
otu_clr_with_meta.m <- rbind(otu_clr.m, t(metadata.df[colnames(otu_clr.m),discrete_variables_to_add_with_counts]))
genus_clr_with_meta.m <- rbind(genus_clr.m, t(metadata.df[colnames(genus_clr.m),discrete_variables_to_add_with_counts]))

# Filter to samples in groups of interest

# Filter taxa, e.g to taxa above a certain prevalence
otu_filtered.m <- otu.m[unique(otu_taxa_summary_filtered.df$OTU.ID),]
otu_clr_filtered.m <- otu_clr.m[unique(otu_taxa_summary_filtered.df$OTU.ID),]

genus_filtered.m <- genus.m[unique(genus_taxa_summary_filtered.df$taxonomy_genus),]
genus_clr_filtered.m <- genus_clr.m[unique(genus_taxa_summary_filtered.df$taxonomy_genus),]

genus_clr_filtered_nose.m <- genus_clr.m[unique(genus_taxa_summary_nose_filtered.df$taxonomy_genus),]
genus_clr_filtered_gold_star.m <- genus_clr.m[unique(genus_taxa_summary_gold_star_filtered.df$taxonomy_genus),]

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# Run Fastspar externally on the raw counts (not the filtered matrices above) and then load the results.
otu_fastspar_cor.m <- as.matrix(read.table("Additional_results/OTU_correlation_fastspar.tsv",
                                                            sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
otu_fastspar_pval.m <- as.matrix(read.table("Additional_results/OTU_pvalues_fastspar.tsv",
                                                             sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))

otu_fastspar_cor_filtered.m  <- otu_fastspar_cor.m[otu_taxa_summary_filtered.df$OTU.ID,
                                                                                 otu_taxa_summary_filtered.df$OTU.ID]

otu_fastspar_pval_filtered.m  <- otu_fastspar_pval.m[otu_taxa_summary_filtered.df$OTU.ID,
                                                                                   otu_taxa_summary_filtered.df$OTU.ID]

genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/Genus_correlation_fastspar.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/Genus_pvalues_fastspar.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))



genus_fastspar_cor_filtered.m  <- genus_fastspar_cor.m[genus_taxa_summary_filtered.df$taxonomy_genus,
                                                                                     genus_taxa_summary_filtered.df$taxonomy_genus]

genus_fastspar_pval_filtered.m  <- genus_fastspar_pval.m[genus_taxa_summary_filtered.df$taxonomy_genus,
                                                                                     genus_taxa_summary_filtered.df$taxonomy_genus]


genus_fastspar_cor_nose_filtered.m  <- genus_fastspar_cor.m[genus_taxa_summary_nose_filtered.df$taxonomy_genus,
                                                                                          genus_taxa_summary_nose_filtered.df$taxonomy_genus]

# rownames(metadata.df[metadata.df$Nose == "0",])
# rownames(metadata.df[metadata.df$Nose == "1",])
# rownames(metadata.df[metadata.df$Gold_Star == "0",])
# rownames(metadata.df[metadata.df$Gold_Star == "1",])

genus_fastspar_pval_nose_filtered.m  <- genus_fastspar_pval.m[genus_taxa_summary_nose_filtered.df$taxonomy_genus,
                                                                                            genus_taxa_summary_nose_filtered.df$taxonomy_genus]

genus_fastspar_cor_gold_star_filtered.m  <- genus_fastspar_cor.m[genus_taxa_summary_gold_star_filtered.df$taxonomy_genus,
                                                                                     genus_taxa_summary_gold_star_filtered.df$taxonomy_genus]

genus_fastspar_pval_gold_star_filtered.m  <- genus_fastspar_pval.m[genus_taxa_summary_gold_star_filtered.df$taxonomy_genus,
                                                                                       genus_taxa_summary_gold_star_filtered.df$taxonomy_genus]

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# Generate correlation network plots

# relabel column and row names to shorter form (have to be unique!)
# colnames(genus_fastspar_cor_filtered.m) <- genus_relabeller_function(colnames(genus_fastspar_cor_filtered.m))
# rownames(genus_fastspar_cor_filtered.m) <- genus_relabeller_function(rownames(genus_fastspar_cor_filtered.m))
# genus_correlation_network.l <- generate_correlation_network(genus_fastspar_cor_filtered.m,
                                                            # relabeller_function = genus_relabeller_function,
                                                            # cor_threshold = 0.1)

length(as.character(unlist(lapply(rownames(otu_fastspar_cor_filtered.m), function(x) otu_taxonomy_map.df[x,]$taxonomy_species))))
unique(as.character(unlist(lapply(rownames(otu_fastspar_cor_filtered.m), function(x) otu_taxonomy_map.df[x,]$taxonomy_species))))

otu_correlation_network.l <- generate_correlation_network(cor_matrix = otu_fastspar_cor_filtered.m,
                                                            p_matrix = otu_fastspar_pval_filtered.m,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.4,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            label_colour = "black",
                                                            label_size = 2,
                                                            plot_height = 8,
                                                            plot_width = 8,
                                                            plot_title = "",
                                                            filename="Result_figures/correlation_analysis/networks/otu_filtered_fastspar_cor_network.pdf")
network_features.v <- (otu_correlation_network.l$network_data %>% activate(nodes) %>% as.data.frame())$name
otu_taxonomy_map.df[network_features.v,]$taxonomy_species

# genus_correlation_network.l$network_plot
# generate_correlation_network(genus_fastspar_cor_filtered.m)
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor_filtered.m,
                                                            p_matrix = genus_fastspar_pval_filtered.m,
                                                            relabeller_function = genus_relabeller_function,
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
                                                            exclude_to_from_df = edges_to_remove.df,
                                                            filename="Result_figures/correlation_analysis/networks/genus_filtered_fastspar_cor_network.pdf",
                                                            myseed = 1, edgetype = "link",show_p_label = F,file_type = "pdf")

genus_correlation_network.l$network_plot



# genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor_gold_star_filtered.m,
#                                                             p_matrix = genus_fastspar_pval_gold_star_filtered.m,
#                                                             p_value_threshold = 0.05,
#                                                             cor_threshold = 0.4,
#                                                             relabeller_function = genus_relabeller_function,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             label_colour = "black",
#                                                             label_size = 2,
#                                                             plot_height = 8,
#                                                             plot_width = 8,
#                                                             plot_title = "",
#                                                             filename="Result_figures/correlation_analysis/networks/genus_gold_star_filtered_fastspar_cor_network.pdf")
# genus_correlation_network.l$network_plot
# genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor_nose_filtered.m,
#                                                             p_matrix = genus_fastspar_pval_nose_filtered.m,
#                                                             p_value_threshold = 0.05,
#                                                             cor_threshold = 0.4,
#                                                             relabeller_function = genus_relabeller_function,
#                                                             node_size = 4,
#                                                             node_colour = "grey20",
#                                                             node_fill = "grey20",
#                                                             label_colour = "black",
#                                                             label_size = 2,
#                                                             plot_height = 8,
#                                                             plot_width = 8,
#                                                             plot_title = "",
#                                                             filename="Result_figures/correlation_analysis/networks/genus_nose_filtered_fastspar_cor_network.pdf")

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# Corrplots
# plot_corrplot(correlation_matrix = otu_fastspar_cor_filtered.m,
#               p_value_matrix = otu_fastspar_pval_filtered.m,
#               p_value_threshold = 0.05,
#               label_size = .3,
#               plot_title_size = 1,
#               # relabeller_function = genus_relabeller_function,
#               plot_title = "FastSpar correlations between all features,\n20% prevalence all samples, p-value < 0.05",
#               filename = "Result_other/correlation_analysis/corrplots/otu_filtered_fastspar_corrplot.pdf")

plot_corrplot(correlation_matrix = genus_fastspar_cor_filtered.m,
              p_value_matrix = genus_fastspar_pval_filtered.m,
              p_value_threshold = 0.05,
              label_size = .3,
              plot_title_size = 1,
              relabeller_function = genus_relabeller_function,
              plot_title = "FastSpar correlations between all genera,\n20% prevalence all samples, p-value < 0.05",
              filename = "Result_figures/correlation_analysis/corrplots/genus_filtered_fastspar_corrplot.pdf")

plot_corrplot(correlation_matrix = genus_fastspar_cor_nose_filtered.m,
              p_value_matrix = genus_fastspar_pval_nose_filtered.m,
              p_value_threshold = 0.05,
              label_size = .3,
              plot_title_size = 1,
              relabeller_function = genus_relabeller_function,
              plot_title = "FastSpar correlations between all genera,\n20% prevalence samples within each Nose group, p-value < 0.05",
              filename = "Result_figures/correlation_analysis/corrplots/genus_filtered_nose_fastspar_corrplot.pdf")

plot_corrplot(correlation_matrix = genus_fastspar_cor_gold_star_filtered.m,
              p_value_matrix = genus_fastspar_pval_gold_star_filtered.m,
              p_value_threshold = 0.05,
              label_size = .3,
              plot_title_size = 1,
              relabeller_function = genus_relabeller_function,
              plot_height = 10,
              plot_width = 10,
              plot_title = "FastSpar correlations between all genera,\n20% prevalence samples within each Gold Star group, p-value < 0.05",
              filename = "Result_figures/correlation_analysis/corrplots/genus_gold_star_filtered_fastspar_corrplot.pdf")

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# By feature
myfeature <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Cardiobacteriales;f__Cardiobacteriaceae;g__Dichelobacter"

source("code/helper_functions.R")
plot_feature_correlations_external(cor_matrix = genus_fastspar_cor_filtered.m,
                                   p_value_matrix = genus_fastspar_pval_filtered.m,
                                   feature = myfeature,
                                   top_n = 25,
                                   plot_width = 10, plot_height = 10)


for (feature in colnames(genus_fastspar_cor_filtered.m)){
  filename <- paste0("Result_figures/correlation_analysis/by_feature/",gsub(";", "_", genus_relabeller_function(feature)), ".pdf")
  plot_feature_correlations_external(cor_matrix = genus_fastspar_cor_filtered.m,
                                     p_value_matrix = genus_fastspar_pval_filtered.m,
                                     feature = feature,
                                     top_n = 25,
                                     filename = filename,
                                     format = "pdf",
                                     plot_width = 80, plot_height = 50,y_label_size = .6)
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# Correlation distributions, feature level

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

# spiec_genus <- spiec.easi(t(genus.m), method='mb', lambda.min.ratio=1e-2,
#                           nlambda=20, pulsar.params=list(rep.num=50))
# 
# 
# d <- ncol(t(genus.m))
# n <- nrow(t(genus.m))
# e <- d
# graph <- SpiecEasi::make_graph('cluster', d, e)
# 
# huge::huge.roc(spiec_genus$est$path, graph, verbose=FALSE)
# stars.pr(getOptMerge(spiec_genus), graph, verbose=FALSE)
# 
# spiec_genus <- spiec.easi(t(genus.m), method='mb', lambda.min.ratio=1e-2,
#                             nlambda=20, pulsar.params=list(rep.num=50)) 


  
