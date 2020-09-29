# Heatmaps for differential abundance results

# Community, Gold_star, Tympanic_membrane and Nose


library(ComplexHeatmap)
library(reshape2)

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


source("code/helper_functions.R")



# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)
metadata.df <- metadata.df[order(metadata.df$Index),]
rownames(metadata.df) <- metadata.df$Index

# metadata.df <- subset(metadata.df, Tympanic_membrane != "Unable to visualise/ Not examined")

# Load abundance data
otu_rel.m <- as.matrix(read.table("Result_tables/relative_abundance_tables/OTU_relative_abundances.csv", sep =",", header =T, row.names = 1))
genus_rel.m <- as.matrix(read.table("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", sep =",", header =T, row.names = 1))

# Load the DESeq results
otu_deseq.df <- read.csv("Result_tables/DESeq_results/OTU_deseq.csv", header =T)
genus_deseq.df <- read.csv("Result_tables/DESeq_results/Genus_deseq.csv", header =T)

variables_of_interest <- c("Community", "Gold_Star","Tympanic_membrane", "Nose","Tympanic_membrane_Gold_Star", 
                           "Otitis_Status", "Otitis_Status_Gold_Star")

# Filter to variables of interest
otu_deseq.df <- subset(otu_deseq.df, Variable %in% variables_of_interest)
genus_deseq.df <- subset(genus_deseq.df, Variable %in% variables_of_interest)

# Get the genus/features that are significant
signficant_otu.v <- unique(otu_deseq.df$OTU)
signficant_otu_community.v <- unique(subset(otu_deseq.df, Variable == "Community")$OTU)
signficant_otu_gold_star.v <- unique(subset(otu_deseq.df, Variable == "Gold_Star")$OTU)
signficant_otu_tympanic_membrane.v <- unique(subset(otu_deseq.df, Variable == "Tympanic_membrane")$OTU)
signficant_otu_nose.v <- unique(subset(otu_deseq.df, Variable == "Nose")$OTU)
signficant_otu_otitis_status.v <- unique(subset(otu_deseq.df, Variable == "Otitis_Status")$OTU)
signficant_otu_tympanic_membrane_gold_star.v <- unique(subset(otu_deseq.df, Variable == "Tympanic_membrane_Gold_Star")$OTU)
signficant_otu_otitis_status_gold_star.v <- unique(subset(otu_deseq.df, Variable == "Otitis_Status_Gold_Star")$OTU)

signficant_genus.v <- unique(genus_deseq.df$Taxonomy)
signficant_genus_community.v <- unique(subset(genus_deseq.df, Variable == "Community")$Taxonomy)
signficant_genus_gold_star.v <- unique(subset(genus_deseq.df, Variable == "Gold_Star")$Taxonomy)
signficant_genus_tympanic_membrane.v <- unique(subset(genus_deseq.df, Variable == "Tympanic_membrane")$Taxonomy)
signficant_genus_nose.v <- unique(subset(genus_deseq.df, Variable == "Nose")$Taxonomy)
signficant_genus_otitis_status.v <- unique(subset(genus_deseq.df, Variable == "Otitis_Status")$Taxonomy)
signficant_genus_tympanic_membrane_gold_star.v <- unique(subset(genus_deseq.df, Variable == "Tympanic_membrane_Gold_Star")$Taxonomy)
signficant_genus_otitis_status_gold_star.v <- unique(subset(genus_deseq.df, Variable == "Otitis_Status_Gold_Star")$Taxonomy)

# Create abundance matrix for all significant features/taxa
otu_rel_all_sig.m <- otu_rel.m[signficant_otu.v,] * 100
otu_rel_community_sig.m <- otu_rel.m[signficant_otu_community.v,] * 100
otu_rel_gold_star_sig.m <- otu_rel.m[signficant_otu_gold_star.v,] * 100
otu_rel_tympanic_membrane_sig.m <- otu_rel.m[signficant_otu_tympanic_membrane.v,] * 100
otu_rel_nose_sig.m <- otu_rel.m[signficant_otu_nose.v,] * 100
otu_rel_otitis_status_sig.m <- otu_rel.m[signficant_otu_otitis_status.v,] * 100
otu_rel_tympanic_membrane_gold_star_sig.m <- otu_rel.m[signficant_otu_tympanic_membrane_gold_star.v,] * 100
otu_rel_otitis_status_gold_star.m <- otu_rel.m[signficant_otu_otitis_status_gold_star.v,] * 100

genus_rel_all_sig.m <- genus_rel.m[signficant_genus.v,] * 100
genus_rel_community_sig.m <- genus_rel.m[signficant_genus_community.v,] * 100
genus_rel_gold_star_sig.m <- genus_rel.m[signficant_genus_gold_star.v,] * 100
genus_rel_tympanic_membrane_sig.m <- genus_rel.m[signficant_genus_tympanic_membrane.v,] * 100
genus_rel_nose_sig.m <- genus_rel.m[signficant_genus_nose.v,] * 100
genus_rel_otitis_status_sig.m <- genus_rel.m[signficant_genus_otitis_status.v,] * 100
genus_rel_tympanic_membrane_gold_star_sig.m <- genus_rel.m[signficant_genus_tympanic_membrane_gold_star.v,] * 100
genus_rel_otitis_status_gold_star.m <- genus_rel.m[signficant_genus_otitis_status_gold_star.v,] * 100


# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#                                                           Make heatmaps

# All features
# new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m), paste0(otu_relabeller_function(rownames(otu_rel_all_sig.m)), "    ", rownames(otu_rel_all_sig.m)))
new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m), paste0(otu_relabeller_function(rownames(otu_rel_all_sig.m)), "    ", rownames(otu_rel_all_sig.m)))
# paste0(lapply(rownames(otu_rel_all_sig.m), combined_otu_labeller), "    ", rownames(otu_rel_all_sig.m))
make_heatmap(myheatmap_matrix = otu_rel_all_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/all_significant_features_heatmap.pdf"),
             
             variables = c("Community", "Gold_Star","Tympanic_membrane", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 15,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
             )


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Features - community
new_row_labels.df <- data.frame(rownames(otu_rel_community_sig.m), 
                                paste0(otu_relabeller_function(rownames(otu_rel_community_sig.m)), "    ", rownames(otu_rel_community_sig.m)))
make_heatmap(myheatmap_matrix = otu_rel_community_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/community_significant_features_heatmap.pdf"),
             variables = c("Community", "Gold_Star","Tympanic_membrane", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 5,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Features - Gold star
new_row_labels.df <- data.frame(rownames(otu_rel_gold_star_sig.m), 
                                paste0(otu_relabeller_function(rownames(otu_rel_gold_star_sig.m)), "    ", rownames(otu_rel_gold_star_sig.m)))
make_heatmap(myheatmap_matrix = otu_rel_gold_star_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/gold_star_significant_features_heatmap.pdf"),
             variables = c("Gold_Star", "Community","Tympanic_membrane", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 5,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Features - Tympanic_membrane
new_row_labels.df <- data.frame(rownames(otu_rel_tympanic_membrane_sig.m), 
                                paste0(otu_relabeller_function(rownames(otu_rel_tympanic_membrane_sig.m)), "    ", rownames(otu_rel_tympanic_membrane_sig.m)))
make_heatmap(myheatmap_matrix = otu_rel_tympanic_membrane_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_significant_features_heatmap.pdf"),
             variables = c("Tympanic_membrane","Community", "Gold_Star", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 7,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Features - Nose
new_row_labels.df <- data.frame(rownames(otu_rel_nose_sig.m), 
                                paste0(otu_relabeller_function(rownames(otu_rel_nose_sig.m)), "    ", 
                                       rownames(otu_rel_nose_sig.m)))
colnames(otu_rel_nose_sig.m) %in% rownames(metadata.df)
rownames(metadata.df) %in% colnames(otu_rel_nose_sig.m) 
make_heatmap(myheatmap_matrix = otu_rel_nose_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/nose_significant_features_heatmap.pdf"),
             variables = c("Nose","Community", "Tympanic_membrane", "Gold_Star"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 4,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Features - Otitis status
new_row_labels.df <- data.frame(rownames(otu_rel_otitis_status_sig.m), 
                                paste0(otu_relabeller_function(rownames(otu_rel_otitis_status_sig.m)), "    ", 
                                       rownames(otu_rel_otitis_status_sig.m)))
# colnames(otu_rel_otitis_status_sig.m) %in% rownames(metadata.df)
# metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status, levels = c("Acute Otitis Media","Effusion", "Perforation", "HxOM","Never OM"))
# temp <- otu_rel_otitis_status_sig.m
# temp[temp == 0] <- NA
make_heatmap(myheatmap_matrix = otu_rel_otitis_status_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/otitis_status_significant_features_heatmap.pdf"),
             variables = c("Otitis_Status","Community", "Nose", "Gold_Star"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 8 ,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Features - Tympanic_membrane_gold_Star
new_row_labels.df <- data.frame(rownames(otu_rel_tympanic_membrane_gold_star_sig.m), 
                                paste0(otu_relabeller_function(rownames(otu_rel_tympanic_membrane_gold_star_sig.m)), "    ", rownames(otu_rel_tympanic_membrane_gold_star_sig.m)))


make_heatmap(myheatmap_matrix = otu_rel_tympanic_membrane_gold_star_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_gold_star_significant_features_heatmap.pdf"),
             variables = c("Tympanic_membrane_Gold_Star","Community", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 13,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Features - Otitis status gold star
new_row_labels.df <- data.frame(rownames(otu_rel_otitis_status_gold_star.m), 
                                paste0(otu_relabeller_function(rownames(otu_rel_otitis_status_gold_star.m)), "    ", 
                                       rownames(otu_rel_otitis_status_gold_star.m)))
# colnames(otu_rel_otitis_status_sig.m) %in% rownames(metadata.df)
# metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status, levels = c("Acute Otitis Media","Effusion", "Perforation", "HxOM","Never OM"))
# temp <- otu_rel_otitis_status_sig.m
# temp[temp == 0] <- NA
make_heatmap(myheatmap_matrix = otu_rel_otitis_status_gold_star.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/otitis_status_gold_star_significant_features_heatmap.pdf"),
             variables = c("Otitis_Status_Gold_Star","Community", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 8 ,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# All genus
new_row_labels.df <- data.frame(rownames(genus_rel_all_sig.m), 
                                as.character(lapply(rownames(genus_rel_all_sig.m), first_resolved_taxonomy)))

make_heatmap(myheatmap_matrix = genus_rel_all_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/all_significant_genus_heatmap.pdf"),
             
             variables = c("Community", "Gold_Star","Tympanic_membrane", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 5,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             show_cell_values = F,
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - community
new_row_labels.df <- data.frame(rownames(genus_rel_community_sig.m), 
                                as.character(lapply(rownames(genus_rel_community_sig.m), first_resolved_taxonomy)))
make_heatmap(myheatmap_matrix = genus_rel_community_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/community_significant_genus_heatmap.pdf"),
             variables = c("Community", "Gold_Star","Tympanic_membrane", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 2.5,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - Gold star
new_row_labels.df <- data.frame(rownames(genus_rel_gold_star_sig.m), 
                                as.character(lapply(rownames(genus_rel_gold_star_sig.m), first_resolved_taxonomy)))
make_heatmap(myheatmap_matrix = genus_rel_gold_star_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/gold_star_significant_genus_heatmap.pdf"),
             variables = c("Gold_Star", "Community","Tympanic_membrane", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 3,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - Tympanic_membrane
new_row_labels.df <- data.frame(rownames(genus_rel_tympanic_membrane_sig.m), 
                                as.character(lapply(rownames(genus_rel_tympanic_membrane_sig.m), first_resolved_taxonomy)))
make_heatmap(myheatmap_matrix = genus_rel_tympanic_membrane_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_significant_genus_heatmap.pdf"),
             variables = c("Tympanic_membrane","Community", "Gold_Star", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 2.5,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - Nose
new_row_labels.df <- data.frame(rownames(genus_rel_nose_sig.m), 
                                as.character(lapply(rownames(genus_rel_nose_sig.m), first_resolved_taxonomy))) 
make_heatmap(myheatmap_matrix = genus_rel_nose_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/nose_significant_genus_heatmap.pdf"),
             variables = c("Nose","Community", "Tympanic_membrane", "Gold_Star"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 2.5,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - Otitis status
new_row_labels.df <- data.frame(rownames(genus_rel_otitis_status_sig.m), 
                                as.character(lapply(rownames(genus_rel_otitis_status_sig.m), first_resolved_taxonomy)))
make_heatmap(myheatmap_matrix = genus_rel_otitis_status_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/otitis_status_significant_genus_heatmap.pdf"),
             variables = c("Otitis_Status","Community", "Nose", "Gold_Star"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 3 ,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - Tympanic_membrane_gold_Star
new_row_labels.df <- data.frame(rownames(genus_rel_tympanic_membrane_gold_star_sig.m), 
                                as.character(lapply(rownames(genus_rel_tympanic_membrane_gold_star_sig.m), first_resolved_taxonomy)))

make_heatmap(myheatmap_matrix = genus_rel_tympanic_membrane_gold_star_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_gold_star_significant_genus_heatmap.pdf"),
             variables = c("Tympanic_membrane_Gold_Star","Community", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 4,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - Otitis status gold star
new_row_labels.df <- data.frame(rownames(genus_rel_otitis_status_gold_star.m), 
                                as.character(lapply(rownames(genus_rel_otitis_status_gold_star.m), first_resolved_taxonomy)))
make_heatmap(myheatmap_matrix = genus_rel_otitis_status_gold_star.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/otitis_status_gold_star_significant_genus_heatmap.pdf"),
             variables = c("Otitis_Status_Gold_Star","Community", "Nose"),
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height =2.8,
             plot_width = 20,
             simple_anno_size = unit(.5,"cm"),
             grid_colour = "grey",
             
             col_name_size = 10,
             # row_name_size = 10,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             palette_choice = 'purple',
             show_column_dend = F,
             show_row_dend = F,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 6,
)



