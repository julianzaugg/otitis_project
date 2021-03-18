# Heatmaps for differential abundance results


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
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)
metadata.df <- metadata.df[order(metadata.df$Index),]
rownames(metadata.df) <- metadata.df$Index

# Remove AOM, just make values NA
metadata.df[metadata.df$Otitis_Status == "Acute Otitis Media","Otitis_Status"] <- NA
metadata.df <- metadata.df[!is.na(metadata.df$Otitis_Status),]

# Set order of Otitis Status
metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status , levels = c("Never OM", "HxOM", "Effusion","Perforation"))
# Set order of Season
metadata.df$Season <- factor(metadata.df$Season, levels = c("Spring", "Winter", "Autumn"))
# Set order of Nose
metadata.df$Nose <- factor(metadata.df$Nose, levels = c("Normal", "Serous", "Purulent"))
# Set order of No_peop_res_discrete
metadata.df$No_peop_res_discrete <- factor(metadata.df$No_peop_res_discrete, levels = c("2 to 3", "4 to 6", "7 to 12", "Unknown"))

variables_of_interest <- c("Community", "Nose", "Otitis_Status", "Season", "No_peop_res_discrete")

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load abundance data
otu_rel.m <- as.matrix(read.table("Result_tables/relative_abundance_tables/OTU_relative_abundances.csv", sep =",", header =T, row.names = 1))
genus_rel.m <- as.matrix(read.table("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", sep =",", header =T, row.names = 1))

genus_rel.m <- genus_rel.m[,colnames(genus_rel.m) %in% metadata.df$Index]

# Load the DESeq results
otu_deseq.df <- read.csv("Result_tables/DESeq_results/OTU_deseq.csv", header =T)
genus_deseq.df <- read.csv("Result_tables/DESeq_results/Genus_deseq.csv", header =T)

otu_deseq_within_community.df <- read.csv("Result_tables/DESeq_results/OTU_deseq_within_community.csv", header =T)
genus_deseq_within_community.df <- read.csv("Result_tables/DESeq_results/Genus_deseq_within_community.csv", header =T)

# otu_within_community_otitis_status_deseq.df <- read.csv("Result_tables/DESeq_results/OTU_deseq_within_Community__Otitis_Status.csv", header =T)
# genus_within_community_otitis_status_deseq.df <- read.csv("Result_tables/DESeq_results/Genus_deseq_within_Community__Otitis_Status.csv", header =T)

# Filter to variables of interest
otu_deseq.df <- subset(otu_deseq.df, Variable %in% variables_of_interest)
genus_deseq.df <- subset(genus_deseq.df, Variable %in% variables_of_interest)

otu_deseq_within_community.df <- subset(otu_deseq_within_community.df, Variable %in% variables_of_interest)
genus_deseq_within_community.df <- subset(genus_deseq_within_community.df, Variable %in% variables_of_interest)

# otu_within_community_otitis_status_deseq.df <- subset(otu_within_community_otitis_status_deseq.df, Variable %in% variables_of_interest)
# genus_within_community_otitis_status_deseq.df <- subset(genus_within_community_otitis_status_deseq.df, Variable %in% variables_of_interest)

# Get the genus/features that are significant
significant_otu.v <- unique(otu_deseq.df$OTU)
significant_otu_community.v <- unique(subset(otu_deseq.df, Variable == "Community")$OTU)
significant_otu_gold_star.v <- unique(subset(otu_deseq.df, Variable == "Gold_Star")$OTU)
significant_otu_tympanic_membrane.v <- unique(subset(otu_deseq.df, Variable == "Tympanic_membrane")$OTU)
significant_otu_nose.v <- unique(subset(otu_deseq.df, Variable == "Nose")$OTU)
significant_otu_otitis_status.v <- unique(subset(otu_deseq.df, Variable == "Otitis_Status")$OTU)
significant_otu_tympanic_membrane_gold_star.v <- unique(subset(otu_deseq.df, Variable == "Tympanic_membrane__Gold_Star")$OTU)
significant_otu_otitis_status_gold_star.v <- unique(subset(otu_deseq.df, Variable == "Otitis_Status__Gold_Star")$OTU)
# significant_otu_within_community_otitis_status.v <- unique(otu_within_community_otitis_status_deseq.df$OTU)
significant_otu_within_community.v <- unique(otu_deseq_within_community.df$OTU)

significant_genus.v <- unique(genus_deseq.df$Taxonomy)
significant_genus_community.v <- unique(subset(genus_deseq.df, Variable == "Community")$Taxonomy) #*
significant_genus_gold_star.v <- unique(subset(genus_deseq.df, Variable == "Gold_Star")$Taxonomy)
significant_genus_tympanic_membrane.v <- unique(subset(genus_deseq.df, Variable == "Tympanic_membrane")$Taxonomy)
significant_genus_nose.v <- unique(subset(genus_deseq.df, Variable == "Nose")$Taxonomy) #*
significant_genus_otitis_status.v <- unique(subset(genus_deseq.df, Variable == "Otitis_Status")$Taxonomy)  #*
significant_genus_tympanic_membrane_gold_star.v <- unique(subset(genus_deseq.df, Variable == "Tympanic_membrane__Gold_Star")$Taxonomy)
significant_genus_otitis_status_gold_star.v <- unique(subset(genus_deseq.df, Variable == "Otitis_Status__Gold_Star")$Taxonomy)
# significant_genus_within_community_otitis_status.v <- unique(genus_within_community_otitis_status_deseq.df$Taxonomy)
significant_genus_within_community.v <- unique(genus_deseq_within_community.df$Taxonomy)

# Create abundance matrix for all significant features/taxa
otu_rel_all_sig.m <- otu_rel.m[significant_otu.v,] * 100
otu_rel_community_sig.m <- otu_rel.m[significant_otu_community.v,] * 100
otu_rel_gold_star_sig.m <- otu_rel.m[significant_otu_gold_star.v,] * 100
otu_rel_tympanic_membrane_sig.m <- otu_rel.m[significant_otu_tympanic_membrane.v,] * 100
otu_rel_nose_sig.m <- otu_rel.m[significant_otu_nose.v,] * 100
otu_rel_otitis_status_sig.m <- otu_rel.m[significant_otu_otitis_status.v,] * 100
otu_rel_tympanic_membrane_gold_star_sig.m <- otu_rel.m[significant_otu_tympanic_membrane_gold_star.v,] * 100
otu_rel_otitis_status_gold_star.m <- otu_rel.m[significant_otu_otitis_status_gold_star.v,] * 100
# otu_rel_within_community_otitis_status.m <- otu_rel.m[significant_otu_within_community_otitis_status.v,] * 100
otu_rel_within_community_all_sig.m <- otu_rel.m[significant_otu_within_community.v,] * 100

genus_rel_all_sig.m <- genus_rel.m[significant_genus.v,] * 100
genus_rel_community_sig.m <- genus_rel.m[significant_genus_community.v,] * 100
genus_rel_gold_star_sig.m <- genus_rel.m[significant_genus_gold_star.v,] * 100
genus_rel_tympanic_membrane_sig.m <- genus_rel.m[significant_genus_tympanic_membrane.v,] * 100
genus_rel_nose_sig.m <- genus_rel.m[significant_genus_nose.v,] * 100
genus_rel_otitis_status_sig.m <- genus_rel.m[significant_genus_otitis_status.v,] * 100
genus_rel_tympanic_membrane_gold_star_sig.m <- genus_rel.m[significant_genus_tympanic_membrane_gold_star.v,] * 100
genus_rel_otitis_status_gold_star.m <- genus_rel.m[significant_genus_otitis_status_gold_star.v,] * 100
# genus_rel_within_community_otitis_status.m <- genus_rel.m[significant_genus_within_community_otitis_status.v,] * 100
genus_rel_within_community_all_sig.m <- genus_rel.m[significant_genus_within_community.v,] * 100

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#                                                           Make heatmaps

# All genus
new_row_labels.df <- data.frame(rownames(genus_rel_all_sig.m), 
                                as.character(lapply(rownames(genus_rel_all_sig.m), first_resolved_taxonomy)))

# temp <- genus_rel_all_sig.m
# temp[temp > 0] <- 1
my_heatmap <- make_heatmap(myheatmap_matrix = genus_rel_all_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/all_significant_genus_heatmap.pdf"),
             
             variables = c("Community","Otitis_Status","Nose", "Season", "No_peop_res_discrete"),
             row_title = "Genus (lowest resolved taxonomy)",
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 3,
             plot_width = 15,
             simple_anno_size = unit(.3,"cm"),
             grid_colour = "grey60",
             grid_thickness = .5,
             show_cell_values = F,
             
             col_name_size = 6,
             row_name_size = 8,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 8,
             col_annotation_label_size = 8,
             col_annotation_title_size = 8,
             col_annotation_legend_grid_height = .4,
             col_annotation_legend_grid_width = .4,
             scale_legend_title_size = 8,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             # legend_labels = c(c(0, 0.001, 0.01,0.05, seq(.1,.5,.1))*100, "> 60"),
             # my_breaks = c(0, 0.001, 0.01,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             # my_palette = colorRampPalette(c("#17468a","#ffdd47","#99113a"))(20),
             palette_choice = 'red',
             show_column_dend = F,
             show_row_dend = F,
             
)

# Heatmap(row_names_gp = gpar(colour = "red"),)
svg(filename = "Result_figures/DESeq_plots/all_significant_genus_heatmap.svg",
    width = 15,
    height = 3)
draw(my_heatmap$heatmap, annotation_legend_list = c(my_heatmap$legend),merge_legends =T)
dev.off()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# All Genus Community only 
new_row_labels.df <- data.frame(rownames(genus_rel_community_sig.m), 
                                as.character(lapply(rownames(genus_rel_community_sig.m), first_resolved_taxonomy)))

my_heatmap <- make_heatmap(myheatmap_matrix = genus_rel_community_sig.m, 
                           mymetadata = metadata.df,
                           filename = paste0("Result_figures/DESeq_plots/all_significant_genus_community_heatmap.pdf"),
                           
                           variables = c("Community","Otitis_Status","Nose", "Season", "No_peop_res_discrete"),
                           row_title = "Genus (lowest resolved taxonomy)",
                           my_row_labels = new_row_labels.df,
                           cluster_columns = F,
                           cluster_rows = F,
                           plot_height = 2,
                           plot_width = 15,
                           simple_anno_size = unit(.3,"cm"),
                           grid_colour = "grey60",
                           grid_thickness = .5,
                           show_cell_values = F,
                           
                           col_name_size = 6,
                           row_name_size = 8,
                           column_title_size = 10,
                           row_title_size = 10,
                           annotation_bar_name_size = 8,
                           col_annotation_label_size = 8,
                           col_annotation_title_size = 8,
                           col_annotation_legend_grid_height = .4,
                           col_annotation_legend_grid_width = .4,
                           scale_legend_title_size = 8,
                           
                           legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
                           my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
                           discrete_legend = T,
                           legend_title = "Relative abundance %",
                           # my_palette = c("white","red"),
                           # my_palette = colorRampPalette(c("#17468a","#ffdd47","#99113a"))(20),
                           palette_choice = 'red',
                           show_column_dend = T,
                           show_row_dend = F,
)
svg(filename = "Result_figures/DESeq_plots/all_significant_genus_community_heatmap.svg",
    width = 15,
    height = 4)
draw(my_heatmap$heatmap, annotation_legend_list = c(my_heatmap$legend),merge_legends =T)
dev.off()
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# All Genus Nose only 
new_row_labels.df <- data.frame(rownames(genus_rel_nose_sig.m), 
                                as.character(lapply(rownames(genus_rel_nose_sig.m), first_resolved_taxonomy)))

my_heatmap <- make_heatmap(myheatmap_matrix = genus_rel_nose_sig.m, 
                           mymetadata = metadata.df,
                           filename = paste0("Result_figures/DESeq_plots/all_significant_genus_nose_heatmap.pdf"),
                           
                           variables = c("Nose","Community","Otitis_Status","Season", "No_peop_res_discrete"),
                           row_title = "Genus (lowest resolved taxonomy)",
                           my_row_labels = new_row_labels.df,
                           cluster_columns = F,
                           cluster_rows = F,
                           plot_height = 2,
                           plot_width = 15,
                           simple_anno_size = unit(.3,"cm"),
                           grid_colour = "grey60",
                           grid_thickness = .5,
                           show_cell_values = F,
                           
                           col_name_size = 6,
                           row_name_size = 8,
                           column_title_size = 10,
                           row_title_size = 10,
                           annotation_bar_name_size = 8,
                           col_annotation_label_size = 8,
                           col_annotation_title_size = 8,
                           col_annotation_legend_grid_height = .4,
                           col_annotation_legend_grid_width = .4,
                           scale_legend_title_size = 8,
                           
                           legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
                           my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
                           discrete_legend = T,
                           legend_title = "Relative abundance %",
                           # my_palette = c("white","red"),
                           # my_palette = colorRampPalette(c("#17468a","#ffdd47","#99113a"))(20),
                           palette_choice = 'red',
                           show_column_dend = T,
                           show_row_dend = F,
)
svg(filename = "Result_figures/DESeq_plots/all_significant_genus_nose_heatmap.svg",
    width = 15,
    height = 2)
draw(my_heatmap$heatmap, annotation_legend_list = c(my_heatmap$legend),merge_legends =T)
dev.off()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# All Genus within community
new_row_labels.df <- data.frame(rownames(genus_rel_within_community_all_sig.m), 
                                as.character(lapply(rownames(genus_rel_within_community_all_sig.m), first_resolved_taxonomy)))

my_heatmap <- make_heatmap(myheatmap_matrix = genus_rel_within_community_all_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/all_significant_genus_within_community_heatmap.pdf"),
             
             variables = c("Community","Otitis_Status","Nose", "Season", "No_peop_res_discrete"),
             row_title = "Genus (lowest resolved taxonomy)",
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 4,
             plot_width = 15,
             simple_anno_size = unit(.3,"cm"),
             grid_colour = "grey60",
             grid_thickness = .5,
             show_cell_values = F,
             
             col_name_size = 6,
             row_name_size = 8,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 8,
             col_annotation_label_size = 8,
             col_annotation_title_size = 8,
             col_annotation_legend_grid_height = .4,
             col_annotation_legend_grid_width = .4,
             scale_legend_title_size = 8,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             # my_palette = colorRampPalette(c("#17468a","#ffdd47","#99113a"))(20),
             palette_choice = 'red',
             show_column_dend = T,
             show_row_dend = F,
)
svg(filename = "Result_figures/DESeq_plots/all_significant_genus_within_community_heatmap.svg",
    width = 15,
    height = 4)
draw(my_heatmap$heatmap, annotation_legend_list = c(my_heatmap$legend),merge_legends =T)
dev.off()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# Features - community
# new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m),
                                # paste0(otu_relabeller_function(rownames(otu_rel_all_sig.m)), "    ", rownames(otu_rel_community_sig.m)))

new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m), 
                                paste0(unlist(lapply(rownames(otu_rel_all_sig.m),combined_otu_labeller)), 
                                       "  ",
                                       unlist(lapply(rownames(otu_rel_all_sig.m), function(x) otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,"ASV_ID"]))
                                       )
                                )
# new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m), 
                                # unlist(lapply(rownames(otu_rel_all_sig.m),combined_otu_labeller)))

my_heatmap <- make_heatmap(myheatmap_matrix = otu_rel_all_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/all_significant_features_heatmap.pdf"),
             variables = c("Community","Otitis_Status","Nose", "Season", "No_peop_res_discrete"),
             row_title = "ASV (lowest resolved taxonomy)",
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 8,
             plot_width = 15,
             simple_anno_size = unit(.3,"cm"),
             grid_colour = "grey60",
             grid_thickness = .5,
             show_cell_values = F,
             
             col_name_size = 6,
             row_name_size = 8,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 8,
             col_annotation_label_size = 8,
             col_annotation_title_size = 8,
             col_annotation_legend_grid_height = .4,
             col_annotation_legend_grid_width = .4,
             scale_legend_title_size = 8,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             # my_palette = colorRampPalette(c("#17468a","#ffdd47","#99113a"))(20),
             palette_choice = 'red',
             show_column_dend = T,
             show_row_dend = F,
)
svg(filename = "Result_figures/DESeq_plots/all_significant_features_heatmap.svg",
    width = 15,
    height = 8)
draw(my_heatmap$heatmap, annotation_legend_list = c(my_heatmap$legend),merge_legends =T)
dev.off()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# Features within community
# new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m),
# paste0(otu_relabeller_function(rownames(otu_rel_all_sig.m)), "    ", rownames(otu_rel_community_sig.m)))
# new_row_labels.df <- data.frame(rownames(otu_rel_within_community_all_sig.m), 
                                # unlist(lapply(rownames(otu_rel_within_community_all_sig.m),combined_otu_labeller)))

new_row_labels.df <- data.frame(rownames(otu_rel_within_community_all_sig.m), 
                                paste0(unlist(lapply(rownames(otu_rel_within_community_all_sig.m),combined_otu_labeller)), 
                                       "  ",
                                       unlist(lapply(rownames(otu_rel_within_community_all_sig.m), function(x) otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,"ASV_ID"]))
                                )
)

my_heatmap <- make_heatmap(myheatmap_matrix = otu_rel_within_community_all_sig.m, 
             mymetadata = metadata.df,
             filename = paste0("Result_figures/DESeq_plots/all_significant_features_within_community_heatmap.pdf"),
             variables = c("Community","Otitis_Status","Nose", "Season", "No_peop_res_discrete"),
             row_title = "ASV (lowest resolved taxonomy)",
             my_row_labels = new_row_labels.df,
             cluster_columns = F,
             cluster_rows = F,
             plot_height = 15,
             plot_width = 15,
             simple_anno_size = unit(.3,"cm"),
             grid_colour = "grey60",
             grid_thickness = .5,
             show_cell_values = F,
             
             col_name_size = 6,
             row_name_size = 8,
             column_title_size = 10,
             row_title_size = 10,
             annotation_bar_name_size = 8,
             col_annotation_label_size = 8,
             col_annotation_title_size = 8,
             col_annotation_legend_grid_height = .4,
             col_annotation_legend_grid_width = .4,
             scale_legend_title_size = 8,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Relative abundance %",
             # my_palette = c("white","red"),
             # my_palette = colorRampPalette(c("#17468a","#ffdd47","#99113a"))(20),
             palette_choice = 'red',
             show_column_dend = T,
             show_row_dend = F,
)

svg(filename = "Result_figures/DESeq_plots/all_significant_features_within_community_heatmap.svg",
    width = 15,
    height = 15)
draw(my_heatmap$heatmap, annotation_legend_list = c(my_heatmap$legend),merge_legends =T)
dev.off()











# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Genus - community
# new_row_labels.df <- data.frame(rownames(genus_rel_community_sig.m), 
#                                 as.character(lapply(rownames(genus_rel_community_sig.m), first_resolved_taxonomy)))
# 
# 
# 
# make_heatmap(myheatmap_matrix = genus_rel_community_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/community_significant_genus_heatmap.pdf"),
#              variables = c("Community"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 2,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Genus - Gold star
# new_row_labels.df <- data.frame(rownames(genus_rel_gold_star_sig.m), 
#                                 as.character(lapply(rownames(genus_rel_gold_star_sig.m), first_resolved_taxonomy)))
# make_heatmap(myheatmap_matrix = genus_rel_gold_star_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/gold_star_significant_genus_heatmap.pdf"),
#              variables = c("Gold_Star"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 3,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Genus - Tympanic_membrane
# new_row_labels.df <- data.frame(rownames(genus_rel_tympanic_membrane_sig.m), 
#                                 as.character(lapply(rownames(genus_rel_tympanic_membrane_sig.m), first_resolved_taxonomy)))
# make_heatmap(myheatmap_matrix = genus_rel_tympanic_membrane_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_significant_genus_heatmap.pdf"),
#              variables = c("Tympanic_membrane"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 2,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Genus - Nose
# new_row_labels.df <- data.frame(rownames(genus_rel_nose_sig.m), 
#                                 as.character(lapply(rownames(genus_rel_nose_sig.m), first_resolved_taxonomy))) 
# make_heatmap(myheatmap_matrix = genus_rel_nose_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/nose_significant_genus_heatmap.pdf"),
#              variables = c("Nose"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 2,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Genus - Otitis status
# new_row_labels.df <- data.frame(rownames(genus_rel_otitis_status_sig.m), 
#                                 as.character(lapply(rownames(genus_rel_otitis_status_sig.m), first_resolved_taxonomy)))
# make_heatmap(myheatmap_matrix = genus_rel_otitis_status_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/otitis_status_significant_genus_heatmap.pdf"),
#              variables = c("Otitis_Status"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 2.5,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Genus - Tympanic_membrane_gold_Star
# new_row_labels.df <- data.frame(rownames(genus_rel_tympanic_membrane_gold_star_sig.m), 
#                                 as.character(lapply(rownames(genus_rel_tympanic_membrane_gold_star_sig.m), first_resolved_taxonomy)))
# 
# make_heatmap(myheatmap_matrix = genus_rel_tympanic_membrane_gold_star_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_gold_star_significant_genus_heatmap.pdf"),
#              variables = c("Tympanic_membrane__Gold_Star"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 4,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Genus - Otitis status gold star
# new_row_labels.df <- data.frame(rownames(genus_rel_otitis_status_gold_star.m), 
#                                 as.character(lapply(rownames(genus_rel_otitis_status_gold_star.m), first_resolved_taxonomy)))
# make_heatmap(myheatmap_matrix = genus_rel_otitis_status_gold_star.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/otitis_status_gold_star_significant_genus_heatmap.pdf"),
#              variables = c("Otitis_Status__Gold_Star"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height =2.5,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Genus - Within community otitis status
# new_row_labels.df <- data.frame(rownames(genus_rel_within_community_otitis_status.m), 
#                                 as.character(lapply(rownames(genus_rel_within_community_otitis_status.m), first_resolved_taxonomy)))
# make_heatmap(myheatmap_matrix = genus_rel_within_community_otitis_status.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/within_community_otitis_status_significant_genus_heatmap.pdf"),
#              variables = c("Community__Otitis_Status", "Nose"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height =2.5,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# 
# 
# 
# 
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# 
# # All features
# # new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m), paste0(otu_relabeller_function(rownames(otu_rel_all_sig.m)), "    ", rownames(otu_rel_all_sig.m)))
# new_row_labels.df <- data.frame(rownames(otu_rel_all_sig.m), paste0(otu_relabeller_function(rownames(otu_rel_all_sig.m)), "    ", rownames(otu_rel_all_sig.m)))
# # paste0(lapply(rownames(otu_rel_all_sig.m), combined_otu_labeller), "    ", rownames(otu_rel_all_sig.m))
# make_heatmap(myheatmap_matrix = otu_rel_all_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/all_significant_features_heatmap.pdf"),
#              
#              variables = c("Community", "Gold_Star","Tympanic_membrane","Otitis_Status", "Nose"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 15,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
#              )
# 
# 
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - community
# new_row_labels.df <- data.frame(rownames(otu_rel_community_sig.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_community_sig.m)), "    ", rownames(otu_rel_community_sig.m)))
# make_heatmap(myheatmap_matrix = otu_rel_community_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/community_significant_features_heatmap.pdf"),
#              variables = c("Community"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 5,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - Gold star
# new_row_labels.df <- data.frame(rownames(otu_rel_gold_star_sig.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_gold_star_sig.m)), "    ", rownames(otu_rel_gold_star_sig.m)))
# make_heatmap(myheatmap_matrix = otu_rel_gold_star_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/gold_star_significant_features_heatmap.pdf"),
#              variables = c("Gold_Star"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 5,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - Tympanic_membrane
# new_row_labels.df <- data.frame(rownames(otu_rel_tympanic_membrane_sig.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_tympanic_membrane_sig.m)), "    ", rownames(otu_rel_tympanic_membrane_sig.m)))
# make_heatmap(myheatmap_matrix = otu_rel_tympanic_membrane_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_significant_features_heatmap.pdf"),
#              variables = c("Tympanic_membrane"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 7,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - Nose
# new_row_labels.df <- data.frame(rownames(otu_rel_nose_sig.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_nose_sig.m)), "    ", 
#                                        rownames(otu_rel_nose_sig.m)))
# colnames(otu_rel_nose_sig.m) %in% rownames(metadata.df)
# rownames(metadata.df) %in% colnames(otu_rel_nose_sig.m) 
# make_heatmap(myheatmap_matrix = otu_rel_nose_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/nose_significant_features_heatmap.pdf"),
#              variables = c("Nose"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 4,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - Otitis status
# new_row_labels.df <- data.frame(rownames(otu_rel_otitis_status_sig.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_otitis_status_sig.m)), "    ", 
#                                        rownames(otu_rel_otitis_status_sig.m)))
# # colnames(otu_rel_otitis_status_sig.m) %in% rownames(metadata.df)
# # metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status, levels = c("Acute Otitis Media","Effusion", "Perforation", "HxOM","Never OM"))
# # temp <- otu_rel_otitis_status_sig.m
# # temp[temp == 0] <- NA
# make_heatmap(myheatmap_matrix = otu_rel_otitis_status_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/otitis_status_significant_features_heatmap.pdf"),
#              variables = c("Otitis_Status"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 8 ,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - Tympanic_membrane_gold_Star
# new_row_labels.df <- data.frame(rownames(otu_rel_tympanic_membrane_gold_star_sig.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_tympanic_membrane_gold_star_sig.m)), "    ", rownames(otu_rel_tympanic_membrane_gold_star_sig.m)))
# 
# 
# make_heatmap(myheatmap_matrix = otu_rel_tympanic_membrane_gold_star_sig.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/tympanic_membrane_gold_star_significant_features_heatmap.pdf"),
#              variables = c("Tympanic_membrane__Gold_Star"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 13,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - Otitis status gold star
# new_row_labels.df <- data.frame(rownames(otu_rel_otitis_status_gold_star.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_otitis_status_gold_star.m)), "    ", 
#                                        rownames(otu_rel_otitis_status_gold_star.m)))
# # colnames(otu_rel_otitis_status_sig.m) %in% rownames(metadata.df)
# # metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status, levels = c("Acute Otitis Media","Effusion", "Perforation", "HxOM","Never OM"))
# # temp <- otu_rel_otitis_status_sig.m
# # temp[temp == 0] <- NA
# make_heatmap(myheatmap_matrix = otu_rel_otitis_status_gold_star.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/otitis_status_gold_star_significant_features_heatmap.pdf"),
#              variables = c("Otitis_Status__Gold_Star"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height = 8 ,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )
# 
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# # Features - Within community otitis status
# new_row_labels.df <- data.frame(rownames(otu_rel_within_community_otitis_status.m), 
#                                 paste0(otu_relabeller_function(rownames(otu_rel_within_community_otitis_status.m)), "    ", 
#                                        rownames(otu_rel_within_community_otitis_status.m)))
# 
# make_heatmap(myheatmap_matrix = otu_rel_within_community_otitis_status.m, 
#              mymetadata = metadata.df,
#              filename = paste0("Result_figures/DESeq_plots/within_community_otitis_status_significant_features_heatmap.pdf"),
#              variables = c("Community__Otitis_Status","Tympanic_membrane__Gold_Star", "Nose"),
#              my_row_labels = new_row_labels.df,
#              cluster_columns = F,
#              cluster_rows = F,
#              plot_height =15,
#              plot_width = 20,
#              simple_anno_size = unit(.5,"cm"),
#              grid_colour = "grey",
#              
#              col_name_size = 10,
#              # row_name_size = 10,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
#              my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
#              discrete_legend = T,
#              legend_title = "Relative abundance %",
#              # my_palette = c("white","red"),
#              palette_choice = 'purple',
#              show_column_dend = F,
#              show_row_dend = F,
#              column_title_size = 10,
#              row_title_size = 10,
#              annotation_bar_name_size = 6,
# )



