library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)


source("Code/helper_functions.R")

# Calculate the significance values for taxa between multiple groups
calculate_taxa_significances_multiple <- function(mydata, variable_column, value_column, taxonomy_column){
  results.df <- data.frame("Taxonomy" = character(),
                           "Variable" = character(),
                           "Group_1" = character(),
                           "Group_2" = character(),
                           "Dunn_pvalue" = character(),
                           "Dunn_padj" = character(),
                           "KrusW_pvalue" = character()
  )
  for (taxa in unique(mydata[,taxonomy_column])){ # For each taxa in the taxonomy column
    taxa_data <- subset(mydata, get(taxonomy_column) == taxa)
    n_groups = length(as.character(unique(taxa_data[,variable_column])))
    if (any(is.na(taxa_data[,variable_column]))){
      return()
    }
    if (n_groups > 2){
      kw <- kruskal.test(get(value_column)~get(variable_column), data = taxa_data)
      dunn <- dunnTest(x = get(value_column)~get(variable_column), data = taxa_data, method = "bh", alpha = 0.05)
      dunn <- separate(dunn$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
      names(dunn) <- c("Group_1","Group_2","Dunn_pvalue","Dunn_padj")
      dunn$Taxonomy <- taxa
      dunn$KrusW_pvalue <- kw$p.value
      dunn$Variable <- variable_column
      results.df <- rbind(results.df, dunn)
    }
  }
  results.df$P_value_label <- as.character(lapply(results.df[,"Dunn_padj"], function(x) ifelse(x <= 0.001, "***", 
                                                                                                 ifelse(x <= 0.01, "**", 
                                                                                                        ifelse(x <= 0.05, "*", "ns")))))
  
  results.df[,c("Taxonomy", "Variable", "Group_1","Group_2", "Dunn_pvalue", "Dunn_padj", "KrusW_pvalue")]
}

process_significances <- function(my_sig_data,my_abundance_data,variable_column,value_column){
  my_sig_data <- my_sig_data[which(my_sig_data[,"Dunn_padj"] <= sig_threshold),]
  
  # Determine the maximum diversity value for the pair of groups being compared
  for (row in 1:nrow(my_sig_data)){
    group_1 <- as.character(my_sig_data[row,"Group_1"])
    group_2 <- as.character(my_sig_data[row,"Group_2"])
    y_max <- max(my_abundance_data[which(my_abundance_data[,variable_column] == group_1),][,value_column],
                 my_abundance_data[which(my_abundance_data[,variable_column] == group_2),][,value_column])
    my_sig_data[row,"y_max"] <- y_max
    my_sig_data[row,"level_index_group_1"] <- which(levels(my_abundance_data[,variable_column]) == group_1)
    my_sig_data[row,"level_index_group_2"] <- which(levels(my_abundance_data[,variable_column]) == group_2)
    my_sig_data[row,"level_distance"] <- abs(my_sig_data[row,"level_index_group_2"] - my_sig_data[row,"level_index_group_1"])
  }
  
  my_sig_data <- my_sig_data[order(my_sig_data$level_distance),]
  scale <- sig_line_starting_scale # starting scale
  for (row in 1:nrow(my_sig_data)){
    my_sig_data[row, "y_position"] <- max(my_sig_data[, "y_max"]) * scale
    scale <- scale + sig_line_scaling_percentage # increase scale value
  }
  
  my_sig_data$Group_1 <- factor(my_sig_data$Group_1, levels = levels(my_abundance_data[,variable_column]))
  my_sig_data$Group_2 <- factor(my_sig_data$Group_2, levels = levels(my_abundance_data[,variable_column]))
  
  my_sig_data$P_value_label <- as.character(lapply(my_sig_data[,"Dunn_padj"], function(x) ifelse(x <= 0.001, "***", 
                                                                                                 ifelse(x <= 0.01, "**", 
                                                                                                        ifelse(x <= 0.05, "*", "ns")))))
  
  my_sig_data
}


# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Ensure rownames are the Index
rownames(metadata.df) <- metadata.df$Index

# Set ordering of variables
metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status, levels = c("Acute Otitis Media", "Perforation","Effusion", "HxOM","Never OM"))

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load abundance data
otu_rel.df <- read.csv("Result_tables/relative_abundance_tables/OTU_relative_abundances.csv", header = T)
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)


top_genus <- abundance_summary(genus_rel.df,20)$taxonomy_genus
genus_palette <- setNames(c(palette_20c[1:length(top_genus)], "grey"), c(top_genus, "Other"))

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# All samples

# Get top genus
top_genus_all_samples <- abundance_summary(genus_rel.df,20)$taxonomy_genus

# Filter to top genus
genus_rel_filtered.df <-genus_rel.df[genus_rel.df$taxonomy_genus %in% top_genus_all_samples,]

# Collapse abundance dataframe
genus_rel_filtered.df <- collapse_abundance_dataframe(genus_rel_filtered.df)

# Combine with metadata
genus_rel_filtered.df <- left_join(genus_rel_filtered.df, metadata.df, by = c("Sample"= "Index"))

genus_rel_filtered.df$Label <- as.character(lapply(genus_rel_filtered.df$taxonomy_genus,first_resolved_taxonomy))


otitis_status_genus_significances.df <- calculate_taxa_significances_multiple(mydata = genus_rel_filtered.df, 
                                                                              variable_column = "Otitis_Status",
                                                                              value_column = "Relative_abundance",
                                                                              taxonomy_column = "taxonomy_genus")
otitis_status_genus_significances.df <- subset(otitis_status_genus_significances.df, Dunn_padj <= 0.05)
otitis_status_genus_significances.df <- otitis_status_genus_significances.df[order(otitis_status_genus_significances.df$Taxonomy),]
write.csv(file = "Result_tables/abundance_analysis_tables/otitis_status_genus_significances.csv",
          x = otitis_status_genus_significances.df, row.names = F, quote = F)



colours.l <- unique(genus_rel_filtered.df[,c("Otitis_Status", "Otitis_Status_colour")])
colours.l <- setNames(as.character(colours.l[,"Otitis_Status_colour"]), 
                      colours.l[,"Otitis_Status"])

genus_rel_filtered.df <- genus_rel_filtered.df %>% filter(taxonomy_genus %in% top_genus_all_samples[1:9])


otitis_top_genus_boxplot <- 
  ggplot(genus_rel_filtered.df, aes(x = Label, 
                                  y = Relative_abundance*100,
                                  fill = Otitis_Status,
                                  shape = Otitis_Status)) +
  geom_boxplot(position = position_dodge(width =.75), 
               outlier.shape = NA, width=.5,lwd =.2) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = 0,
                                              dodge.width = .75)) +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = .8,
               position = position_dodge(width = .75),show.legend = F) +
  xlab("Taxonomy") +
  ylab("Relative abundance %") +
  scale_shape_manual(values = c(25,24,23,22,21),name = "Otitis status") +
  scale_fill_manual(values = colours.l, name = "Otitis status") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 10)) + 
  # facet_wrap(~Community) +
  theme(
    # axis.text.x = element_text(angle = 90, size = 6,hjust = 0, vjust = 0.5),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black"),#face = "italic"
    axis.title = element_text(size = 7,face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 5),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_blank()
  )
ggsave(filename = "Result_figures/abundance_analysis_plots/otitis_top_genus_boxplot.pdf",
       plot = otitis_top_genus_boxplot,
       width = 25,
       height = 7,
       units = "cm")



otitis_top_genus_mean_sd_plot <-  ggplot(genus_rel_filtered.df, aes(x = Label, 
                                  y = Relative_abundance*100,
                                  colour = Otitis_Status,
                                  fill = Otitis_Status,
                                  shape = Otitis_Status)) +
  
  # scale_y_log10() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),colour="black",
               geom="errorbar", width=0.5,lwd=.8, position = position_dodge(width = 0.75)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0.4,lwd=.3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),colour = "grey10",
               geom="point", position = position_dodge(width = 0.75)) +
  xlab("Taxonomy") +
  ylab("Relative abundance %") +
  scale_shape_manual(values = c(25,24,23,22,21),name = "Otitis status") +
  scale_fill_manual(values = colours.l, name = "Otitis status") +
  scale_colour_manual(values = colours.l, name = "Otitis status") +
  scale_y_continuous(breaks = seq(0,70, by = 10),expand = c(0,0)) +
  coord_cartesian(ylim = c(0,70)) +
  theme(
    # axis.text.x = element_text(angle = 90, size = 6,hjust = 0, vjust = 0.5),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(colour = "grey50",size = .3,linetype = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black"),#face = "italic"
    axis.title = element_text(size = 7,face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 5),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_blank()
  )

ggsave(filename = "Result_figures/abundance_analysis_plots/otitis_top_genus_mean_sd_plot.pdf",
       plot = otitis_top_genus_mean_sd_plot,
       width = 25,
       height = 7,
       units = "cm")

