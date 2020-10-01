# Diversity calculations Alpha (Shannon, Chao1, Simpson) and Beta (CLR transformed counts) diveristy for each sample.
# Significance tests of each discrete group comparing diversity indices.
# Generate boxplots for discrete 
# Comparison of continuous variables vs diversity indices

# "Alpha diversity measures the diversity within a single sample and is generally 
# based on the number and relative abundance of taxa at some rank (e.g. species or OTUs).
# Beta diversity also uses the number of relative abundance of taxa at some rank, but measures 
# variation between samples. In other words, an alpha diversity statistic describes a single 
# sample and a beta diversity statistic describes how two samples compare."

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()
library(vegan)
library(reshape2)
library(ggplot2)
library(dplyr)

library(FSA)
library(phyloseq)
# library(nlme)

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


generate_p_labels <- function(sig_table){
  for (sig_column in c("Chao1_Dunn_padj", "Shannon_Dunn_padj", "Simpson_Dunn_padj")){
    metric = strsplit(sig_column, "_")[[1]][1]
    sig_table[,paste0(metric, "_p_label")] <-
      as.character(lapply(sig_table[,sig_column], 
                          function(x) ifelse(x <= 0.001, "***", 
                                             ifelse(x <= 0.01, "**",
                                                    ifelse(x <= 0.05, "*", "ns")))))
  }
  sig_table
}


setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")



# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Set the Index to be the row name
rownames(metadata.df) <- metadata.df$Index

# Define the discrete variables
discrete_variables <- c("Nose","Tympanic_membrane", "Otitis_Status",
                        "Season","Community","Gold_Star",
                        "H.influenzae_culture","M.catarrhalis_culture","S.pneumoniae_culture",
                        "Otitis_Status__Gold_Star", "Tympanic_membrane__Gold_Star",
                        "Community__Season","Community__Gold_Star","Community__Otitis_Status",
                        "H.Influenzae_qPCR", "M.catarrhalis_qPCR", "S.pneumoniae_qPCR",
                        "Corynebacterium_pseudodiphtheriticum","Dolosigranulum_pigrum","N_HRV")
                        # "N_Adeno","N_WUKI","N_BOCA","N_COV_OC43","N_COV_NL63",
                        # "N_HKU_1","N_ENT","N_hMPV","N_PARA_1","N_PARA_2","N_RSV_A","N_RSV_B","N_HRV","N_FLU_B","N_FLU_A","Virus_any")

metadata.df$Tympanic_membrane[metadata.df$Tympanic_membrane == "Unable to visualise/Not examined"] <- NA
# Combined with Community
for (x in discrete_variables){
  if (!grepl("Community", x) & !grepl("__", x)){
    metadata.df[,paste0(x,'__Community')] <-  with(metadata.df, paste0(get(x), "__",Community))
  }
}
discrete_variables <- c(discrete_variables, grep(".*__Community", colnames(metadata.df), value =T))

# Load the counts
otu.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts.csv", header =T, row.names = 1))
genus.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts.csv", header =T, row.names = 1))

# Order the matrices and metadata to be the same order
otu.m <- otu.m[,rownames(metadata.df)]
genus.m <- genus.m[,rownames(metadata.df)]

# Create the rarefied matrix
otu_rare.m <- t(rrarefy(t(otu.m[,colSums(otu.m) >= 2000]), 2000))
genus_rare.m <- t(rrarefy(t(genus.m[,colSums(genus.m) >= 2000]), 2000))

# Create phyloseq object
otu_rare_phyloseq <- otu_table(otu_rare.m, taxa_are_rows=TRUE)
genus_rare_phyloseq <- otu_table(genus_rare.m, taxa_are_rows=TRUE)

# Estimate alpha diversities
otu_rare_alpha.df <- estimate_richness(otu_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_rare_alpha.df <- otu_rare_alpha.df[rownames(metadata.df),]

genus_rare_alpha.df <- estimate_richness(genus_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
genus_rare_alpha.df <- genus_rare_alpha.df[rownames(metadata.df),]

# ---------------------------
# Combine with metadata
otu_rare_alpha.df <- left_join(metadata.df[c("Index", 
                                             discrete_variables, 
                                             grep("_colour", names(metadata.df), value = T))],m2df(otu_rare_alpha.df, "Index"), by = "Index")
genus_rare_alpha.df <- left_join(metadata.df[c("Index",
                                                   discrete_variables, 
                                                   grep("_colour", names(metadata.df), value = T))],m2df(genus_rare_alpha.df, "Index"), by = "Index")

# Remove colour columns
otu_rare_alpha.df <- otu_rare_alpha.df[,!grepl("_colour", names(otu_rare_alpha.df))]
genus_rare_alpha.df <- genus_rare_alpha.df[,!grepl("_colour", names(genus_rare_alpha.df))]
# ---------------------------

# Calculate summary for each variable
otu_alpha_diversity_summary.df <- summarise_diversities_each_variable(otu_rare_alpha.df, 
                                                                        variables = discrete_variables)

genus_alpha_diversity_summary.df <- summarise_diversities_each_variable(genus_rare_alpha.df, 
                                                                        variables = discrete_variables)


# Calculate summary for each variable for each COMMUNITY
# genus_alpha_diversity_summary_per_community.df <- data.frame()
# for (community in as.character(unique(metadata.df$Community))){
#   discrete_variables_subset <- discrete_variables[!discrete_variables %in% c("Community")]
#   genus_rare_alpha_subset.df <- subset(genus_rare_alpha.df, Community == "Remote")
#   temp <- summarise_diversities_each_variable(genus_rare_alpha_subset.df, 
#                                               variables = discrete_variables_subset)
#   temp$Community <- community
#   genus_alpha_diversity_summary_per_community.df <- rbind(genus_alpha_diversity_summary_per_community.df,temp)
# }

# Calculate significances for each variable (and for each COMMUNITY)
otu_alpha_mann_significances.df <- data.frame()
otu_alpha_dunn_significances_multiple.df <- data.frame()

genus_alpha_mann_significances.df <- data.frame()
genus_alpha_dunn_significances_multiple.df <- data.frame()
# genus_alpha_dunn_significances_per_community.df <- data.frame()
# genus_alpha_dunn_significances_multiple_per_community.df <- data.frame()
for (variable in discrete_variables){
  n_groups <- unique(genus_rare_alpha.df[[variable]])
  if (length(n_groups) < 3){
    print(paste0("Pair: ", variable))
    otu_alpha_mann_significances.df <- rbind(otu_alpha_mann_significances.df, 
                                               calculate_alpha_diversity_significance(otu_rare_alpha.df,
                                                                                      variable = variable))
    
    genus_alpha_mann_significances.df <- rbind(genus_alpha_mann_significances.df, 
                                               calculate_alpha_diversity_significance(genus_rare_alpha.df,
                                                                                      variable = variable))
  } else{
    print(paste0("Multiple: ", variable))
    otu_alpha_dunn_significances_multiple.df <- rbind(otu_alpha_dunn_significances_multiple.df, 
                                                        calculate_alpha_diversity_significance_multiple(otu_rare_alpha.df,
                                                                                                        variable))
    genus_alpha_dunn_significances_multiple.df <- rbind(genus_alpha_dunn_significances_multiple.df, 
                                                        calculate_alpha_diversity_significance_multiple(genus_rare_alpha.df,
                                                                                                        variable))
  }
}
# Remove entries that are not (near) significant
threshold <- 0.1
otu_alpha_mann_significances.df <- otu_alpha_mann_significances.df[apply(otu_alpha_mann_significances.df[,c("Shannon_MannW_padj", "Simpson_MannW_padj", "Chao1_MannW_padj")],1,min) <= threshold,]
genus_alpha_mann_significances.df <- genus_alpha_mann_significances.df[apply(genus_alpha_mann_significances.df[,c("Shannon_MannW_padj", "Simpson_MannW_padj", "Chao1_MannW_padj")],1,min) <= threshold,]

otu_alpha_dunn_significances_multiple.df <- otu_alpha_dunn_significances_multiple.df[apply(otu_alpha_dunn_significances_multiple.df[,c("Shannon_Dunn_padj", "Simpson_Dunn_padj", "Chao1_Dunn_padj")], 1, min) <= threshold,]
genus_alpha_dunn_significances_multiple.df <- genus_alpha_dunn_significances_multiple.df[apply(genus_alpha_dunn_significances_multiple.df[,c("Shannon_Dunn_padj", "Simpson_Dunn_padj", "Chao1_Dunn_padj")], 1, min) <= threshold,]

# Write per-sample diversities to file
write.csv(otu_rare_alpha.df,
          "Result_tables/diversity_analysis/otu/sample_otu_alpha_diversities.csv", quote = F, row.names = F)
write.csv(genus_rare_alpha.df,
          "Result_tables/diversity_analysis/genus/sample_genus_alpha_diversities.csv", quote = F, row.names = F
)
# Write variable summaries to file 
write.csv(otu_alpha_diversity_summary.df,
          "Result_tables/diversity_analysis/otu/otu_alpha_diversities_summary.csv", quote = F, row.names = F
)
write.csv(genus_alpha_diversity_summary.df,
          "Result_tables/diversity_analysis/genus/genus_alpha_diversities_summary.csv", quote = F, row.names = F
)
# Write pair significances
write.csv(otu_alpha_mann_significances.df,
          "Result_tables/diversity_analysis/otu/otu_alpha_diversities_signficances_mannwhitney.csv", quote = F, row.names = F
)
write.csv(genus_alpha_mann_significances.df,
          "Result_tables/diversity_analysis/genus/genus_alpha_diversities_signficances_mannwhitney.csv", quote = F, row.names = F
)
# Write Dunn significances
write.csv(otu_alpha_dunn_significances_multiple.df,
          "Result_tables/diversity_analysis/otu/otu_alpha_diversities_signficances_dunn.csv", quote = F, row.names = F
)
write.csv(genus_alpha_dunn_significances_multiple.df,
          "Result_tables/diversity_analysis/genus/genus_alpha_diversities_signficances_dunn.csv", quote = F, row.names = F
)




# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
# Calculate summaries and signifiance tests for each group for each variable WITHIN EACH COMMUNITY
# significance tests and p-value adjustments only apply to non-NA groups

# alpha_diversity_summary.df <- NULL
# alpha_diversity_summary_genus.df <- NULL
# 
# alpha_diversity_significances.df <- NULL
# alpha_diversity_significances_genus.df <- NULL
# 
# for (myvar in discrete_variables){
#   if (myvar == "Remote_Community") {next}
#   for (community in unique(otu_rare_alpha.df$Remote_Community)){
#     data_subset <- subset(otu_rare_alpha.df, Remote_Community == community)
#     data_subset_decontaminated <- subset(otu_decontaminated_rare_alpha.df, Remote_Community == community)
# 
#     data_subset_genus <- subset(genus_rare_alpha.df, Remote_Community == community)
#     data_subset_genus_decontaminated <- subset(genus_decontaminated_rare_alpha.df, Remote_Community == community)
#     
#     if (is.null(alpha_diversity_summary.df)){
#       # Normal, OTU level
#       alpha_diversity_summary.df <- summarise_alpha_diversities(data_subset, myvar)
#       alpha_diversity_summary.df <- melt(alpha_diversity_summary.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
#       alpha_diversity_summary.df$Remote_Community <- community
#       
#       alpha_diversity_significances.df <- calculate_alpha_diversity_significance(data_subset, myvar)
#       alpha_diversity_significances.df$Remote_Community <- community
#       
#       # Normal, Genus level
#       alpha_diversity_summary_genus.df <- summarise_alpha_diversities(data_subset_genus, myvar)
#       alpha_diversity_summary_genus.df <- melt(alpha_diversity_summary_genus.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
#       alpha_diversity_summary_genus.df$Remote_Community <- community
#       
#       alpha_diversity_significances_genus.df <- calculate_alpha_diversity_significance(data_subset_genus, myvar)
#       alpha_diversity_significances_genus.df$Remote_Community <- community
#       
#       
#     } else{
#       # Normal, OTU level
#       temp_summary <- melt(summarise_alpha_diversities(data_subset, myvar), measure.vars = myvar, value.name = "Group", variable.name = "Variable")
#       temp_summary$Remote_Community <- community
#       alpha_diversity_summary.df <- rbind(alpha_diversity_summary.df, temp_summary)
#       
#       temp_significance <- calculate_alpha_diversity_significance(data_subset, myvar)
#       temp_significance$Remote_Community <- community
#       alpha_diversity_significances.df <- rbind(alpha_diversity_significances.df, temp_significance)
#       
#       # Normal, Genus level
#       temp_summary <- melt(summarise_alpha_diversities(data_subset_genus, myvar),measure.vars = myvar, value.name = "Group", variable.name = "Variable")
#       temp_summary$Remote_Community <- community
#       alpha_diversity_summary_genus.df <- rbind(alpha_diversity_summary_genus.df, temp_summary)
#       
#       temp_significance <- calculate_alpha_diversity_significance(data_subset_genus, myvar)
#       temp_significance$Remote_Community <- community
#       alpha_diversity_significances_genus.df <- rbind(alpha_diversity_significances_genus.df, temp_significance)
#       
#       # Decontaminated, OTU level
#       temp_summary <- melt(summarise_alpha_diversities(data_subset_decontaminated, myvar),measure.vars = myvar, value.name = "Group", variable.name = "Variable")
#       temp_summary$Remote_Community <- community
#       alpha_diversity_summary_decontaminated.df <- rbind(alpha_diversity_summary_decontaminated.df, temp_summary)
#       
#       temp_significance <- calculate_alpha_diversity_significance(data_subset_decontaminated, myvar)
#       temp_significance$Remote_Community <- community
#       alpha_diversity_significances_decontaminated.df <- rbind(alpha_diversity_significances_decontaminated.df, temp_significance)
#       
#       # Decontaminated, Genus level
#       temp_summary <- melt(summarise_alpha_diversities(data_subset_genus_decontaminated, myvar),measure.vars = myvar, value.name = "Group", variable.name = "Variable")
#       temp_summary$Remote_Community <- community
#       alpha_diversity_summary_genus_decontaminated.df <- rbind(alpha_diversity_summary_genus_decontaminated.df, temp_summary)
#       
#       temp_significance <- calculate_alpha_diversity_significance(data_subset_genus_decontaminated, myvar)
#       temp_significance$Remote_Community <- community
#       alpha_diversity_significances_genus_decontaminated.df <- rbind(alpha_diversity_significances_genus_decontaminated.df, temp_significance)
#       
#     }
#   }
# }
# # Normal, OTU level
# write.csv(x = alpha_diversity_summary.df, 
#           file = paste0("Result_tables/diversity_analysis/otu/otu_alpha_diversities_summary_within_community.csv"), 
#           quote = F, row.names = F)
# 
# write.csv(x = alpha_diversity_significances.df, 
#           file = paste0("Result_tables/diversity_analysis/otu/otu_alpha_diversities_significance_within_community.csv"), 
#           quote = F, row.names = F)
# 
# # Normal, Genus level
# write.csv(x = alpha_diversity_summary_genus.df, 
#           file = paste0("Result_tables/diversity_analysis/genus/genus_alpha_diversities_summary_within_community.csv"), 
#           quote = F, row.names = F)
# 
# write.csv(x = alpha_diversity_significances_genus.df, 
#           file = paste0("Result_tables/diversity_analysis/genus/genus_alpha_diversities_significance_within_community.csv"), 
#           quote = F, row.names = F)
# 
# # Decontaminated, OTU level
# write.csv(x = alpha_diversity_summary_decontaminated.df, 
#           file = paste0("Result_tables/diversity_analysis/otu_decontaminated/otu_alpha_diversities_summary_within_community_decontaminated.csv"), 
#           quote = F, row.names = F)
# 
# write.csv(x = alpha_diversity_significances_decontaminated.df, 
#           file = paste0("Result_tables/diversity_analysis/otu_decontaminated/otu_alpha_diversities_significance_within_community_decontaminated.csv"), 
#           quote = F, row.names = F)
# 
# # Decontaminated, Genus level
# write.csv(x = alpha_diversity_summary_decontaminated.df, 
#           file = paste0("Result_tables/diversity_analysis/genus_decontaminated/genus_alpha_diversities_summary_within_community_decontaminated.csv"), 
#           quote = F, row.names = F)
# 
# write.csv(x = alpha_diversity_significances_genus_decontaminated.df, 
#           file = paste0("Result_tables/diversity_analysis/genus_decontaminated/genus_alpha_diversities_significance_within_community_decontaminated.csv"), 
#           quote = F, row.names = F)
# 

# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
# Make figures
# generate_diversity_boxplot <- function(mydata, variable, metric, variable_colours_available = T, fill_palette = NULL){
#   internal_data.df <- mydata[!is.na(mydata[variable]),]
#   internal_data.df[[variable]] <- factor(internal_data.df[[variable]])
#   variable_values <- factor(as.character(unique(internal_data.df[[variable]])))
#   if (variable_colours_available == T){
#     color_col_name <- paste0(variable, "_colour")
#     variable_colours <- setNames(as.character(unique(internal_data.df[[color_col_name]])), as.character(unique(internal_data.df[[variable]])))
#   } else{
#     if (is.null(fill_palette)){
#       internal_colour_palette <- my_colour_palette_206_distinct
#     } else{
#       internal_colour_palette <- fill_palette
#     }
#     variable_colours <- setNames(internal_colour_palette[1:length(variable_values)], variable_values)  
#   }
#   myplot <- ggplot(internal_data.df, aes(x = get(variable), y = get(metric))) +
#     geom_boxplot(outlier.shape = NA, aes(fill = get(variable))) +
#     scale_fill_manual(values = variable_colours, name = variable) +
#     geom_jitter(size=0.5, width = 0.10, height=0) +
#     guides(fill=FALSE) +
#     # scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5,.5)) +
#     xlab(variable) +
#     ylab(metric)  +
#     common_theme +
#     theme(panel.border = element_blank(), 
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black", size = 0.5),
#           panel.background = element_blank(),
#           strip.background = element_rect(fill = "white", colour = "white", size = 1),
#           strip.text = element_text(size = 6),
#           legend.key=element_blank(),
#           legend.direction="vertical",
#           legend.background = element_rect(colour ="white", size = .3),
#           legend.text.align = 0,
#           legend.title = element_text(size=10, face="bold"),
#           legend.title.align = 0.5,
#           legend.margin = margin(c(2,2,2,2)),
#           legend.key.height=unit(.4,"cm"),
#           legend.text = element_text(size = 8),
#           axis.text = element_text(size = 9, colour = "black"),
#           axis.text.x = element_text(angle = 90, vjust = .5),
#           axis.title = element_text(size = 10,face = "bold"),
#           complete = F,
#           plot.title = element_text(size = 6,hjust = 0.5))
#   myplot
# }
# 
# # generate_diversity_boxplot(otu_rare_alpha.df, variable = "Remote_Community_OM_Classification",fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
# #   guides(fill = F, color = F) + 
# #   ggtitle("Chao1") +
# #   scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
# 
# 
# for (myvar in discrete_variables){
#   
#   # Normal, otu level
#   myplot <- generate_diversity_boxplot(otu_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Chao1") +
#     scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#   
#   ggsave(filename = paste0("Result_figures/diversity_analysis/otu/otu_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(otu_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Shannon") +
#     scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/otu/otu_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(otu_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#     guides(fill = F, color = F) +
#     ggtitle("Simpson") +
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/otu/otu_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   # Normal, genus level
#   myplot <- generate_diversity_boxplot(genus_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Chao1") +
#     scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#   
#   ggsave(filename = paste0("Result_figures/diversity_analysis/genus/genus_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(genus_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Shannon") +
#     scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/genus/genus_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(genus_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#     guides(fill = F, color = F) +
#     ggtitle("Simpson") +
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/genus/genus_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   # Decontaminated, otu level
#   myplot <- generate_diversity_boxplot(otu_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Chao1") +
#     scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#   
#   ggsave(filename = paste0("Result_figures/diversity_analysis/otu_decontaminated/otu_decontaminated_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(otu_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Shannon") +
#     scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/otu_decontaminated/otu_decontaminated_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(otu_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#     guides(fill = F, color = F) +
#     ggtitle("Simpson") +
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/otu_decontaminated/otu_decontaminated_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   # Decontaminated, genus level
#   myplot <- generate_diversity_boxplot(genus_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Chao1") +
#     scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#   
#   ggsave(filename = paste0("Result_figures/diversity_analysis/genus_decontaminated/genus_decontaminated_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(genus_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#     guides(fill = F, color = F) + 
#     ggtitle("Shannon") +
#     scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/genus_decontaminated/genus_decontaminated_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#   
#   myplot <- generate_diversity_boxplot(genus_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#     guides(fill = F, color = F) +
#     ggtitle("Simpson") +
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#   ggsave(filename = paste0("Result_figures/diversity_analysis/genus_decontaminated/genus_decontaminated_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
# }
# 
# 
# for (myvar in discrete_variables){
#   if (myvar == "Remote_Community") {next}
#   for (community in unique(otu_rare_alpha.df$Remote_Community)){
#     data_subset <- subset(otu_rare_alpha.df, Remote_Community == community)
#     data_subset_decontaminated <- subset(otu_decontaminated_rare_alpha.df, Remote_Community == community)
#     
#     data_subset_genus <- subset(genus_rare_alpha.df, Remote_Community == community)
#     data_subset_genus_decontaminated <- subset(genus_decontaminated_rare_alpha.df, Remote_Community == community)
#     # Normal, otu level
#     myplot <- generate_diversity_boxplot(data_subset, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Chao1, Community ", community)) +
#       scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#     
#     ggsave(filename = paste0("Result_figures/diversity_analysis/otu_within_community/otu_within_community_",community,"_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Shannon, Community ", community)) +
#       scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/otu_within_community/otu_within_community_",community,"_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#       guides(fill = F, color = F) +
#       ggtitle(paste0("Simpson, Community ", community)) +
#       scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/otu_within_community/otu_within_community_",community,"_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     # Normal, genus level
#     myplot <- generate_diversity_boxplot(data_subset_genus, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Chao1, Community ", community)) +
#       scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#     
#     ggsave(filename = paste0("Result_figures/diversity_analysis/genus_within_community/genus_within_community_",community,"_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset_genus, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Shannon, Community ", community)) +
#       scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/genus_within_community/genus_within_community_",community,"_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset_genus, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#       guides(fill = F, color = F) +
#       ggtitle(paste0("Simpson, Community ", community)) +
#       scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/genus_within_community/genus_within_community_",community,"_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     # Decontaminated, otu level
#     myplot <- generate_diversity_boxplot(data_subset_decontaminated, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Chao1, Community ", community)) +
#       scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#     
#     ggsave(filename = paste0("Result_figures/diversity_analysis/otu_within_community_decontaminated/otu_within_community_",community,"_decontaminated_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset_decontaminated, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Shannon, Community ", community)) +
#       scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/otu_within_community_decontaminated/otu_within_community_",community,"_decontaminated_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset_decontaminated, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#       guides(fill = F, color = F) +
#       ggtitle(paste0("Simpson, Community ", community)) +
#       scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/otu_within_community_decontaminated/otu_within_community_",community,"_decontaminated_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     # Decontaminated, genus level
#     myplot <- generate_diversity_boxplot(data_subset_genus_decontaminated, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Chao1, Community ", community)) +
#       scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
#     
#     ggsave(filename = paste0("Result_figures/diversity_analysis/genus_within_community_decontaminated/genus_within_community_",community,"_decontaminated_",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset_genus_decontaminated, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
#       guides(fill = F, color = F) + 
#       ggtitle(paste0("Shannon, Community ", community)) +
#       scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/genus_within_community_decontaminated/genus_within_community_",community,"_decontaminated_",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
#     
#     myplot <- generate_diversity_boxplot(data_subset_genus_decontaminated, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
#       guides(fill = F, color = F) +
#       ggtitle(paste0("Simpson, Community ", community)) +
#       scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
#     ggsave(filename = paste0("Result_figures/diversity_analysis/genus_within_community_decontaminated/genus_within_community_",community,"_decontaminated_",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
#   }
# }


# ---------------------------------------------
# Calculate the beta-diversity for each variable 
# Centre-log transform the counts first and use a euclidean distance. This should be equivalent or superior to 
# a bray curtis transform/distance used on counts. 
# As far as I can tell in the literature, e.g.
# Calle M Luz. (2019). Statistical Analysis of Metagenomics Data. Genomics Inform, 17(1), e6–.
# the euclidean distance between CLR values is an appropriate beta diversity measure




