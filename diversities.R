# Diversity calculations (Shannon, Chao1, Simpson) for each sample.
# Significance tests of each discrete group comparing diversity indices.
# Generate boxplots for discrete 
# Comparison of continuous variables vs diversity indices

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

# library(FSA)
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

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project/")
source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index
rownames(metadata_decontaminated.df) <- metadata_decontaminated.df$Index

# Load the counts
otu.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts.csv", header =T, row.names = 1))
genus.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts.csv", header =T, row.names = 1))

otu_decontaminated.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts_decontaminated.csv", header =T, row.names = 1))
genus_decontaminated.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts_decontaminated.csv", header =T, row.names = 1))

# Order the matrices and metadata to be the same order
otu.m <- otu.m[,rownames(metadata.df)]
genus.m <- genus.m[,rownames(metadata.df)]

otu_decontaminated.m <- otu_decontaminated.m[,rownames(metadata_decontaminated.df)]
genus_decontaminated.m <- genus_decontaminated.m[,rownames(metadata_decontaminated.df)]

# Create the rarefied matrix
otu_rare.m <- t(rrarefy(t(otu.m[,colSums(otu.m) >= 2000]), 2000))
genus_rare.m <- t(rrarefy(t(genus.m[,colSums(genus.m) >= 2000]), 2000))

otu_decontaminated_rare.m <- t(rrarefy(t(otu_decontaminated.m[,colSums(otu_decontaminated.m) >= 2000]), 2000))
genus_decontaminated_rare.m <- t(rrarefy(t(genus_decontaminated.m[,colSums(genus_decontaminated.m) >= 2000]), 2000))

# Define the discrete variables
# discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season","Nose", 
                        # "Otitis_status_OM_6mo","Remote_Community_Otitis_status","OM_6mo_Type_OM","Remote_Community_Season")
discrete_variables <- c("Remote_Community","Gold_Star","OM_6mo","Season","Nose","OM_Classification", "Remote_Community_Season",
                        "Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae",
                        "Remote_Community_OM_Classification")

# create phyloseq object
otu_rare_phyloseq <- otu_table(otu_rare.m, taxa_are_rows=TRUE)
genus_rare_phyloseq <- otu_table(genus_rare.m, taxa_are_rows=TRUE)

otu_decontaminated_rare_phyloseq <- otu_table(otu_decontaminated_rare.m, taxa_are_rows=TRUE)
genus_decontaminated_rare_phyloseq <- otu_table(genus_decontaminated_rare.m, taxa_are_rows=TRUE)

# Estimate alpha diversities
otu_rare_alpha.df <- estimate_richness(otu_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_rare_alpha.df <- otu_rare_alpha.df[rownames(metadata.df),]

otu_decontaminated_rare_alpha.df <- estimate_richness(otu_decontaminated_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_decontaminated_rare_alpha.df <- otu_decontaminated_rare_alpha.df[rownames(metadata_decontaminated.df),]

genus_rare_alpha.df <- estimate_richness(genus_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
genus_rare_alpha.df <- genus_rare_alpha.df[rownames(metadata.df),]

genus_decontaminated_rare_alpha.df <- estimate_richness(genus_decontaminated_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
genus_decontaminated_rare_alpha.df <- genus_decontaminated_rare_alpha.df[rownames(metadata_decontaminated.df),]


# Combine with metadata

otu_rare_alpha.df <- left_join(metadata.df[c("Index", 
                                             discrete_variables, 
                                             grep("_colour", names(metadata.df), value = T))],m2df(otu_rare_alpha.df, "Index"), by = "Index")
genus_rare_alpha.df <- left_join(metadata.df[c("Index",
                                                   discrete_variables, 
                                                   grep("_colour", names(metadata.df), value = T))],m2df(genus_rare_alpha.df, "Index"), by = "Index")

otu_decontaminated_rare_alpha.df <- left_join(metadata_decontaminated.df[c("Index", 
                                                                           discrete_variables, 
                                                                           grep("_colour", names(metadata.df), value = T))],m2df(otu_decontaminated_rare_alpha.df, "Index"), by = "Index")
genus_decontaminated_rare_alpha.df <- left_join(metadata_decontaminated.df[c("Index", 
                                                                                 discrete_variables, 
                                                                                 grep("_colour", names(metadata.df), value = T))],m2df(genus_decontaminated_rare_alpha.df, "Index"), by = "Index")


# Write per-sample diversities to file
write.csv(otu_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_otu_alpha_diversities.csv", quote = F, row.names = F
)
write.csv(genus_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_genus_alpha_diversities.csv", quote = F, row.names = F
)

write.csv(otu_decontaminated_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_otu_alpha_diversities_decontaminated.csv", quote = F, row.names = F
)
write.csv(genus_decontaminated_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_genus_alpha_diversities_decontaminated.csv", quote = F, row.names = F
)

# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
# temp <- summarise_alpha_diversities(otu_rare_alpha.df, "Gold_Star")
# melt(temp, measure.vars = "Gold_Star", value.name = "Group", variable.name = "Variable")
# otu_rare_alpha.df %>% dplyr::group_by(Gold_Star) %>% dplyr::summarise(Shannon_Mean =mean(Shannon))

# Calculate summaries and signifiance tests for each group for each variable
# significance tests and p-value adjustments only apply to non-NA groups
alpha_diversity_summary.df <- NULL
alpha_diversity_summary_decontaminated.df <- NULL

alpha_diversity_summary_genus.df <- NULL
alpha_diversity_summary_genus_decontaminated.df <- NULL


alpha_diversity_significances.df <- NULL
alpha_diversity_significances_decontaminated.df <- NULL

alpha_diversity_significances_genus.df <- NULL
alpha_diversity_significances_genus_decontaminated.df <- NULL

for (myvar in discrete_variables){
  if (is.null(alpha_diversity_summary.df)){
    # Normal, OTU level
    alpha_diversity_summary.df <- summarise_alpha_diversities(otu_rare_alpha.df, myvar)
    alpha_diversity_summary.df <- melt(alpha_diversity_summary.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")

    alpha_diversity_significances.df <- calculate_alpha_diversity_significance(otu_rare_alpha.df, myvar)
    
    # Normal, Genus level
    alpha_diversity_summary_genus.df <- summarise_alpha_diversities(genus_rare_alpha.df, myvar)
    alpha_diversity_summary_genus.df <- melt(alpha_diversity_summary_genus.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
    
    alpha_diversity_significances_genus.df <- calculate_alpha_diversity_significance(genus_rare_alpha.df, myvar)
    
    # Decontaminated, OTU level
    alpha_diversity_summary_decontaminated.df <- summarise_alpha_diversities(otu_decontaminated_rare_alpha.df, myvar)
    alpha_diversity_summary_decontaminated.df <- melt(alpha_diversity_summary_decontaminated.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
    
    alpha_diversity_significances_decontaminated.df <- calculate_alpha_diversity_significance(otu_decontaminated_rare_alpha.df, myvar)
    
    # Decontaminated, Genus level
    alpha_diversity_summary_genus_decontaminated.df <- summarise_alpha_diversities(genus_decontaminated_rare_alpha.df, myvar)
    alpha_diversity_summary_genus_decontaminated.df <- melt(alpha_diversity_summary_genus_decontaminated.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
    
    alpha_diversity_significances_genus_decontaminated.df <- calculate_alpha_diversity_significance(genus_decontaminated_rare_alpha.df, myvar)
    
  } else{
    # Normal, OTU level
    alpha_diversity_summary.df <- rbind(alpha_diversity_summary.df,
                                        melt(summarise_alpha_diversities(otu_rare_alpha.df, myvar  ),
                                             measure.vars = myvar, value.name = "Group", variable.name = "Variable"))
    alpha_diversity_significances.df <- rbind(alpha_diversity_significances.df, 
                                              calculate_alpha_diversity_significance(otu_rare_alpha.df, myvar))
    # Normal, Genus level
    alpha_diversity_summary_genus.df <- rbind(alpha_diversity_summary_genus.df,
                                              melt(summarise_alpha_diversities(genus_rare_alpha.df, myvar  ),
                                                   measure.vars = myvar, value.name = "Group", variable.name = "Variable"))
    alpha_diversity_significances_genus.df <- rbind(alpha_diversity_significances_genus.df, 
                                                    calculate_alpha_diversity_significance(genus_rare_alpha.df, myvar))
    # Decontaminated, OTU level
    alpha_diversity_summary_decontaminated.df <- rbind(alpha_diversity_summary_decontaminated.df,
                                                       melt(summarise_alpha_diversities(otu_decontaminated_rare_alpha.df, myvar),
                                                            measure.vars = myvar, value.name = "Group", variable.name = "Variable"))
    alpha_diversity_significances_decontaminated.df <- rbind(alpha_diversity_significances_decontaminated.df, 
                                                             calculate_alpha_diversity_significance(otu_decontaminated_rare_alpha.df, myvar))
    
    # Decontaminated, Genus level
    alpha_diversity_summary_genus_decontaminated.df <- rbind(alpha_diversity_summary_genus_decontaminated.df,
                                                             melt(summarise_alpha_diversities(genus_decontaminated_rare_alpha.df, myvar),
                                                                  measure.vars = myvar, value.name = "Group", variable.name = "Variable"))
    alpha_diversity_significances_genus_decontaminated.df <- rbind(alpha_diversity_significances_genus_decontaminated.df, 
                                                                   calculate_alpha_diversity_significance(genus_decontaminated_rare_alpha.df, myvar))
  }
}
# Normal, OTU level
write.csv(x = alpha_diversity_summary.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_significance.csv"), 
          quote = F, row.names = F)

# Normal, Genus level
write.csv(x = alpha_diversity_summary_genus.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_summary.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_significance.csv"), 
          quote = F, row.names = F)

# Decontaminated, OTU level
write.csv(x = alpha_diversity_summary_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary_decontaminated.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_significance_decontaminated.csv"), 
          quote = F, row.names = F)

# Decontaminated, Genus level
write.csv(x = alpha_diversity_summary_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_summary_decontaminated.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_significance_decontaminated.csv"), 
          quote = F, row.names = F)


# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
# Calculate summaries and signifiance tests for each group for each variable WITHIN EACH COMMUNITY
# significance tests and p-value adjustments only apply to non-NA groups

alpha_diversity_summary.df <- NULL
alpha_diversity_summary_decontaminated.df <- NULL

alpha_diversity_summary_genus.df <- NULL
alpha_diversity_summary_genus_decontaminated.df <- NULL

alpha_diversity_significances.df <- NULL
alpha_diversity_significances_decontaminated.df <- NULL

alpha_diversity_significances_genus.df <- NULL
alpha_diversity_significances_genus_decontaminated.df <- NULL

for (myvar in discrete_variables){
  if (myvar == "Remote_Community") {next}
  for (community in unique(otu_rare_alpha.df$Remote_Community)){
    data_subset <- subset(otu_rare_alpha.df, Remote_Community == community)
    data_subset_decontaminated <- subset(otu_decontaminated_rare_alpha.df, Remote_Community == community)

    data_subset_genus <- subset(genus_rare_alpha.df, Remote_Community == community)
    data_subset_genus_decontaminated <- subset(genus_decontaminated_rare_alpha.df, Remote_Community == community)
    
    if (is.null(alpha_diversity_summary.df)){
      # Normal, OTU level
      alpha_diversity_summary.df <- summarise_alpha_diversities(data_subset, myvar)
      alpha_diversity_summary.df <- melt(alpha_diversity_summary.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      alpha_diversity_summary.df$Remote_Community <- community
      
      alpha_diversity_significances.df <- calculate_alpha_diversity_significance(data_subset, myvar)
      alpha_diversity_significances.df$Remote_Community <- community
      
      # Normal, Genus level
      alpha_diversity_summary_genus.df <- summarise_alpha_diversities(data_subset_genus, myvar)
      alpha_diversity_summary_genus.df <- melt(alpha_diversity_summary_genus.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      alpha_diversity_summary_genus.df$Remote_Community <- community
      
      alpha_diversity_significances_genus.df <- calculate_alpha_diversity_significance(data_subset_genus, myvar)
      alpha_diversity_significances_genus.df$Remote_Community <- community
      
      # Decontaminated, OTU level
      alpha_diversity_summary_decontaminated.df <- summarise_alpha_diversities(data_subset_decontaminated, myvar)
      alpha_diversity_summary_decontaminated.df <- melt(alpha_diversity_summary_decontaminated.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      alpha_diversity_summary_decontaminated.df$Remote_Community <- community
      
      alpha_diversity_significances_decontaminated.df <- calculate_alpha_diversity_significance(data_subset_decontaminated, myvar)
      alpha_diversity_significances_decontaminated.df$Remote_Community <- community
      
      # Decontaminated, Genus level
      alpha_diversity_summary_genus_decontaminated.df <- summarise_alpha_diversities(data_subset_genus_decontaminated, myvar)
      alpha_diversity_summary_genus_decontaminated.df <- melt(alpha_diversity_summary_genus_decontaminated.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      alpha_diversity_summary_genus_decontaminated.df$Remote_Community <- community
      
      alpha_diversity_significances_genus_decontaminated.df <- calculate_alpha_diversity_significance(data_subset_genus_decontaminated, myvar)
      alpha_diversity_significances_genus_decontaminated.df$Remote_Community <- community
      
    } else{
      # Normal, OTU level
      temp_summary <- melt(summarise_alpha_diversities(data_subset, myvar), measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      temp_summary$Remote_Community <- community
      alpha_diversity_summary.df <- rbind(alpha_diversity_summary.df, temp_summary)
      
      temp_significance <- calculate_alpha_diversity_significance(data_subset, myvar)
      temp_significance$Remote_Community <- community
      alpha_diversity_significances.df <- rbind(alpha_diversity_significances.df, temp_significance)
      
      # Normal, Genus level
      temp_summary <- melt(summarise_alpha_diversities(data_subset_genus, myvar),measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      temp_summary$Remote_Community <- community
      alpha_diversity_summary_genus.df <- rbind(alpha_diversity_summary_genus.df, temp_summary)
      
      temp_significance <- calculate_alpha_diversity_significance(data_subset_genus, myvar)
      temp_significance$Remote_Community <- community
      alpha_diversity_significances_genus.df <- rbind(alpha_diversity_significances_genus.df, temp_significance)
      
      # Decontaminated, OTU level
      temp_summary <- melt(summarise_alpha_diversities(data_subset_decontaminated, myvar),measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      temp_summary$Remote_Community <- community
      alpha_diversity_summary_decontaminated.df <- rbind(alpha_diversity_summary_decontaminated.df, temp_summary)
      
      temp_significance <- calculate_alpha_diversity_significance(data_subset_decontaminated, myvar)
      temp_significance$Remote_Community <- community
      alpha_diversity_significances_decontaminated.df <- rbind(alpha_diversity_significances_decontaminated.df, temp_significance)
      
      # Decontaminated, Genus level
      temp_summary <- melt(summarise_alpha_diversities(data_subset_genus_decontaminated, myvar),measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      temp_summary$Remote_Community <- community
      alpha_diversity_summary_genus_decontaminated.df <- rbind(alpha_diversity_summary_genus_decontaminated.df, temp_summary)
      
      temp_significance <- calculate_alpha_diversity_significance(data_subset_genus_decontaminated, myvar)
      temp_significance$Remote_Community <- community
      alpha_diversity_significances_genus_decontaminated.df <- rbind(alpha_diversity_significances_genus_decontaminated.df, temp_significance)
      
    }
  }
}
# Normal, OTU level
write.csv(x = alpha_diversity_summary.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary_within_community.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_significance_within_community.csv"), 
          quote = F, row.names = F)

# Normal, Genus level
write.csv(x = alpha_diversity_summary_genus.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_summary_within_community.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_significance_within_community.csv"), 
          quote = F, row.names = F)

# Decontaminated, OTU level
write.csv(x = alpha_diversity_summary_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary_within_community_decontaminated.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_significance_within_community_decontaminated.csv"), 
          quote = F, row.names = F)

# Decontaminated, Genus level
write.csv(x = alpha_diversity_summary_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_summary_within_community_decontaminated.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/genus_alpha_diversities_significance_within_community_decontaminated.csv"), 
          quote = F, row.names = F)


# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
# Make figures
generate_diversity_boxplot <- function(mydata, variable, metric, variable_colours_available = T, fill_palette = NULL){
  internal_data.df <- mydata[!is.na(mydata[variable]),]
  internal_data.df[[variable]] <- factor(internal_data.df[[variable]])
  variable_values <- factor(as.character(unique(internal_data.df[[variable]])))
  if (variable_colours_available == T){
    color_col_name <- paste0(variable, "_colour")
    variable_colours <- setNames(as.character(unique(internal_data.df[[color_col_name]])), as.character(unique(internal_data.df[[variable]])))
  } else{
    if (is.null(fill_palette)){
      internal_colour_palette <- my_colour_palette_206_distinct
    } else{
      internal_colour_palette <- fill_palette
    }
    variable_colours <- setNames(internal_colour_palette[1:length(variable_values)], variable_values)  
  }
  myplot <- ggplot(internal_data.df, aes(x = get(variable), y = get(metric))) +
    geom_boxplot(outlier.shape = NA, aes(fill = get(variable))) +
    scale_fill_manual(values = variable_colours, name = variable) +
    geom_jitter(size=0.5, width = 0.10, height=0) +
    guides(fill=FALSE) +
    # scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5,.5)) +
    xlab(variable) +
    ylab(metric)  +
    common_theme +
    theme(panel.border = element_blank(), 
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
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.title = element_text(size = 10,face = "bold"),
          complete = F,
          plot.title = element_text(size = 6,hjust = 0.5))
  myplot
}

# generate_diversity_boxplot(otu_rare_alpha.df, variable = "Remote_Community_OM_Classification",fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
#   guides(fill = F, color = F) + 
#   ggtitle("Chao1") +
#   scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))


for (myvar in discrete_variables){
  
  # Normal, otu level
  myplot <- generate_diversity_boxplot(otu_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
    guides(fill = F, color = F) + 
    ggtitle("Chao1") +
    scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
  
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
  
  myplot <- generate_diversity_boxplot(otu_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
    guides(fill = F, color = F) + 
    ggtitle("Shannon") +
    scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
  
  myplot <- generate_diversity_boxplot(otu_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
    guides(fill = F, color = F) +
    ggtitle("Simpson") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
  
  # Decontaminated, otu level
  myplot <- generate_diversity_boxplot(otu_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Chao1",variable_colours_available = T) + 
    guides(fill = F, color = F) + 
    ggtitle("Chao1") +
    scale_y_continuous(limits = c(0,200), breaks = seq(0,200,50))
  
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_decontaminated_Chao1.pdf"),myplot, width = 10, height = 8,units = "cm")
  
  myplot <- generate_diversity_boxplot(otu_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Shannon",variable_colours_available = T) + 
    guides(fill = F, color = F) + 
    ggtitle("Shannon") +
    scale_y_continuous(limits = c(0,5), breaks = seq(0,5,.5))
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_decontaminated_Shannon.pdf"),myplot, width = 10, height = 8,units = "cm")
  
  myplot <- generate_diversity_boxplot(otu_decontaminated_rare_alpha.df, variable = myvar,fill_palette = my_colour_palette_10_distinct,metric = "Simpson",variable_colours_available = T) +
    guides(fill = F, color = F) +
    ggtitle("Simpson") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))
  ggsave(filename = paste0("Result_figures/diversity_analysis/",myvar,"_decontaminated_Simpson.pdf"),myplot, width = 10, height = 8,units = "cm")
}
