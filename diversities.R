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
otu_genus.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts.csv", header =T, row.names = 1))

otu_decontaminated.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts_decontaminated.csv", header =T, row.names = 1))
otu_genus_decontaminated.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts_decontaminated.csv", header =T, row.names = 1))


# Order the matrices and metadata to be the same order
otu.m <- otu.m[,rownames(metadata.df)]
otu_genus.m <- otu_genus.m[,rownames(metadata.df)]

otu_decontaminated.m <- otu_decontaminated.m[,rownames(metadata_decontaminated.df)]
otu_genus_decontaminated.m <- otu_genus_decontaminated.m[,rownames(metadata_decontaminated.df)]

# Create the rarefied matrix
otu_rare.m <- t(rrarefy(t(otu.m[,colSums(otu.m) >= 2000]), 2000))
otu_genus_rare.m <- t(rrarefy(t(otu_genus.m[,colSums(otu_genus.m) >= 2000]), 2000))

otu_decontaminated_rare.m <- t(rrarefy(t(otu_decontaminated.m[,colSums(otu_decontaminated.m) >= 2000]), 2000))
otu_genus_decontaminated_rare.m <- t(rrarefy(t(otu_genus_decontaminated.m[,colSums(otu_genus_decontaminated.m) >= 2000]), 2000))

# Define the discrete variables
discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season","Nose", 
                        "Otitis_status_OM_6mo","Remote_Community_Otitis_status","OM_6mo_Type_OM","Remote_Community_Season")


# create phyloseq object
otu_rare_phyloseq <- otu_table(otu_rare.m, taxa_are_rows=TRUE)
otu_genus_rare_phyloseq <- otu_table(otu_genus_rare.m, taxa_are_rows=TRUE)

otu_decontaminated_rare_phyloseq <- otu_table(otu_decontaminated_rare.m, taxa_are_rows=TRUE)
otu_genus_decontaminated_rare_phyloseq <- otu_table(otu_genus_decontaminated_rare.m, taxa_are_rows=TRUE)

# Estimate alpha diversities
otu_rare_alpha.df <- estimate_richness(otu_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_rare_alpha.df <- otu_rare_alpha.df[rownames(metadata.df),]

otu_decontaminated_rare_alpha.df <- estimate_richness(otu_decontaminated_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_decontaminated_rare_alpha.df <- otu_decontaminated_rare_alpha.df[rownames(metadata_decontaminated.df),]

otu_genus_rare_alpha.df <- estimate_richness(otu_genus_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_genus_rare_alpha.df <- otu_genus_rare_alpha.df[rownames(metadata.df),]

otu_genus_decontaminated_rare_alpha.df <- estimate_richness(otu_genus_decontaminated_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_genus_decontaminated_rare_alpha.df <- otu_genus_decontaminated_rare_alpha.df[rownames(metadata_decontaminated.df),]


# Combine with metadata
otu_rare_alpha.df <- left_join(metadata.df[c("Index", discrete_variables)],m2df(otu_rare_alpha.df, "Index"), by = "Index")
otu_genus_rare_alpha.df <- left_join(metadata.df[c("Index", discrete_variables)],m2df(otu_genus_rare_alpha.df, "Index"), by = "Index")

otu_decontaminated_rare_alpha.df <- left_join(metadata_decontaminated.df[c("Index", discrete_variables)],m2df(otu_decontaminated_rare_alpha.df, "Index"), by = "Index")
otu_genus_decontaminated_rare_alpha.df <- left_join(metadata_decontaminated.df[c("Index", discrete_variables)],m2df(otu_genus_decontaminated_rare_alpha.df, "Index"), by = "Index")


# Write per-sample diversities to file
write.csv(otu_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_otu_alpha_diversities.csv", quote = F, row.names = F
)
write.csv(otu_genus_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_genus_alpha_diversities.csv", quote = F, row.names = F
)

write.csv(otu_decontaminated_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_otu_alpha_diversities_decontaminated.csv", quote = F, row.names = F
)
write.csv(otu_genus_decontaminated_rare_alpha.df,
          "Result_tables/diversity_analysis/sample_genus_alpha_diversities_decontaminated.csv", quote = F, row.names = F
)


# temp <- summarise_alpha_diversities(otu_rare_alpha.df, "Gold_Star")
# melt(temp, measure.vars = "Gold_Star", value.name = "Group", variable.name = "Variable")
# otu_rare_alpha.df %>% dplyr::group_by(Gold_Star) %>% dplyr::summarise(Shannon_Mean =mean(Shannon))

# Calculate summaries and signifiance tests for each group for each variable
# significance tests and p-value adjustments only apply to non-NA groups
alpha_diversity_summary.df <- NULL
alpha_diversity_summary_decontaminated.df <- NULL

alpha_diversity_significances.df <- NULL
alpha_diversity_significances_decontaminated.df <- NULL

for (myvar in discrete_variables){
  if (is.null(alpha_diversity_summary.df)){
    alpha_diversity_summary.df <- summarise_alpha_diversities(otu_rare_alpha.df, myvar)
    alpha_diversity_summary.df <- melt(alpha_diversity_summary.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
    
    alpha_diversity_summary_decontaminated.df <- summarise_alpha_diversities(otu_decontaminated_rare_alpha.df, myvar)
    alpha_diversity_summary_decontaminated.df <- melt(alpha_diversity_summary_decontaminated.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
    
    alpha_diversity_significances.df <- calculate_alpha_diversity_significance(otu_rare_alpha.df, myvar)
    alpha_diversity_significances_decontaminated.df <- calculate_alpha_diversity_significance(otu_decontaminated_rare_alpha.df, myvar)

  } else{
    alpha_diversity_summary.df <- rbind(alpha_diversity_summary.df,
                                        melt(summarise_alpha_diversities(otu_rare_alpha.df, myvar  ),
                                             measure.vars = myvar, value.name = "Group", variable.name = "Variable"))
    alpha_diversity_summary_decontaminated.df <- rbind(alpha_diversity_summary_decontaminated.df,
                                                      melt(summarise_alpha_diversities(otu_decontaminated_rare_alpha.df, myvar),
                                                           measure.vars = myvar, value.name = "Group", variable.name = "Variable"))
    alpha_diversity_significances.df <- rbind(alpha_diversity_significances.df, 
                                              calculate_alpha_diversity_significance(otu_rare_alpha.df, myvar))
    alpha_diversity_significances_decontaminated.df <- rbind(alpha_diversity_significances_decontaminated.df, 
                                                             calculate_alpha_diversity_significance(otu_decontaminated_rare_alpha.df, myvar))
  }
}

write.csv(x = alpha_diversity_summary.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_summary_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary_decontaminated.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_significance.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_significance_decontaminated.csv"), 
          quote = F, row.names = F)


# Calculate summaries and signifiance tests for each group for each variable WITHIN EACH COMMUNITY
# significance tests and p-value adjustments only apply to non-NA groups

alpha_diversity_summary.df <- NULL
alpha_diversity_summary_decontaminated.df <- NULL

alpha_diversity_significances.df <- NULL
alpha_diversity_significances_decontaminated.df <- NULL

for (myvar in discrete_variables){
  if (myvar == "Remote_Community") {next}
  for (community in unique(otu_rare_alpha.df$Remote_Community)){
    data_subset <- subset(otu_rare_alpha.df, Remote_Community == community)
    data_subset_decontaminated <- subset(otu_decontaminated_rare_alpha.df, Remote_Community == community)
    if (is.null(alpha_diversity_summary.df)){
      
      alpha_diversity_summary.df <- summarise_alpha_diversities(data_subset, myvar)
      alpha_diversity_summary.df <- melt(alpha_diversity_summary.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")

      alpha_diversity_summary_decontaminated.df <- summarise_alpha_diversities(data_subset_decontaminated, myvar)
      alpha_diversity_summary_decontaminated.df <- melt(alpha_diversity_summary_decontaminated.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      
      alpha_diversity_significances.df <- calculate_alpha_diversity_significance(data_subset, myvar)
      alpha_diversity_significances_decontaminated.df <- calculate_alpha_diversity_significance(data_subset_decontaminated, myvar)
      
      alpha_diversity_summary.df$Remote_Community <- community
      alpha_diversity_summary.df$Remote_Community <- community
      alpha_diversity_summary_decontaminated.df$Remote_Community <- community
      alpha_diversity_significances.df$Remote_Community <- community
      alpha_diversity_significances_decontaminated.df$Remote_Community <- community
      
    } else{
      temp_summary <- melt(summarise_alpha_diversities(data_subset, myvar),
                           measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      temp_summary$Remote_Community <- community
      alpha_diversity_summary.df <- rbind(alpha_diversity_summary.df, temp_summary)
      
      temp_summary <- melt(summarise_alpha_diversities(data_subset_decontaminated, myvar),
                           measure.vars = myvar, value.name = "Group", variable.name = "Variable")
      temp_summary$Remote_Community <- community
      alpha_diversity_summary_decontaminated.df <- rbind(alpha_diversity_summary_decontaminated.df, temp_summary)
      
      temp_significance <- calculate_alpha_diversity_significance(data_subset, myvar)
      temp_significance$Remote_Community <- community
      alpha_diversity_significances.df <- rbind(alpha_diversity_significances.df, temp_significance)
      
      temp_significance <- calculate_alpha_diversity_significance(data_subset_decontaminated, myvar)
      temp_significance$Remote_Community <- community
      alpha_diversity_significances_decontaminated.df <- rbind(alpha_diversity_significances_decontaminated.df, temp_significance)
      
    }
  }
}
write.csv(x = alpha_diversity_summary.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary_within_community.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_summary_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_summary_within_community_decontaminated.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_within_community_significance.csv"), 
          quote = F, row.names = F)

write.csv(x = alpha_diversity_significances_decontaminated.df, 
          file = paste0("Result_tables/diversity_analysis/otu_alpha_diversities_significance_within_community_decontaminated.csv"), 
          quote = F, row.names = F)



# for (myvar in discrete_variables){
#   # Normal
#   write.csv(summarise_alpha_diversities(otu_rare_alpha.df, myvar), 
#             file = paste0("Result_tables/diversity_analysis/", myvar, "_otu_alpha_diversities_summary.csv"), quote = F, row.names = F)
#   write.csv(calculate_alpha_diversity_significance(otu_rare_alpha.df, myvar), 
#             file = paste0("Result_tables/diversity_analysis/", myvar, "_otu_alpha_diversities_significance.csv"), quote = F, row.names = F)
#   
#   write.csv(summarise_alpha_diversities(otu_genus_rare_alpha.df, myvar), 
#             file = paste0("Result_tables/diversity_analysis/", myvar, "_genus_alpha_diversities_summary.csv"), quote = F, row.names = F)
#   write.csv(calculate_alpha_diversity_significance(otu_genus_rare_alpha.df, myvar), 
#             file = paste0("Result_tables/diversity_analysis/", myvar, "_genus_alpha_diversities_significance.csv"), quote = F, row.names = F)
#   
#   # Decontaminated data
#   write_decontaminated.csv(summarise_alpha_diversities(otu_decontaminated_rare_alpha.df, myvar), 
#                            file = paste0("Result_tables/diversity_analysis/", myvar, "_otu_alpha_diversities_summary_decontaminated.csv"), quote = F, row.names = F)
#   write_decontaminated.csv(calculate_alpha_diversity_significance(otu_decontaminated_rare_alpha.df, myvar), 
#                            file = paste0("Result_tables/diversity_analysis/", myvar, "_otu_alpha_diversities_significance_decontaminated.csv"), quote = F, row.names = F)
#   
#   write_decontaminated.csv(summarise_alpha_diversities(otu_genus_decontaminated_rare_alpha.df, myvar), 
#                            file = paste0("Result_tables/diversity_analysis/", myvar, "_genus_alpha_diversities_summary_decontaminated.csv"), quote = F, row.names = F)
#   write_decontaminated.csv(calculate_alpha_diversity_significance(otu_genus_decontaminated_rare_alpha.df, myvar), 
#                            file = paste0("Result_tables/diversity_analysis/", myvar, "_genus_alpha_diversities_significance_decontaminated.csv"), quote = F, row.names = F)
# }




generate_diversity_boxplot <- function(mydata, variable, metric, variable_colours_available = T){
  internal_data.df <- mydata[!is.na(mydata[variable]),]
  variable_values <- factor(as.character(unique(internal_data.df[[variable]])))
  if (variable_colours_available == T){
    color_col_name <- paste0(variable, "_colour")
    variable_colours <- setNames(as.character(unique(internal_data.df[[color_col_name]])), as.character(unique(internal_data.df[[variable]])))
  } else{
    variable_colours <- setNames(my_colour_palette_206_distinct[1:length(variable_values)], variable_values)  
  }
  myplot <- ggplot(internal_data.df, aes(x = get(variable), y = get(metric))) +
    geom_boxplot(outlier.shape = NA, aes(fill = get(variable))) +
    scale_fill_manual(values = variable_colours, name = variable) +
    geom_jitter(size=0.5, width = 0.10, height=0) +
    guides(fill=FALSE) +
    # scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5,.5)) +
    xlab("") +
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
          plot.title = element_text(size = 6))
  myplot
}


# Add the colours back
# temp <- summarise_alpha_diversities(otu_genus_rare_alpha.df, "Commodity")
# temp <- left_join(temp,metadata.df[c("Commodity", grep("_colour", names(metadata.df), value = T))])

# FEATURE (OTU) LEVEL
temp <- otu_rare_alpha.df
temp <- left_join(temp,metadata.df[c("Index", grep("_colour", names(metadata.df), value = T))], by = "Index")

ggsave(plot = generate_diversity_boxplot(temp, "Commodity", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Commodity_otu_chao1.pdf"),width = 10, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Commodity", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Commodity_otu_shannon.pdf"),width = 10, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Commodity", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Commodity_otu_simpson.pdf"),width = 10, height = 12, units = "cm")

ggsave(plot = generate_diversity_boxplot(temp, "Sample_type", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_type_otu_chao1.pdf"),width = 15, height = 15, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_type", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_type_otu_shannon.pdf"),width = 15, height = 15, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_type", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_type_otu_simpson.pdf"),width = 15, height = 15, units = "cm")

ggsave(plot = generate_diversity_boxplot(temp, "Sample_treatment", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_treatment_otu_chao1.pdf"),width = 8, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_treatment", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_treatment_otu_shannon.pdf"),width = 8, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_treatment", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_treatment_otu_simpson.pdf"),width = 8, height = 12, units = "cm")

ggsave(plot = generate_diversity_boxplot(temp, "study_accession", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Study_accession_otu_chao1.pdf"),width = 18, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "study_accession", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Study_accession_otu_shannon.pdf"),width = 18, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "study_accession", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Study_accession_otu_simpson.pdf"),width = 18, height = 12, units = "cm")



# GENUS LEVEL
temp <- otu_genus_rare_alpha.df
temp <- left_join(temp,metadata.df[c("Index", grep("_colour", names(metadata.df), value = T))], by = "Index")

ggsave(plot = generate_diversity_boxplot(temp, "Commodity", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Commodity_genus_chao1.pdf"),width = 10, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Commodity", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Commodity_genus_shannon.pdf"),width = 10, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Commodity", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Commodity_genus_simpson.pdf"),width = 10, height = 12, units = "cm")

ggsave(plot = generate_diversity_boxplot(temp, "Sample_type", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_type_genus_chao1.pdf"),width = 15, height = 15, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_type", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_type_genus_shannon.pdf"),width = 15, height = 15, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_type", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_type_genus_simpson.pdf"),width = 15, height = 15, units = "cm")

ggsave(plot = generate_diversity_boxplot(temp, "Sample_treatment", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_treatment_genus_chao1.pdf"),width = 8, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_treatment", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_treatment_genus_shannon.pdf"),width = 8, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "Sample_treatment", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Sample_treatment_genus_simpson.pdf"),width = 8, height = 12, units = "cm")

ggsave(plot = generate_diversity_boxplot(temp, "study_accession", "Chao1", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Study_accession_genus_chao1.pdf"),width = 18, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "study_accession", "Shannon", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Study_accession_genus_shannon.pdf"),width = 18, height = 12, units = "cm")
ggsave(plot = generate_diversity_boxplot(temp, "study_accession", "Simpson", variable_colours_available = T),
       filename = paste0("Result_figures/diversity_analysis/Study_accession_genus_simpson.pdf"),width = 18, height = 12, units = "cm")


# ------------------------------------------------------------------------------------------------------------------------------
# Takes awhile to calculate, uncomment to run

# Calculate the beta-diversity for each variable.
# Centre-log transform the counts first and use a euclidean distance. This should be equivalent or superior to 
# a bray curtis distance used on transformed counts.
# As far as I can tell in the literature, the euclidean distance between CLR values is an appropriate beta diversity measure
# and will be equivalent to PCA.

beta_diversity_significances <- data.frame("Variable" = character(),
                                           "P_value" = numeric(),
                                           "R_value" = numeric()
)
for (myvar in c("Commodity","Sample_type","Sample_treatment", "study_accession")){
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  otu_genus_rare_subset.m <- otu_genus_rare.m[,rownames(metadata_subset.df)]
  # print(myvar)
  temp <- with(metadata_subset.df, anosim(t(clr(otu_genus_rare_subset.m)),get(myvar), distance = "euclidean",permutations = 999))
  beta_diversity_significances <- rbind(beta_diversity_significances, data.frame("Variable" = myvar,
                                                                                 "P_value" = temp$signif,
                                                                                 "R_value" = temp$statistic))
}

beta_diversity_significances$padj <- round(p.adjust(beta_diversity_significances$P_value,method = "BH"),6)
write.csv(beta_diversity_significances, file = "Result_tables/diversity_analysis/beta_diversity_significance.csv", row.names = F, quote = F)


temp <- with(metadata.df, anosim(t(clr(otu_genus_rare_subset.m)),Commodity, distance = "euclidean",permutations = 999))
summary(temp)
plot(temp)

pairwise.wilcox.test(temp$dis.rank, temp$class.vec, p.adj = "BH")
