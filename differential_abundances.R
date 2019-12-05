library(DESeq2)
# library(BiocParallel)




############################################################
######################## Functions #########################

# For each rowname (OTU), get the corresponding taxonomy_species
# Assumes "OTU.ID" and "taxonomy_species" columns in the provided map dataframe
assign_taxonomy_to_otu <- function(otutable, taxon_map){
  taxonomies <- c()
  for (otuid in rownames(otutable)){
    taxonomies <- c(taxonomies, as.character(taxon_map[taxon_map$OTU.ID == otuid,]$taxonomy_species))
  }
  return(taxonomies)
}

# Filter and sort DESeq result tables
filter_and_sort_dds_results <- function(x, p_value_threshold = 0.05){
  filtered_table <- x
  filtered_table <- filtered_table[!is.na(filtered_table$padj),]
  filtered_table <- filtered_table[filtered_table$padj <= p_value_threshold,]
  filtered_table <- filtered_table[order(filtered_table$padj),]
  return(filtered_table)
}


filter_matrix_rows <- function(my_matrix, row_max){
  rows_before <- dim(my_matrix)[1]
  rows_after <- dim(my_matrix[apply(my_matrix,1,max) >= row_max,])[1]
  print(paste0("Rows before = ", rows_before))
  print(paste0("Rows after = ", rows_after))
  print(paste0("Lost % = ", round((rows_before-rows_after)/rows_before*100, 2), "%"))
  return(my_matrix[apply(my_matrix,1,max) >= row_max,])
}

############################################################

# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project/")
source("Code/helper_functions.R")

# Load count table at the OTU level. These are the counts for OTUs that were above our abundance thresholds
# otu_decontaminated.m <- as.matrix(read.table("Result_tables/count_tables/OTU_counts_rarefied.csv", sep =",", header =T, row.names = 1))
# genus_rare.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts_rarefied.csv", sep =",", header =T, row.names = 1))

otu.m <- as.matrix(read.table("Result_tables/count_tables/OTU_counts.csv", sep =",", header =T, row.names = 1))
genus.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts.csv", sep =",", header =T, row.names = 1))

otu_decontaminated.m <- as.matrix(read.table("Result_tables/count_tables/OTU_counts_decontaminated.csv", sep =",", header =T, row.names = 1))
genus_decontaminated.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts_decontaminated.csv", sep =",", header =T, row.names = 1))

# Filter out OTUs/species/genera that do not have at # reads in at least one sample
# head(melt(sort(colSums(otu_decontaminated.m))))
# otu_decontaminated.m <- filter_matrix_rows(otu_decontaminated.m,50)
# head(melt(sort(colSums(otu_decontaminated.m))))
# genus_rare.m <- filter_matrix_rows(genus_rare.m,0)

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv", sep =",", header = T)

# Define the discrete variables
# discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season","Nose", 
#                         "Otitis_status_OM_6mo","Remote_Community_Otitis_status","OM_6mo_Type_OM","Remote_Community_Season")
discrete_variables <- c("Remote_Community","Gold_Star","OM_6mo","Season","Nose","OM_Classification", "Remote_Community_Season",
                        "Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae")

# Only keep columns (samples) in the metadata
# otu_decontaminated.m <- otu_decontaminated.m[,colnames(otu_decontaminated.m) %in% as.character(metadata_decontaminated.df$Index)]
# genus_rare.m <- genus_rare.m[,colnames(genus_rare.m) %in% as.character(metadata_decontaminated.df$Index)]
otu.m <- otu.m[,colnames(otu.m) %in% as.character(metadata.df$Index)]
genus.m <- genus.m[,colnames(genus.m) %in% as.character(metadata.df$Index)]

otu_decontaminated.m <- otu_decontaminated.m[,colnames(otu_decontaminated.m) %in% as.character(metadata_decontaminated.df$Index)]
genus_decontaminated.m <- genus_decontaminated.m[,colnames(genus_decontaminated.m) %in% as.character(metadata_decontaminated.df$Index)]

# Order the metadata_decontaminated.df by the index value
metadata.df <- metadata.df[order(metadata.df$Index),]
metadata_decontaminated.df <- metadata_decontaminated.df[order(metadata_decontaminated.df$Index),]

# Since we likely removed samples from the count matrix
# in the main script, remove them from the metadata_decontaminated.df here
# samples_removed <- metadata_decontaminated.df$Index[!metadata_decontaminated.df$Index %in% colnames(otu_decontaminated.m)]
# metadata_decontaminated.df <- metadata_decontaminated.df[! metadata_decontaminated.df$Index %in% samples_removed,]

# Rownames should match the sample columns in the otu table
rownames(metadata.df) <- metadata.df$Index
rownames(metadata_decontaminated.df) <- metadata_decontaminated.df$Index

# Order the otu_tables the same order as the metadata
otu.m <- otu.m[,rownames(metadata.df)]
genus.m <- genus.m[,rownames(metadata.df)]

otu_decontaminated.m <- otu_decontaminated.m[,rownames(metadata_decontaminated.df)]
genus_decontaminated.m <- genus_decontaminated.m[,rownames(metadata_decontaminated.df)]

dim(otu_decontaminated.m)
dim(genus_decontaminated.m)
dim(metadata_decontaminated.df)


# ---------------------------------------------------------------------------------------------------------
# Perform differential abundance calculations at the OTU level and genus level, 
# comparing between the groups within variables of interest

# DESeq requires a count matrix ('countData'), a corresponding metadata_decontaminated.df ('colData') and a 'design' formula. The formula expresses
# how the counts for each OTU/genus depend on the variables defined in the 'colData'. See help(DESeqDataSetFromMatrix) for more information.
# The first column of the metadata_decontaminated.df ('colData') must match the ordering of the columns of the countData

# Ensure names of the otu / genus count matrices match the order of the metadata_decontaminated.df!
# Assumes number of samples in metadata_decontaminated.df and count data are the same
all(colnames(otu_decontaminated.m) == metadata_decontaminated.df$Index) # Should be 'True'
all(colnames(otu_decontaminated.m) == rownames(metadata_decontaminated.df)) # Should be 'True'



# Convert variables to factors
metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)
metadata_decontaminated.df[discrete_variables] <- lapply(metadata_decontaminated.df[discrete_variables], factor)


compare_groups_deseq <- function(mydata.m, mymetadata.df, myvariables, assign_taxonomy = T){
  # Compare groups for all variables
  combined_results_ordered.df <- data.frame()
  for (myvar in myvariables){
    # Get all non-NA entries in the metadata
    mymetadata_filtered.df <- mymetadata.df[!is.na(mymetadata.df[,myvar]),]
    
    # Ensure factored variable
    mymetadata_filtered.df[,myvar] <- factor(mymetadata_filtered.df[,myvar])
    
    # Extract corresponding entries from data
    mydata_filtered.m <- mydata.m[,rownames(mymetadata_filtered.df)]

    # Run DESeq
    dds <- DESeqDataSetFromMatrix(countData = mydata_filtered.m, colData = mymetadata_filtered.df, design = as.formula(paste0("~", myvar)))
    geoMeans <- apply(counts(dds), 1, gm_mean)
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    dds <- try(DESeq(dds, test = "Wald", fitType = "parametric", parallel = T))
    group_combinations <- combn(sort(unique(mymetadata_filtered.df[,myvar])),2)
    
    for (i in 1:ncol(group_combinations)){
      group_1 <- as.character(group_combinations[1,i])
      group_2 <- as.character(group_combinations[2,i])
      
      # group_1_meta <- subset(full, get(myvar) == group_1)
      # group_2_meta <- subset(full, get(myvar) == group_2)
      # n_group_1 <- dim(group_1_meta)[1]
      # n_group_2 <- dim(group_2_meta)[1]
      
      n_group_1 <- dim(subset(mymetadata_filtered.df, get(myvar) == group_1))[1]
      n_group_2 <- dim(subset(mymetadata_filtered.df, get(myvar) == group_2))[1]
      
      # Extract results for contrasted groups
      resMFSource <- results(dds, contrast = c(myvar,group_1,group_2), alpha=0.01, independentFiltering = F, cooksCutoff = F,parallel = T)
      resMFSource$Group_1 <- group_1
      resMFSource$Group_2 <- group_2
      resMFSource$Variable <- myvar
      resMFSource$N_Group_1 <- n_group_1
      resMFSource$N_Group_2 <- n_group_2
      
      # Assign the taxonomy to the results
      if (assign_taxonomy == T){
        resMFSource$Taxonomy <- assign_taxonomy_to_otu(resMFSource, otu_taxonomy_map.df)   
        # Convert to dataframe
        resMFSource <- m2df(resMFSource, "OTU")
      } else{
        # Convert to dataframe
        resMFSource <- m2df(resMFSource, "Taxonomy")
      }
      # print(resMFSource)
      resMFSource <- filter_and_sort_dds_results(resMFSource, 0.01)
      combined_results_ordered.df <- rbind(combined_results_ordered.df, resMFSource)
    }
  }
  combined_results_ordered.df
}
otu_group_comparison.df <- compare_groups_deseq(mydata.m = otu.m, mymetadata.df = metadata.df, myvariables = discrete_variables, assign_taxonomy = T)
genus_group_comparison.df <- compare_groups_deseq(mydata.m = genus.m, mymetadata.df = metadata.df, myvariables = discrete_variables, assign_taxonomy = F)
otu_decontaminated_group_comparison.df <- compare_groups_deseq(mydata.m = otu_decontaminated.m, mymetadata.df = metadata_decontaminated.df, myvariables = discrete_variables, assign_taxonomy = T)
genus_decontaminated_group_comparison.df <- compare_groups_deseq(mydata.m = genus_decontaminated.m, mymetadata.df = metadata_decontaminated.df, myvariables = discrete_variables, assign_taxonomy = F)

write.csv(x =otu_group_comparison.df,file ="Result_tables/DESeq_results/OTU_deseq.csv",quote = F, row.names =F)
write.csv(x =genus_group_comparison.df,file ="Result_tables/DESeq_results/Genus_deseq.csv",quote = F, row.names =F)
write.csv(x =otu_decontaminated_group_comparison.df,file ="Result_tables/DESeq_results/OTU_deseq_decontaminated.csv",quote = F, row.names =F)
write.csv(x =genus_decontaminated_group_comparison.df,file ="Result_tables/DESeq_results/Genus_deseq_decontaminated.csv",quote = F, row.names =F)


compare_groups_deseq_within_group <- function(mydata.m, mymetadata.df, myvariables, within_group_variable, assign_taxonomy = F){
  combined_results.df <- data.frame()
  for (myvar_value in unique(metadata.df[,within_group_variable])){
    temp <- compare_groups_deseq(mydata.m = mydata.m, 
                                 mymetadata.df = subset(mymetadata.df, get(within_group_variable) == myvar_value), 
                                 myvariables = myvariables, 
                                 assign_taxonomy = assign_taxonomy)
    temp[,within_group_variable] <- myvar_value
    combined_results.df <- rbind(combined_results.df, temp)
  }
  combined_results.df
}

reduced_variables <- discrete_variables[which(!discrete_variables == "Remote_Community")]
otu_group_comparison_within_community.df <- compare_groups_deseq_within_group(mydata.m = otu.m, 
                                                                                mymetadata.df = metadata.df, 
                                                                                myvariables = reduced_variables, 
                                                                                within_group_variable = "Remote_Community", 
                                                                                assign_taxonomy = T)

genus_group_comparison_within_community.df <- compare_groups_deseq_within_group(mydata.m = genus.m, 
                                                                                mymetadata.df = metadata.df, 
                                                                                myvariables = reduced_variables, 
                                                                                within_group_variable = "Remote_Community", 
                                                                                assign_taxonomy = F)

otu_decontaminated_group_comparison_within_community.df <- compare_groups_deseq_within_group(mydata.m = otu_decontaminated.m, 
                                                                                             mymetadata.df = metadata_decontaminated.df, 
                                                                                             myvariables = reduced_variables, 
                                                                                             within_group_variable = "Remote_Community", 
                                                                                             assign_taxonomy = T)

genus_decontaminated_group_comparison_within_community.df <- compare_groups_deseq_within_group(mydata.m = genus_decontaminated.m, 
                                                                                               mymetadata.df = metadata_decontaminated.df, 
                                                                                               myvariables = reduced_variables, 
                                                                                               within_group_variable = "Remote_Community", 
                                                                                               assign_taxonomy = F)
write.csv(x =otu_group_comparison_within_community.df,file ="Result_tables/DESeq_results/OTU_deseq_within_community.csv",quote = F, row.names =F)
write.csv(x =genus_group_comparison_within_community.df,file ="Result_tables/DESeq_results/Genus_deseq_within_community.csv",quote = F, row.names =F)
write.csv(x =otu_decontaminated_group_comparison_within_community.df,file ="Result_tables/DESeq_results/OTU_deseq_within_community_decontaminated.csv",quote = F, row.names =F)
write.csv(x =genus_decontaminated_group_comparison_within_community.df,file ="Result_tables/DESeq_results/Genus_deseq_within_community_decontaminated.csv",quote = F, row.names =F)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Heatmap of abundances for significant taxa / features

# Genus level, within community
genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv",header = T)
genus_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata_decontaminated.csv",header = T)

# significant_taxa.v <- subset(genus_decontaminated_group_comparison_within_community.df, Variable == "Nose")$Taxonomy
significant_taxa.v <- genus_group_comparison_within_community.df$Taxonomy
significant_taxa_decontaminated.v <- genus_decontaminated_group_comparison_within_community.df$Taxonomy

# Generate matrix for heatmap
heatmap.m <- genus_data.df[c("Sample", "taxonomy_genus","Relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_genus %in% significant_taxa.v,]
heatmap.m <- heatmap.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- metadata.df[colnames(heatmap.m),]

heatmap_decontaminated.m <- genus_data_decontaminated.df[c("Sample", "taxonomy_genus","Relative_abundance")]
heatmap_decontaminated.m <- heatmap_decontaminated.m[heatmap_decontaminated.m$taxonomy_genus %in% significant_taxa_decontaminated.v,]
heatmap_decontaminated.m <- heatmap_decontaminated.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap_decontaminated.m <- df2matrix(heatmap_decontaminated.m)
heatmap_metadata_decontaminated.df <- metadata_decontaminated.df[colnames(heatmap_decontaminated.m),]

source("code/helper_functions.R")
make_heatmap(myheatmap_matrix = heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_genus_relative_abundance_within_community_heatmap.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Nose"),
             column_title = "Sample",
             row_title = "Genus",
             plot_height = 4,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


make_heatmap(myheatmap_matrix = heatmap_decontaminated.m*100, 
             mymetadata = heatmap_metadata_decontaminated.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_genus_relative_abundance_within_community_heatmap_decontaminated.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Nose"),
             column_title = "Sample",
             row_title = "Genus",
             plot_height = 4,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


# Feature level, within community
otu_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata.csv",header = T)
otu_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata_decontaminated.csv",header = T)

significant_taxa.v <- otu_group_comparison_within_community.df$OTU
# significant_taxa.v <- subset(otu_group_comparison_within_community.df, Variable == "Gold_Star")$OTU
significant_taxa_decontaminated.v <- otu_decontaminated_group_comparison_within_community.df$OTU
# significant_taxa_decontaminated.v <- subset(otu_decontaminated_group_comparison_within_community.df, Variable == "Gold_Star")$OTU

# Generate matrix for heatmap
heatmap.m <- otu_data.df[c("Sample", "OTU.ID","Relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$OTU %in% significant_taxa.v,]
heatmap.m <- heatmap.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- metadata.df[colnames(heatmap.m),]

heatmap_decontaminated.m <- otu_data_decontaminated.df[c("Sample", "OTU.ID","Relative_abundance")]
heatmap_decontaminated.m <- heatmap_decontaminated.m[heatmap_decontaminated.m$OTU %in% significant_taxa_decontaminated.v,]
heatmap_decontaminated.m <- heatmap_decontaminated.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap_decontaminated.m <- df2matrix(heatmap_decontaminated.m)
heatmap_metadata_decontaminated.df <- metadata_decontaminated.df[colnames(heatmap_decontaminated.m),]

make_heatmap(myheatmap_matrix = heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_OTU_relative_abundance_within_community_heatmap.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Gold_Star"),
             column_title = "Sample",
             row_title = "Sequence variant",
             plot_height = 20,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


make_heatmap(myheatmap_matrix = heatmap_decontaminated.m*100, 
             mymetadata = heatmap_metadata_decontaminated.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_OTU_relative_abundance_within_community_heatmap_decontaminated.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Gold_Star"),
             column_title = "Sample",
             row_title = "Sequence variant",
             plot_height = 20,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


# ----
# Genus level

significant_taxa.v <- genus_group_comparison.df$Taxonomy
significant_taxa_decontaminated.v <- genus_decontaminated_group_comparison.df$Taxonomy
# significant_taxa_decontaminated.v <- subset(genus_decontaminated_group_comparison.df, Variable == "Remote_Community")$Taxonomy

# Generate matrix for heatmap
heatmap.m <- genus_data.df[c("Sample", "taxonomy_genus","Relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$taxonomy_genus %in% significant_taxa.v,]
heatmap.m <- heatmap.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- metadata.df[colnames(heatmap.m),]

heatmap_decontaminated.m <- genus_data_decontaminated.df[c("Sample", "taxonomy_genus","Relative_abundance")]
heatmap_decontaminated.m <- heatmap_decontaminated.m[heatmap_decontaminated.m$taxonomy_genus %in% significant_taxa_decontaminated.v,]
heatmap_decontaminated.m <- heatmap_decontaminated.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap_decontaminated.m <- df2matrix(heatmap_decontaminated.m)
heatmap_metadata_decontaminated.df <- metadata_decontaminated.df[colnames(heatmap_decontaminated.m),]

source("code/helper_functions.R")
make_heatmap(myheatmap_matrix = heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_genus_relative_abundance_heatmap.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Nose"),
             column_title = "Sample",
             row_title = "Genus",
             plot_height = 4,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


make_heatmap(myheatmap_matrix = heatmap_decontaminated.m*100, 
             mymetadata = heatmap_metadata_decontaminated.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_genus_relative_abundance_heatmap_decontaminated.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Nose"),
             column_title = "Sample",
             row_title = "Genus",
             plot_height = 4,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


# Feature level
significant_taxa.v <- otu_group_comparison.df$OTU
# significant_taxa.v <- subset(otu_group_comparison.df, Variable == "Gold_Star")$OTU
significant_taxa_decontaminated.v <- otu_decontaminated_group_comparison.df$OTU
# significant_taxa_decontaminated.v <- subset(otu_decontaminated_group_comparison.df, Variable == "Gold_Star")$OTU

# Generate matrix for heatmap
heatmap.m <- otu_data.df[c("Sample", "OTU.ID","Relative_abundance")]
heatmap.m <- heatmap.m[heatmap.m$OTU %in% significant_taxa.v,]
heatmap.m <- heatmap.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap.m <- df2matrix(heatmap.m)
heatmap_metadata.df <- metadata.df[colnames(heatmap.m),]

heatmap_decontaminated.m <- otu_data_decontaminated.df[c("Sample", "OTU.ID","Relative_abundance")]
heatmap_decontaminated.m <- heatmap_decontaminated.m[heatmap_decontaminated.m$OTU %in% significant_taxa_decontaminated.v,]
heatmap_decontaminated.m <- heatmap_decontaminated.m %>% spread(Sample, Relative_abundance,fill = 0)
heatmap_decontaminated.m <- df2matrix(heatmap_decontaminated.m)
heatmap_metadata_decontaminated.df <- metadata_decontaminated.df[colnames(heatmap_decontaminated.m),]

make_heatmap(myheatmap_matrix = heatmap.m*100, 
             mymetadata = heatmap_metadata.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_OTU_relative_abundance_heatmap.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Gold_Star"),
             column_title = "Sample",
             row_title = "Sequence variant",
             plot_height = 10,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


make_heatmap(myheatmap_matrix = heatmap_decontaminated.m*100, 
             mymetadata = heatmap_metadata_decontaminated.df,
             filename = paste0("Result_figures/heatmaps/Sample_DESeq_OTU_relative_abundance_heatmap_decontaminated.pdf"),
             variables = grep("_Season", discrete_variables, value =T,invert = T),
             # variables = c("Remote_Community", "Gold_Star"),
             column_title = "Sample",
             row_title = "Sequence variant",
             plot_height = 10,
             plot_width = 20,
             cluster_columns = F,
             cluster_rows = T,
             column_title_size = 10,
             row_title_size = 10,
             annotation_name_size = 6,
             my_annotation_palette = my_colour_palette_15,
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, "> 60"),
             my_breaks = c(0, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,
             discrete_legend = T,
             legend_title = "Mean relative abundance %",
             palette_choice = 'purple',
             row_dend_width = unit(3, "cm"),
             simple_anno_size = unit(.25, "cm"),
             show_cell_values = F,
             cell_fun_value_col_threshold = 15
)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# WORK IN PROGRESS - heatmap of comparisons vs taxa

# x axis (col) = group_vs_group
# y axis (row) = taxa
# fill = log change
# Annotation = variable (+ Remote community if applicable)

# From deseq results, extract columns of interest
temp <- genus_decontaminated_group_comparison.df[c("Group_1", "Group_2", "Variable","Taxonomy", "log2FoldChange", "padj")]

# Create unique label for group comparison and the corresponding variable
temp$column_label <- with(temp, paste0(Group_1, " vs ", Group_2))
temp$unique_label <- with(temp, paste0(Group_1, " vs ", Group_2, " (",Variable,")"))
heatmap_metadata.df <- unique(temp[c("unique_label", "Variable","column_label")])

rownames(heatmap_metadata.df) <- heatmap_metadata.df$unique_label

# Create heatmap matrix
# unique_label vs Taxonomy, fill is log2FoldChange
heatmap.m <- temp[c("unique_label", "Taxonomy", "log2FoldChange")]
heatmap.m <- df2matrix(heatmap.m %>% spread(unique_label, log2FoldChange,fill = 0))

heatmap_metadata.df[c("unique_label", "column_label")]
source("code/helper_functions.R")
temp <- make_heatmap(myheatmap_matrix = heatmap.m, 
             mymetadata = heatmap_metadata.df,
             filename = "test.pdf",
             variables = c("Variable"),
             cluster_columns = F,
             row_title = "Genus",
             plot_height = 3,
             plot_width = 10,
             # legend_labels = seq(-30, 30, by =5),
             my_breaks = seq(-30, 30, by =5),
             discrete_legend = T,
             my_palette = c("darkred", "white","royalblue"),
             my_col_labels = heatmap_metadata.df[c("unique_label", "column_label")],
             my_annotation_palette = my_colour_palette_10_distinct,
             annotation_name_size = 0,
             # column_split = heatmap_metadata.df$unique_label
             )
# TODO option to exclude labels so heatmaps can be combined
temp2 <- temp
test <- c(temp$heatmap, temp2$heatmap)
test
draw(object = test[[1]] + test[[2]])
draw(object = temp$heatmap + temp2$heatmap, annotation_legend_list = list(temp$legend))

heatmap.m[, rownames(heatmap_metadata.df)] 
colnames(heatmap.m[, rownames(subset(heatmap_metadata.df, Variable == myvar)),drop = F])

my_heatmaps <- c()
for (myvar in unique(heatmap_metadata.df$Variable)){
  temp <- make_heatmap(myheatmap_matrix = heatmap.m, 
               mymetadata = subset(heatmap_metadata.df, Variable == myvar),
               # filename = "test.pdf",
               variables = c("Variable"),
               cluster_columns = F,
               row_title = "Genus",
               plot_height = 3,
               plot_width = 10,
               # legend_labels = seq(-30, 30, by =5),
               my_breaks = seq(-30, 30, by =5),
               discrete_legend = T,
               my_palette = c("darkred", "white","royalblue"),
               my_col_labels = heatmap_metadata.df[c("unique_label", "column_label")],
               my_annotation_palette = my_colour_palette_10_distinct,
               annotation_name_size = 0,
  )
  my_heatmaps <- c(my_heatmaps, temp$heatmap)
}
draw(lapply(my_heatmaps, function(x) x$heatmap))

HeatmapAnnotation()
Heatmap(column_split = T, show_column_dend = T)

# Columns of heatmap matrix must match rows of metadata, hence a unique label is required
# For deseq results, all we want to annotated is the variable each group belongs to. This requires us to map each label to the variable
hm <- Heatmap(heatmap.m,
        cluster_rows = T,
        cluster_columns = F,
        clustering_method_columns = "average",
        clustering_method_rows = "average",
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 6),
        show_column_dend = F,
        show_row_dend = F)
draw(hm)
#draw(hm, annotation_legend_list = c(hm_legend))
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------











# 
# # outfilename <- paste("Result_tables/DESeq_results/OTU_variable_groups.csv", sep= "")
# # write.csv(combined_results_ordered.df, file=outfilename, quote = F, row.names = F)
# 
# # Compare groups for all variables
# combined_results_ordered.df <- data.frame()
# for (myvar in discrete_variables){
#   metadata_filtered.df <- metadata_decontaminated.df[!is.na(metadata_decontaminated.df[,myvar]),]
#   otu_decontaminated_filtered.m <- otu_decontaminated.m[,rownames(metadata_filtered.df)]
#     
#   dds <- DESeqDataSetFromMatrix(countData = otu_decontaminated_filtered.m, colData = metadata_filtered.df, design = as.formula(paste0("~", myvar)))
#   geoMeans <- apply(counts(dds), 1, gm_mean)
#   dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
#   dds <- try(DESeq(dds, test = "Wald", fitType = "parametric", parallel = T))
#   group_combinations <- combn(sort(unique(metadata_filtered.df[,myvar])),2)
#   for (i in 1:ncol(group_combinations)){
#     group_1 <- as.character(group_combinations[1,i])
#     group_2 <- as.character(group_combinations[2,i])
# 
#     # group_1_meta <- subset(full, get(myvar) == group_1)
#     # group_2_meta <- subset(full, get(myvar) == group_2)
#     # n_group_1 <- dim(group_1_meta)[1]
#     # n_group_2 <- dim(group_2_meta)[1]
#     
#     n_group_1 <- dim(subset(metadata_filtered.df, get(myvar) == group_1))[1]
#     n_group_2 <- dim(subset(metadata_filtered.df, get(myvar) == group_2))[1]
#     
#     resMFSource <- results(dds, contrast = c(myvar,group_1,group_2), alpha=0.01, independentFiltering = F, cooksCutoff = F,parallel = T)
#     resMFSource$Group_1 <- group_1
#     resMFSource$Group_2 <- group_2
#     resMFSource$Variable <- myvar
#     resMFSource$N_Group_1 <- n_group_1
#     resMFSource$N_Group_2 <- n_group_2
#     resMFSource$Taxonomy <- assign_taxonomy_to_otu(resMFSource, otu_taxonomy_map.df)
#     resMFSource <- m2df(resMFSource, "OTU")
#     resMFSource <- filter_and_sort_dds_results(resMFSource, 0.01)
#     combined_results_ordered.df <- rbind(combined_results_ordered.df, resMFSource)
#   }
# }
# outfilename <- paste("Result_tables/DESeq_results/OTU_variable_groups.csv", sep= "")
# write.csv(combined_results_ordered.df, file=outfilename, quote = F, row.names = F)
# 
# 
# group_combinations
# grepl("Remote_Community", "Remote_Community_")
# combn(unique(metadata_decontaminated.df[,"Gold_Star"]),2)
# combn(unique(subset(metadata_decontaminated.df, Remote_Community == "0")[,"Gold_Star"]),2)
# combn(unique(subset(metadata_decontaminated.df, Remote_Community == "1")[,"Gold_Star"]),2)
# 
# # Compare groups for all variables within each community
# combined_results_ordered.df <- data.frame()
# for (community in unique(metadata_decontaminated.df$Remote_Community)){
#   for (myvar in discrete_variables){
#     # for (myvar in c("Gold_Star")){
#     if (myvar == "Remote_Community") {next}
#     if (grepl("Remote_Community", myvar)) {next}
#     metadata_filtered.df <- subset(metadata_decontaminated.df, Remote_Community == community)
#     metadata_filtered.df <- metadata_filtered.df[!is.na(metadata_filtered.df[,myvar]),]
#     otu_decontaminated_filtered.m <- otu_decontaminated.m[,rownames(metadata_filtered.df)]
#     
#     metadata_filtered.df[,myvar] <- factor(metadata_filtered.df[,myvar], levels = sort(unique(as.character(metadata_filtered.df[,myvar]))))
#     
#     dds <- DESeqDataSetFromMatrix(countData = otu_decontaminated_filtered.m, colData = metadata_filtered.df, design = as.formula(paste0("~", myvar)))
#     geoMeans <- apply(counts(dds), 1, gm_mean)
#     dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
#     dds <- try(DESeq(dds, test = "Wald", fitType = "parametric", parallel = T))
#     group_combinations <- combn(sort(unique(metadata_filtered.df[,myvar])),2)
#     for (i in 1:ncol(group_combinations)){
#       group_1 <- as.character(group_combinations[1,i])
#       group_2 <- as.character(group_combinations[2,i])
#       
#       n_group_1 <- dim(subset(metadata_filtered.df, get(myvar) == group_1))[1]
#       n_group_2 <- dim(subset(metadata_filtered.df, get(myvar) == group_2))[1]
#       
#       resMFSource <- results(dds, contrast = c(myvar,group_1,group_2), alpha=0.01, independentFiltering = F, cooksCutoff = F,parallel = T)
#       resMFSource$Remote_Community <- community
#       resMFSource$Group_1 <- group_1
#       resMFSource$Group_2 <- group_2
#       resMFSource$Variable <- myvar
#       resMFSource$N_Group_1 <- n_group_1
#       resMFSource$N_Group_2 <- n_group_2
#       resMFSource$Taxonomy <- assign_taxonomy_to_otu(resMFSource, otu_taxonomy_map.df)
#       resMFSource <- m2df(resMFSource, "OTU")
#       resMFSource <- filter_and_sort_dds_results(resMFSource, 0.01)
#       combined_results_ordered.df <- rbind(combined_results_ordered.df, resMFSource)
#     }
#   }
# }
# 
# outfilename <- paste("Result_tables/DESeq_results/OTU_variable_groups_within_community.csv", sep= "")
# write.csv(combined_results_ordered.df, file=outfilename, quote = F, row.names = F)
