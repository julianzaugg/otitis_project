library(DESeq2)
# library(BiocParallel)




############################################################
######################## Functions #########################

# Function that calculates the geometric mean with some error-protection bits. 
# DESeq2 does not appear to work (will throw an error) if every OTU (or genus or genome etc.) 
# contains at least one count of zero in every row of the count data.
# Specifically, the function "dds<-DESeq(dds, betaPrior = FALSE)" will fail
# One way to address this is to use the function below as input to DESeq2 to transform the data.
# Calculate the geometric means prior to estimating the size factors
gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

# Center log ratio transform
clr = function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

matrix2df <- function(mymatrix, column_name){
  out <- as.data.frame(mymatrix)
  out_names <- colnames(out)
  out$placeholder <- rownames(out)
  rownames(out) <- NULL
  names(out)[length(names(out))] <- column_name
  out <- out[,c(column_name, out_names)]
  return(out)
}

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


# Load count table at the OTU level. These are the counts for OTUs that were above our abundance thresholds
otu_rare.m <- as.matrix(read.table("Result_tables/count_tables/OTU_counts_rarefied.csv", sep =",", header =T, row.names = 1))
otu_genus_rare.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts_rarefied.csv", sep =",", header =T, row.names = 1))



# Filter out OTUs/species/genera that do not have at # reads in at least one sample
head(melt(sort(colSums(otu_rare.m))))
otu_rare.m <- filter_matrix_rows(otu_rare.m,50)
head(melt(sort(colSums(otu_rare.m))))
# otu_genus_rare.m <- filter_matrix_rows(otu_genus_rare.m,0)


# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)



# Only keep columns (samples) in the metadata
otu_rare.m <- otu_rare.m[,colnames(otu_rare.m) %in% as.character(metadata.df$Index)]
# otu_genus_rare.m <- otu_genus_rare.m[,colnames(otu_genus_rare.m) %in% as.character(metadata.df$Index)]


# Order the metadata.df by the index value
metadata.df <- metadata.df[order(metadata.df$Index),]


# Since we likely removed samples from the count matrix
# in the main script, remove them from the metadata.df here
samples_removed <- metadata.df$Index[!metadata.df$Index %in% colnames(otu_rare.m)]
metadata.df <- metadata.df[! metadata.df$Index %in% samples_removed,]

# Order the metadata.df by the index value
metadata.df <- metadata.df[order(metadata.df$Index),]

# Rownames should match the sample columns in the otu table
rownames(metadata.df) <- metadata.df$Index

# Order the otu_tables the same order as the metadata
otu_rare.m <- otu_rare.m[,rownames(metadata.df)]
# otu_genus_rare.m <- otu_genus_rare.m[,rownames(metadata.df)]

dim(otu_rare.m)
dim(metadata.df)


# ---------------------------------------------------------------------------------------------------------
# Perform differential abundance calculations at the OTU level and genus level, 
# comparing between the groups within variables of interest

# DESeq requires a count matrix ('countData'), a corresponding metadata.df ('colData') and a 'design' formula. The formula expresses
# how the counts for each OTU/genus depend on the variables defined in the 'colData'. See help(DESeqDataSetFromMatrix) for more information.
# The first column of the metadata.df ('colData') must match the ordering of the columns of the countData

# Ensure names of the otu / genus count matrices match the order of the metadata.df!
# Assumes number of samples in metadata.df and count data are the same
all(colnames(otu_rare.m) == metadata.df$Index) # Should be 'True'
all(colnames(otu_rare.m) == rownames(metadata.df)) # Should be 'True'


# Define the discrete variables
discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season","Nose", 
                        "Otitis_status_OM_6mo","Remote_Community_Otitis_status","OM_6mo_Type_OM","Remote_Community_Season")

# Convert variables to factors
metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)
# for (dv in discrete_variables){
  # metadata.df[,dv]
# }
# dv
# metadata.df$Gold_Star

# factor(re)
# ?factor

# Compare groups for all variables
# as.formula(paste0(custom_formula)))
# design = as.formula(paste0("~", variable))) 
combined_results_ordered.df <- data.frame()
for (myvar in discrete_variables){
  metadata_filtered.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  otu_rare_filtered.m <- otu_rare.m[,rownames(metadata_filtered.df)]
    
  dds <- DESeqDataSetFromMatrix(countData = otu_rare_filtered.m, colData = metadata_filtered.df, design = as.formula(paste0("~", myvar)))
  geoMeans <- apply(counts(dds), 1, gm_mean)
  dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
  dds <- try(DESeq(dds, test = "Wald", fitType = "parametric", parallel = T))
  group_combinations <- combn(sort(unique(metadata_filtered.df[,myvar])),2)
  for (i in 1:ncol(group_combinations)){
    group_1 <- as.character(group_combinations[1,i])
    group_2 <- as.character(group_combinations[2,i])

    # group_1_meta <- subset(full, get(myvar) == group_1)
    # group_2_meta <- subset(full, get(myvar) == group_2)
    # n_group_1 <- dim(group_1_meta)[1]
    # n_group_2 <- dim(group_2_meta)[1]
    
    n_group_1 <- dim(subset(metadata_filtered.df, get(myvar) == group_1))[1]
    n_group_2 <- dim(subset(metadata_filtered.df, get(myvar) == group_2))[1]
    
    resMFSource <- results(dds, contrast = c(myvar,group_1,group_2), alpha=0.01, independentFiltering = F, cooksCutoff = F,parallel = T)
    resMFSource$Group_1 <- group_1
    resMFSource$Group_2 <- group_2
    resMFSource$Variable <- myvar
    resMFSource$N_Group_1 <- n_group_1
    resMFSource$N_Group_2 <- n_group_2
    resMFSource$Taxonomy <- assign_taxonomy_to_otu(resMFSource, otu_taxonomy_map.df)
    resMFSource <- matrix2df(resMFSource, "OTU")
    resMFSource <- filter_and_sort_dds_results(resMFSource)
    combined_results_ordered.df <- rbind(combined_results_ordered.df, resMFSource)
  }
}
outfilename <- paste("Result_tables/DESeq_results/OTU_variable_groups.csv", sep= "")
write.csv(combined_results_ordered.df, file=outfilename, quote = F, row.names = F)


group_combinations
grepl("Remote_Community", "Remote_Community_")
combn(unique(metadata.df[,"Gold_Star"]),2)
combn(unique(subset(metadata.df, Remote_Community == "0")[,"Gold_Star"]),2)
combn(unique(subset(metadata.df, Remote_Community == "1")[,"Gold_Star"]),2)

# Compare groups for all variables within each community
combined_results_ordered.df <- data.frame()
for (community in unique(metadata.df$Remote_Community)){
  for (myvar in discrete_variables){
    # for (myvar in c("Gold_Star")){
    if (myvar == "Remote_Community") {next}
    if (grepl("Remote_Community", myvar)) {next}
    metadata_filtered.df <- subset(metadata.df, Remote_Community == community)
    metadata_filtered.df <- metadata_filtered.df[!is.na(metadata_filtered.df[,myvar]),]
    otu_rare_filtered.m <- otu_rare.m[,rownames(metadata_filtered.df)]
    
    metadata_filtered.df[,myvar] <- factor(metadata_filtered.df[,myvar], levels = sort(unique(as.character(metadata_filtered.df[,myvar]))))
    
    dds <- DESeqDataSetFromMatrix(countData = otu_rare_filtered.m, colData = metadata_filtered.df, design = as.formula(paste0("~", myvar)))
    geoMeans <- apply(counts(dds), 1, gm_mean)
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    dds <- try(DESeq(dds, test = "Wald", fitType = "parametric", parallel = T))
    group_combinations <- combn(sort(unique(metadata_filtered.df[,myvar])),2)
    for (i in 1:ncol(group_combinations)){
      group_1 <- as.character(group_combinations[1,i])
      group_2 <- as.character(group_combinations[2,i])
      
      n_group_1 <- dim(subset(metadata_filtered.df, get(myvar) == group_1))[1]
      n_group_2 <- dim(subset(metadata_filtered.df, get(myvar) == group_2))[1]
      
      resMFSource <- results(dds, contrast = c(myvar,group_1,group_2), alpha=0.01, independentFiltering = F, cooksCutoff = F,parallel = T)
      resMFSource$Remote_Community <- community
      resMFSource$Group_1 <- group_1
      resMFSource$Group_2 <- group_2
      resMFSource$Variable <- myvar
      resMFSource$N_Group_1 <- n_group_1
      resMFSource$N_Group_2 <- n_group_2
      resMFSource$Taxonomy <- assign_taxonomy_to_otu(resMFSource, otu_taxonomy_map.df)
      resMFSource <- matrix2df(resMFSource, "OTU")
      resMFSource <- filter_and_sort_dds_results(resMFSource)
      combined_results_ordered.df <- rbind(combined_results_ordered.df, resMFSource)
    }
  }
}

outfilename <- paste("Result_tables/DESeq_results/OTU_variable_groups_within_community.csv", sep= "")
write.csv(combined_results_ordered.df, file=outfilename, quote = F, row.names = F)
