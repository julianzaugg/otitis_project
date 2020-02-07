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


setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T, row.names = "Sequence_file_ID_clean")
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv", sep =",", header = T,row.names = "Sequence_file_ID_clean")

# Load count matrices
otu_decontaminated.m <-  as.matrix(read.table("Result_tables/count_tables/OTU_counts_decontaminated.csv", sep =",", header =T, row.names = 1))
genus_decontaminated.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts_decontaminated.csv", sep =",", header =T, row.names = 1))

# Load combined data (counts, abundances and metadata)
otu_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata_decontaminated.csv", header = T)
genus_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata_decontaminated.csv", header = T)

discrete_variables <- c("Remote_Community","Gold_Star","OM_6mo","Season","Nose","OM_Classification", "Remote_Community_Season",
                        "Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae",
                        "Remote_Community_OM_Classification")
discrete_variables_to_add_with_counts <- c("Remote_Community","Gold_Star","OM_6mo",
                                           "Season","Nose","OM_Classification")

# Generate taxonomy summaries
otu_decontaminated_taxa_summary.df <- generate_taxa_summary(mydata = otu_data_decontaminated.df, taxa_column = "OTU.ID")
genus_decontaminated_taxa_summary.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus")
genus_decontaminated_taxa_summary_nose.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus", group_by_columns = c("Nose"))
genus_decontaminated_taxa_summary_gold_star.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus", group_by_columns = c("Gold_Star"))
# nose_taxa_summary.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus", group_by_columns = c("Sample","Nose"))
# write.csv(genus_decontaminated_taxa_summary.df,file = "genus_nose_summary.csv", quote = F, row.names = F)
# write.csv(nose_taxa_summary.df,file = "genus_nose_summary_per_sample.csv", quote = F, row.names = F)

# Filter taxa summaries by prevalence
otu_decontaminated_taxa_summary_filtered.df <- otu_decontaminated_taxa_summary.df %>% filter(Percent_group_samples > 20)
genus_decontaminated_taxa_summary_filtered.df <- genus_decontaminated_taxa_summary.df %>% filter(Percent_group_samples > 20)
genus_decontaminated_taxa_summary_nose_filtered.df <- genus_decontaminated_taxa_summary_nose.df %>% filter(Percent_group_samples > 20)
genus_decontaminated_taxa_summary_gold_star_filtered.df <- genus_decontaminated_taxa_summary_gold_star.df %>% filter(Percent_group_samples > 20)

# Transform read counts
otu_decontaminated_clr.m <-  clr(otu_decontaminated.m)
genus_decontaminated_clr.m <-  clr(genus_decontaminated.m)

# Add metadata to count matrix (optional)
otu_decontaminated_clr_with_meta.m <- rbind(otu_decontaminated_clr.m, t(metadata_decontaminated.df[colnames(otu_decontaminated_clr.m),discrete_variables_to_add_with_counts]))
genus_decontaminated_clr_with_meta.m <- rbind(genus_decontaminated_clr.m, t(metadata_decontaminated.df[colnames(genus_decontaminated_clr.m),discrete_variables_to_add_with_counts]))

# Filter to samples in groups of interest

# Filter taxa, e.g to taxa above a certain prevalence
otu_decontaminated_filtered.m <- otu_decontaminated.m[unique(otu_decontaminated_taxa_summary_filtered.df$OTU.ID),]
otu_decontaminated_clr_filtered.m <- otu_decontaminated_clr.m[unique(otu_decontaminated_taxa_summary_filtered.df$OTU.ID),]

genus_decontaminated_filtered.m <- genus_decontaminated.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),]
genus_decontaminated_clr_filtered.m <- genus_decontaminated_clr.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),]

genus_decontaminated_clr_filtered_nose.m <- genus_decontaminated_clr.m[unique(genus_decontaminated_taxa_summary_nose_filtered.df$taxonomy_genus),]
genus_decontaminated_clr_filtered_gold_star.m <- genus_decontaminated_clr.m[unique(genus_decontaminated_taxa_summary_gold_star_filtered.df$taxonomy_genus),]

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# Run Fastspar externally on the raw counts (not the filtered matrices above) and then load the results.
otu_decontaminated_fastspar_cor.m <- as.matrix(read.table("Additional_results/OTU_decontaminated_correlation.tsv",
                                                            sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
otu_decontaminated_fastspar_pval.m <- as.matrix(read.table("Additional_results/OTU_decontaminated_pvalues.tsv",
                                                             sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))

otu_decontaminated_fastspar_cor_filtered.m  <- otu_decontaminated_fastspar_cor.m[otu_decontaminated_taxa_summary_filtered.df$OTU.ID,
                                                                                 otu_decontaminated_taxa_summary_filtered.df$OTU.ID]

otu_decontaminated_fastspar_pval_filtered.m  <- otu_decontaminated_fastspar_pval.m[otu_decontaminated_taxa_summary_filtered.df$OTU.ID,
                                                                                   otu_decontaminated_taxa_summary_filtered.df$OTU.ID]

genus_decontaminated_fastspar_cor.m <- as.matrix(read.table("Additional_results/Genus_decontaminated_correlation.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_decontaminated_fastspar_pval.m <- as.matrix(read.table("Additional_results/Genus_decontaminated_pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))



genus_decontaminated_fastspar_cor_filtered.m  <- genus_decontaminated_fastspar_cor.m[genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus,
                                                                                     genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus]

genus_decontaminated_fastspar_pval_filtered.m  <- genus_decontaminated_fastspar_pval.m[genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus,
                                                                                     genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus]

genus_decontaminated_fastspar_cor_nose_filtered.m  <- genus_decontaminated_fastspar_cor.m[genus_decontaminated_taxa_summary_nose_filtered.df$taxonomy_genus,
                                                                                          genus_decontaminated_taxa_summary_nose_filtered.df$taxonomy_genus]

genus_decontaminated_fastspar_pval_nose_filtered.m  <- genus_decontaminated_fastspar_pval.m[genus_decontaminated_taxa_summary_nose_filtered.df$taxonomy_genus,
                                                                                            genus_decontaminated_taxa_summary_nose_filtered.df$taxonomy_genus]

genus_decontaminated_fastspar_cor_gold_star_filtered.m  <- genus_decontaminated_fastspar_cor.m[genus_decontaminated_taxa_summary_gold_star_filtered.df$taxonomy_genus,
                                                                                     genus_decontaminated_taxa_summary_gold_star_filtered.df$taxonomy_genus]

genus_decontaminated_fastspar_pval_gold_star_filtered.m  <- genus_decontaminated_fastspar_pval.m[genus_decontaminated_taxa_summary_gold_star_filtered.df$taxonomy_genus,
                                                                                       genus_decontaminated_taxa_summary_gold_star_filtered.df$taxonomy_genus]

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# Generate correlation network plots
genus_relabeller_function <- function(my_labels){
  unlist(lapply(my_labels, 
                function(x) {
                  phylostring <- unlist(strsplit(x, split = ";"))
                  # paste(phylostring[2],phylostring[3], phylostring[6], sep = ";")
                  # paste(phylostring[3], phylostring[6], sep = ";")
                  paste(phylostring[3], phylostring[5], phylostring[6], sep = ";")
                }))
}

# relabel column and row names to shorter form (have to be unique!)
# colnames(genus_decontaminated_fastspar_cor_filtered.m) <- genus_relabeller_function(colnames(genus_decontaminated_fastspar_cor_filtered.m))
# rownames(genus_decontaminated_fastspar_cor_filtered.m) <- genus_relabeller_function(rownames(genus_decontaminated_fastspar_cor_filtered.m))
# genus_correlation_network.l <- generate_correlation_network(genus_decontaminated_fastspar_cor_filtered.m,
                                                            # relabeller_function = genus_relabeller_function,
                                                            # cor_threshold = 0.1)
# genus_correlation_network.l$network_plot
# generate_correlation_network(genus_decontaminated_fastspar_cor_filtered.m)
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_decontaminated_fastspar_cor_filtered.m,
                                                      p_matrix = genus_decontaminated_fastspar_pval_filtered.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.4,
                                                      relabeller_function = genus_relabeller_function,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 2,
                                                      plot_height = 8,
                                                      plot_width = 8,
                                                      plot_title = "",
                                                      filename="Result_other/correlation_analysis/networks/genus_decontaminated_filtered_fastspar_cor_network.pdf")
genus_correlation_network.l$network_plot

genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_decontaminated_fastspar_cor_gold_star_filtered.m,
                                                            p_matrix = genus_decontaminated_fastspar_pval_gold_star_filtered.m,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.4,
                                                            relabeller_function = genus_relabeller_function,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            label_colour = "black",
                                                            label_size = 2,
                                                            plot_height = 8,
                                                            plot_width = 8,
                                                            plot_title = "",
                                                            filename="Result_other/correlation_analysis/networks/genus_decontaminated_gold_star_filtered_fastspar_cor_network.pdf")
genus_correlation_network.l$network_plot
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_decontaminated_fastspar_cor_nose_filtered.m,
                                                            p_matrix = genus_decontaminated_fastspar_pval_nose_filtered.m,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.4,
                                                            relabeller_function = genus_relabeller_function,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            label_colour = "black",
                                                            label_size = 2,
                                                            plot_height = 8,
                                                            plot_width = 8,
                                                            plot_title = "",
                                                            filename="Result_other/correlation_analysis/networks/genus_decontaminated_nose_filtered_fastspar_cor_network.pdf")

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# Corrplots
plot_corrplot(correlation_matrix = genus_decontaminated_fastspar_cor_filtered.m,
              p_value_matrix = genus_decontaminated_fastspar_pval_filtered.m,
              p_value_threshold = 0.05,
              label_size = .3,
              relabeller_function = genus_relabeller_function,
              filename = "Result_other/correlation_analysis/corrplots/genus_decontaminated_filtered_fastspar_corrplot.pdf")

plot_corrplot(correlation_matrix = genus_decontaminated_fastspar_cor_nose_filtered.m,
              p_value_matrix = genus_decontaminated_fastspar_pval_nose_filtered.m,
              p_value_threshold = 0.05,
              label_size = .3,
              relabeller_function = genus_relabeller_function,
              filename = "Result_other/correlation_analysis/corrplots/genus_decontaminated_filtered_nose_fastspar_corrplot.pdf")

plot_corrplot(correlation_matrix = genus_decontaminated_fastspar_cor_gold_star_filtered.m,
              p_value_matrix = genus_decontaminated_fastspar_pval_gold_star_filtered.m,
              p_value_threshold = 0.05,
              label_size = .3,
              relabeller_function = genus_relabeller_function,
              filename = "Result_other/correlation_analysis/corrplots/genus_decontaminated_gold_star_filtered_fastspar_corrplot.pdf")

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# By feature
myfeature <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Cardiobacteriales;f__Cardiobacteriaceae;g__Dichelobacter"

source("code/helper_functions.R")
plot_feature_correlations_external(cor_matrix = genus_decontaminated_fastspar_cor_filtered.m,
                                   p_value_matrix = genus_decontaminated_fastspar_pval_filtered.m,
                                   feature = myfeature,
                                   top_n = 25,
                                   plot_width = 10, plot_height = 10)

# for (feature in colnames(genus_decontaminated_fastspar_cor_filtered.m)){
#   # TODO shorten feature name for filename, also remove ;
# }




# spiec_genus <- spiec.easi(t(genus_decontaminated.m), method='mb', lambda.min.ratio=1e-2,
#                           nlambda=20, pulsar.params=list(rep.num=50))
# 
# 
# d <- ncol(t(genus_decontaminated.m))
# n <- nrow(t(genus_decontaminated.m))
# e <- d
# graph <- SpiecEasi::make_graph('cluster', d, e)
# 
# huge::huge.roc(spiec_genus$est$path, graph, verbose=FALSE)
# stars.pr(getOptMerge(spiec_genus), graph, verbose=FALSE)
# 
# spiec_genus <- spiec.easi(t(genus_decontaminated.m), method='mb', lambda.min.ratio=1e-2,
#                             nlambda=20, pulsar.params=list(rep.num=50)) 


  
