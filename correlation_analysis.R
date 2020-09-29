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


setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T, row.names = "Sequence_file_ID_clean")


# Define descrete variables
discrete_variables <- c("Community","Gold_Star","Tympanic_membrane","Otitis_Status", "Season","Nose","Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae","N_HRV")
discrete_variables_to_add_with_counts <- c("Community","Gold_Star","Season","Nose")

# Load feature taxonomy map
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)
rownames(otu_taxonomy_map.df) <- otu_taxonomy_map.df$OTU.ID

# Load count matrices
otu.m <-  as.matrix(read.table("Result_tables/count_tables/OTU_counts.csv", sep =",", header =T, row.names = 1))
genus.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts.csv", sep =",", header =T, row.names = 1))

# Load combined data (counts, abundances and metadata)
# otu_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata.csv", header = T)
# genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv", header = T)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
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


  
