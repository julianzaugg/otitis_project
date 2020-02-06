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


common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
  panel.background = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white", size = 1),
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


# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T, row.names = "Sequence_file_ID_clean")
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv", sep =",", header = T,row.names = "Sequence_file_ID_clean")

# Load counts
otu_decontaminated.m <-  as.matrix(read.table("Result_tables/count_tables/OTU_counts_decontaminated.csv", sep =",", header =T, row.names = 1))
genus_decontaminated.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts_decontaminated.csv", sep =",", header =T, row.names = 1))

# Load combined data
otu_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata_decontaminated.csv", header = T)
genus_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata_decontaminated.csv", header = T)

# discrete_variables <- c("Remote_Community","Gold_Star","OM_6mo","Season","Nose","OM_Classification", "Remote_Community_Season",
#                         "Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae",
#                         "Remote_Community_OM_Classification")
discrete_variables_to_add_with_counts <- c("Remote_Community","Gold_Star","OM_6mo",
                                           "Season","Nose","OM_Classification")

# Generate taxonomy summary
otu_decontaminated_taxa_summary.df <- generate_taxa_summary(mydata = otu_data_decontaminated.df, taxa_column = "OTU.ID")
genus_decontaminated_taxa_summary.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus")
# nose_taxa_summary.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus", group_by_columns = c("Sample","Nose"))
# write.csv(genus_decontaminated_taxa_summary.df,file = "genus_nose_summary.csv", quote = F, row.names = F)
# write.csv(nose_taxa_summary.df,file = "genus_nose_summary_per_sample.csv", quote = F, row.names = F)

# Transform read counts
otu_decontaminated_clr.m <-  clr(otu_decontaminated.m)
genus_decontaminated_clr.m <-  clr(genus_decontaminated.m)

# Add metadata to count matrix (optional)
genus_decontaminated_clr_with_meta.m <- rbind(genus_decontaminated_clr.m, t(metadata_decontaminated.df[colnames(genus_decontaminated_clr.m),discrete_variables_to_add_with_counts]))
# genus_decontaminated.m <- rbind(genus_decontaminated.m, Nose = metadata.df[colnames(genus_decontaminated.m),]$Nose)

# plot_feature_correlations(genus_decontaminated_clr_with_meta.m, feature = "Nose",top_n = 25)
# calculate_feature_correlations(genus_decontaminated_clr_with_meta.m, feature = "Nose")

# Filter to samples in groups of interest

# Filter by prevalence
otu_decontaminated_taxa_summary_filtered.df <- otu_decontaminated_taxa_summary.df %>% filter(Percent_group_samples > 40)
genus_decontaminated_taxa_summary_filtered.df <- genus_decontaminated_taxa_summary.df %>% filter(Percent_group_samples > 40)

# Filter the matrix
otu_decontaminated_filtered.m <- otu_decontaminated.m[unique(otu_decontaminated_taxa_summary_filtered.df$OTU.ID),]
otu_decontaminated_clr_filtered.m <- otu_decontaminated_clr.m[unique(otu_decontaminated_taxa_summary_filtered.df$OTU.ID),]

genus_decontaminated_filtered.m <- genus_decontaminated.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),]
genus_decontaminated_clr_filtered.m <- genus_decontaminated_clr.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),]
# genus_decontaminated_clr_with_meta_filtered.m <- genus_decontaminated_clr_with_meta.m[c(as.character(unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus)),discrete_variables_to_add_with_counts),]


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


# correlation_network.l <- generate_correlation_network(cor_matrix = otu_decontaminated_fastspar_cor.m,
#                                                       p_matrix = otu_decontaminated_fastspar_pval.m,
#                                                       p_value_threshold = 0.05,
#                                                       cor_threshold = .3)
# correlation_network.l$network_plot


genus_relabeller_function <- function(my_labels){
  unlist(lapply(my_labels, 
                function(x) {
                  phylostring <- unlist(strsplit(x, split = ";"))
                  # paste(phylostring[2],phylostring[3], phylostring[6], sep = ";")
                  paste(phylostring[3], phylostring[6], sep = ";")
                }))
}
colnames(genus_decontaminated_fastspar_cor_filtered.m) <- genus_relabeller_function(colnames(genus_decontaminated_fastspar_cor_filtered.m))
rownames(genus_decontaminated_fastspar_cor_filtered.m) <- genus_relabeller_function(rownames(genus_decontaminated_fastspar_cor_filtered.m))

correlation_network.l <- generate_correlation_network(cor_matrix = genus_decontaminated_fastspar_cor_filtered.m,
                                                      p_matrix = genus_decontaminated_fastspar_pval_filtered.m,
                                                      p_value_threshold = 0.001,
                                                      cor_threshold = 0.05)
correlation_network.l$network_plot

plot_corrplot(genus_decontaminated_fastspar_cor_filtered.m)

corrplot(corr = genus_decontaminated_fastspar_cor.m,
         outline = T,
         tl.col = "black",
         tl.cex = 0.1,
         col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", 
                                      "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                      "#4393C3", "#2166AC", "#053061")))(200),
         type = "lower",
         diag = F,
         order = "hclust",
         hclust.method = "average",
         p.mat = genus_decontaminated_fastspar_pval.m,
         sig.level = 0.05,
         insig = "blank",
         pch = 4,
         pch.cex = .1,
         pch.col = "black",
         cl.pos = 'r')




spiec_genus <- spiec.easi(t(genus_decontaminated.m), method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))


d <- ncol(t(genus_decontaminated.m))
n <- nrow(t(genus_decontaminated.m))
e <- d
graph <- SpiecEasi::make_graph('cluster', d, e)

huge::huge.roc(spiec_genus$est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(spiec_genus), graph, verbose=FALSE)

spiec_genus <- spiec.easi(t(genus_decontaminated.m), method='mb', lambda.min.ratio=1e-2,
                            nlambda=20, pulsar.params=list(rep.num=50)) 
sparcc_genus <- sparcc()


  
  
  
  
  
# Generate dataframe for graph generation
correlation_results <- calculate_correlation_matrix(genus_decontaminated_clr_filtered.m)
cor.m <- correlation_results$cor_matrix
cor_pval.m <- correlation_results$cor_pval_matrix
graph.df <- melt(cor.m)
graph.df$pvalue <- melt(cor_pval.m)$value
names(graph.df) <- c("Variable_1", "Variable_2", "Correlation", "P_value")
graph.df <- graph.df[graph.df$P_value <= 0.001,]
graph.df <- graph.df[abs(graph.df$Correlation) >= .3,]

# Generate graph object and remove looped edges and isolated nodes
# browseVignettes("ggraph") 
graph.df <- as_tbl_graph(graph.df) %>%
  # Remove loops
  activate(edges) %>%
  filter(!edge_is_loop()) %>%
  # Remove isolated nodes
  activate(nodes) %>%
  filter(!node_is_isolated())

ggraph(graph.df, layout = 'kk') +
  # geom_edge_link() +
  geom_edge_link(aes(colour = Correlation), show.legend = T, width = .5, alpha = 1) +
  scale_edge_color_gradient2(low = "darkblue", high = "red", mid = "lightyellow", limits = c(-1,1)) +
  geom_node_point(aes(size = P_value))

correlation_graph_plot <- ggraph(graph.df, layout = 'kk') +
  geom_edge_link(aes(colour = Correlation), show.legend = T, width = .5, alpha = 1) +
  scale_edge_color_gradient2(low = "darkblue", high = "red", mid = "lightyellow", limits = c(-1,1)) +
  geom_node_point(colour = "black", fill = "darkolivegreen3", pch = 21,size =6) +
  geom_node_text(aes(label = name), size = 2, nudge_y = -.01,repel = T, fontface = "bold") +
  ggtitle("") +
  theme_graph()
correlation_graph_plot

# geom_node_point(aes(shape=association, size=importance, colour=genus)) +
#   geom_edge_link(aes(colour = corCLR))
# 
# ggraph(gr, 'circlepack', weight=genus) +
#   geom_node_circle(aes(shape=association, size=importance, colour=genus)) +
#   geom_edge_link(aes(colour = corCLR))
# , weight = 'size'





genus_decontaminated_filtered.m <- genus_decontaminated.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),]
genus_decontaminated_clr_filtered.m <- genus_decontaminated_clr.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),]
# correlation_results <- calculate_correlation_matrix(genus_decontaminated_clr_filtered.m, method = "pearson", adjust = "BH")

spiec_genus <- spiec.easi(t(genus_decontaminated.m), method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50)) 
sparcc_genus <- sparcc(t(genus_decontaminated.m))

colnames(sparcc_genus$Cor) <- colnames(t(genus_decontaminated.m))
rownames(sparcc_genus$Cor) <- colnames(t(genus_decontaminated.m))
calculate_correlation_matrix_stats(sparcc_genus$Cor)

cor.test(sparcc_genus$Cor)

genus_decontaminated_clr.m <-  clr(genus_decontaminated.m) # transform
genus_decontaminated_clr_cor.m <- cor(t(genus_decontaminated_clr.m), method = "pearson") # generate correlation matrix
genus_decontaminated_clr_cor_test.m <- corr.test(t(genus_decontaminated_clr.m),adjust = "BH",method = "pearson") # calculate adjusted p-values

# filter to significant comparisons
genus_decontaminated_clr_cor_lower.m <- genus_decontaminated_clr_cor.m
genus_decontaminated_clr_cor_lower.m[which(upper.tri(x = genus_decontaminated_clr_cor_lower.m, diag = T))] <- NA
inds.select.df <- which(!is.na(genus_decontaminated_clr_cor_lower.m) & genus_decontaminated_clr_cor_test.m$p < 1e-5, arr.ind=TRUE)

# Taxa (row / columns) that are significant
rnames.select = rownames(genus_decontaminated_clr_cor_lower.m)[inds.select.df[,1]]
cnames.select = colnames(genus_decontaminated_clr_cor_lower.m)[inds.select.df[,2]]

# Correlation table
corr_network.df <- data.frame(row.names = paste(rnames.select,cnames.select))
corr_network.df[,c('Genus_1','Genus_2')] <- cbind(rnames.select,cnames.select)
corr_network.df[,'correlation_CLR'] <- round(genus_decontaminated_clr_cor_lower.m[inds.select.df], digits = 2)
corr_network.df[,'pvalue_adjusted_CLR'] <- genus_decontaminated_clr_cor_test.m$p[inds.select.df]
write.csv(corr_network.df, file = "genus_correlations_pvalues_nose.csv", quote = F, row.names = F)


# significant_genus <- unique(c(corr_network.df[abs(corr_network.df$corCLR) > 0,]$Genus_1, 
                              # corr_network.df[abs(corr_network.df$corCLR) > 0,]$Genus_2))
significant_genus <- unique(corr_network.df$Genus_1, corr_network.df$Genus_2)

pdf("correlations_nose.pdf", height = 4)
plot_correlations(genus_clr_cor.m[significant_genus, significant_genus],50)
dev.off()
plot_correlations(genus_clr_cor.m,50)

temp <- subset(genus_data_decontaminated.df, taxonomy_genus %in% rownames(genus_clr_cor.m) & Nose %in% c(0,1))
ggplot(temp, aes(x = factor(Nose), y = Read_count_logged)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(size = .4) +
  facet_wrap(~taxonomy_genus)

# ************************************************************

Heatmap(genus_clr_cor.m,show_column_names = F,show_row_names = F)

spiec_genus <- spiec.easi(t(genus.m), method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50)) 
sparcc_genus <- sparcc()

cor.mat <- cor(t(genus.m),method='spear')

# Limit to samples in the same 

sparcc_genus_graph <- abs(sparcc_genus$Cor) >= 0.3
diag(sparcc_genus_graph) <- 0
sparcc_genus_graph <- Matrix(sparcc_genus_graph, sparse=TRUE)

ig.mb     <- adj2igraph(getRefit(spiec_genus))
ig.sparcc <- adj2igraph(sparcc_genus_graph)

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(t(genus.m), 1), 1)+6
am.coord <- layout.fruchterman.reingold(ig.mb)


par(mfrow=c(1,2))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")





# Abundance
for (myvar in discrete_variables){
  for (taxa in unique(genus_data.df$taxonomy_genus)){
    genus_data.df
    break
  }  
}

# group_combinations <- combn(as.character(unique(mydata[,variable])), 2)
# 
# for (i in 1:ncol(group_combinations)) {
#   group_1 <- group_combinations[1,i]
#   group_2 <- group_combinations[2,i]
#   group_1_meta <- subset(mydata, get(variable) == group_1)
#   group_2_meta <- subset(mydata, get(variable) == group_2)
#   if (is.na(group_1) | is.na(group_2)) {next}
#   if (dim(group_1_meta)[1] == 0 | dim(group_2_meta)[1] == 0) {next}
#   N_samples_group_1 <- dim(group_1_meta)[1]
#   N_samples_group_2 <- dim(group_2_meta)[1]
#   