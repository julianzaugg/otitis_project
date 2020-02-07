# Calculate the correlations on metadata

library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
# library(igraph)
# install.packages("psych")
library(psych)
# install.packages("ggraph")
library(ggraph)
library(tidygraph)

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")


# Load the data
data.m <- as.matrix(read.table("data/Culture_Viral_qPCR_binary_results.tsv", sep ="\t", header = T,row.names = 1))

metadata.df <- read.table("data/metadata.tsv", sep ="\t", header = T)
metadata.df$Sample_ID_2 <- gsub("-", ".", metadata.df$Sample_ID)
rownames(metadata.df) <- metadata.df$Sample_ID_2
data.m <- rbind(data.m, Gold_Star = metadata.df[colnames(data.m),]$Gold_Star)

# data.m["N_RSV_B",]
# data.m["N_FLU_B",]
# data.m["N_FLU_A",]

zv <- apply(data.m, 1, function(x) length(unique(x)) == 1)
data.m <- data.m[!zv, ]

# ?logit()

# temp <- calculate_feature_correlations(data.m, feature = "Gold_Star",method = "spearman")
# head(temp)
# temp <- calculate_feature_correlations(data.m, feature = "Gold_Star",method = "pearson")
# head(temp)
# plot_correlations(data.m, feature = "Gold_Star",top_n = 25)
# max(temp$`p-value`[temp$pval_adj_BH < 0.05])
# max(temp$`p-value`)
# max(temp$`pval_adj_BH`)
# dim(data.m)
# dim(temp)
# rownames(data.m)[!rownames(data.m) %in% rownames(temp)]

# plot_feature_correlations(data.m, feature = "Bacillus_cereus",top_n = 25, method = "pearson")

# Calculate correlation results and generate correlation plots for each variable
for (feature in rownames(data.m)){
  feature_filename <- gsub("/_", "_", feature)
  feature_filename <- gsub("/", "_", feature_filename)
  calculate_feature_correlations(data.m,
                                 feature = feature,
                                 filename = paste0("culture_paper_analysis/results/",
                                                   feature_filename, "_correlation_results.csv"),
                                 method = "pearson")
  
  pdf(paste0("culture_paper_analysis/results/", feature_filename, "_correlations.pdf"),height=6,width=6)
  plot_feature_correlations(data.m, feature = feature,top_n = 25, method = "pearson")
  dev.off()
}


# Calculate correlation matrix and p-values
correlation_results <- calculate_correlation_matrix(data.m, method = "pearson", adjust = "BH")
correlation_results <- calculate_correlation_matrix(data.m, method = "pearson", adjust = "BH")

correlation_results_g0 <- calculate_correlation_matrix(data.m[,data.m["Gold_Star",] == 0], method = "pearson", adjust = "BH")
correlation_results_g1 <- calculate_correlation_matrix(data.m[,data.m["Gold_Star",] == 1], method = "pearson", adjust = "BH")

cor.m <- correlation_results$cor_matrix
cor_pval.m <- correlation_results$cor_pval_matrix

cor_g0.m <- correlation_results_g0$cor_matrix
cor_pval_g0.m <- correlation_results_g0$cor_pval_matrix

cor_g1.m <- correlation_results_g1$cor_matrix
cor_pval_g1.m <- correlation_results_g1$cor_pval_matrix

source("code/helper_functions.R")

Mcat <- as.data.frame(t(combn(grep("Mcat", rownames(data.m),value =T),2)))
H_inf <- as.data.frame(t(combn(grep("H.inf", rownames(data.m),value =T),2)))
Spn <- as.data.frame(t(combn(grep("Spn", rownames(data.m),value =T),2)))

edges_to_remove.df <- rbind(Mcat, H_inf,Spn)

correlation_network.l <- generate_correlation_network(cor_matrix = cor.m,
                                                      p_matrix = cor_pval.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.4,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      plot_title = "Correlations between culture variables (p-value <= 0.05, |correlation| >= 0.4)",
                                                      edge_width_min = .5,
                                                      edge_width_max = 1.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/correlation_graph.pdf")

correlation_network.l <- generate_correlation_network(cor_matrix = cor_g0.m,
                                                      p_matrix = cor_pval_g0.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.4,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      plot_title = "Correlations between culture variables, Gold star 0 (p-value <= 0.05, |correlation| >= 0.4)",
                                                      edge_width_min = .5,
                                                      edge_width_max = 1.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/correlation_graph_gold_star_0.pdf")


# cor_g1.m["H.influnezae_>3rd_IQR"]
# plot_feature_correlations(data.m[,data.m["Gold_Star",] == 1], "H.influnezae_>3rd_IQR")

correlation_network.l <- generate_correlation_network(cor_matrix = cor_g1.m,
                                                      p_matrix = cor_pval_g1.m,
                                                      p_value_threshold = 0.05,
                                                      cor_threshold = 0.4,
                                                      node_size = 4,
                                                      node_colour = "grey20",
                                                      node_fill = "grey20",
                                                      label_colour = "black",
                                                      label_size = 4,
                                                      plot_height = 10,
                                                      plot_width = 10,
                                                      plot_title = "Correlations between culture variables, Gold star 1 (p-value <= 0.05, |correlation| >= 0.4)",
                                                      edge_width_min = .5,
                                                      edge_width_max = 1.5,
                                                      network_layout = "fr",
                                                      exclude_to_from_df = edges_to_remove.df,
                                                      filename="culture_paper_analysis/results/correlation_graph_gold_star_1.pdf")



library(corrplot)
pdf(file = "culture_paper_analysis/corrplot_metadata_v2.pdf",width = 8,height =8)
corrplot(corr = cor.m,
         outline = T,
         tl.col = "black",
         tl.cex = .3,
         col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", 
                                   "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                   "#4393C3", "#2166AC", "#053061")))(200),
         type = "lower",
         diag = F,
         # order = "FPC",
         order = "hclust",
         hclust.method = "average",
         p.mat = cor_pval.m,
         sig.level = 0.05,
         insig = "blank",
         # addgrid.col = NA,
         # insig = "pch",
         # insig = "p-value",
         # insig = "n",
         pch = 4,
         pch.cex = .1,
         pch.col = "black",
         cl.pos = 'r')
dev.off()
# Generate dataframe for graph generation
graph.df <- melt(cor_g1.m)
graph.df$pvalue <- melt(cor_pval_g1.m)$value

names(graph.df) <- c("Variable_1", "Variable_2", "Correlation", "P_value")
graph.df <- graph.df[graph.df$P_value <= 0.01,]
graph.df <- graph.df[abs(graph.df$Correlation) >= .4,]
write.csv(graph.df,file = "culture_paper_analysis/results/correlation_data_network.csv", quote =F, row.names = F)


# Generate graph object and remove looped edges and isolated nodes
# browseVignettes("ggraph") 
graph.df <- as_tbl_graph(graph.df) %>%
  # Remove loops
  activate(edges) %>%
  filter(!edge_is_loop()) %>%
  # Remove isolated nodes
  activate(nodes) %>%
  filter(!node_is_isolated())
ggraph(graph.df, layout = "kk") +
  # geom_edge_fan(aes(colour = Correlation), show.legend = T, width = .8, alpha = 1) +
  # geom_edge_elbow(aes(colour = Correlation), show.legend = T, width = .8, alpha = 1,strength = 1) +
  # geom_edge_bend(aes(colour = Correlation), show.legend = T, width = .8, alpha = 1,strength = 1) +
  # geom_edge_hive(aes(colour = Correlation), show.legend = T, width = .8, alpha = 1,strength = 1) +
  geom_node_point(colour = "black", fill = "grey20", pch = 21,size = 4) +
  scale_edge_colour_gradientn(colours = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                                               "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                               "#4393C3", "#2166AC", "#053061")))(200),
                              limits = c(-1,1), breaks = seq(-1,1,.2),
                              guide = guide_edge_colourbar(barwidth = 0.5, barheight = 10))

correlation_graph_plot <- ggraph(graph.df, layout = "kk") +
  geom_edge_link(aes(colour = Correlation), show.legend = T, width = .8, alpha = 1) +
  scale_edge_color_gradient2(low = "darkblue", high = "red", mid = "lightyellow", limits = c(-1,1)) +
  scale_edge_size_continuous(range = c(-1,.1,1)) +
  geom_node_point(colour = "black", fill = "darkolivegreen3", pch = 21,size =6) +
  geom_node_text(aes(label = name), size = 2, nudge_y = -.01,repel = T, fontface = "bold") +
  ggtitle("Gold star 1, p-value <= 0.01 and |correlation| >= 0.4") +
  theme_graph()
correlation_graph_plot
# pdf(file = "culture_paper_analysis/results/correlation_graph_plot.pdf", width = 20,height = 20)
# correlation_graph_plot
# dev.off()

ggsave(filename = "culture_paper_analysis/results/correlation_graph_plot_additional_filtering_gold_star_1.jpeg", 
       plot = correlation_graph_plot, width = 20, height = 20,units = "cm",dpi = 300)



# net.grph <- graph.adjacency(correlation_results$cor_matrix, mode="undirected", weighted=TRUE, diag = F)
# net.grph <- igraph::simplify(net.grph,remove.multiple=TRUE, remove.loops=TRUE)
# 
# temp <- melt(correlation_results$p_values)
# 
# E(net.grph)
# # Colour negative correlation edges as blue
# E(net.grph)[which(E(net.grph)$weight<0)]$color <- "royalblue"
# 
# # Colour positive correlation edges as red
# E(net.grph)[which(E(net.grph)$weight>0)]$color <- "darkred"
# 
# # Convert edge weights to absolute values
# E(net.grph)$weight <- abs(E(net.grph)$weight)
# 
# # Remove edges below absolute Pearson correlation 0.8
# net.grph <- delete_edges(net.grph, E(net.grph)[which(E(net.grph)$weight<0.4)])
# 
# 
# # Remove any vertices remaining that have no edges
# net.grph <- delete.vertices(net.grph, degree(net.grph)==0)
# 
# # Assign names to the graph vertices (optional)
# V(net.grph)$name <- V(net.grph)$name
# 
# # Change shape of graph vertices
# V(net.grph)$shape <- "circle"
# 
# # Change colour of graph vertices
# V(net.grph)$color <- "steelblue"
# 
# # Change colour of vertex frames
# V(net.grph)$vertex.frame.color <- "white"
# 
# # Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# # Multiply scaled vales by a factor of 10
# # scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # vSizes <- (scale01(apply(estrogenMainEffects, 1, mean)) + 1.0) * 10
# 
# # Amplify or decrease the width of the edges
# edgeweights <- E(net.grph)$weight * 10.0
# 
# # Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
# # mst <- mst(net.grph, algorithm="prim")
# 
# # Plot the tree object
# plot(
#   net.grph,
#   layout=layout.fruchterman.reingold,
#   edge.curved=F,
#   # vertex.size=vSizes,
#   vertex.label.dist=-0.5,
#   vertex.label.color="black",
#   vertex.label.cex=.8,
#   # asp=FALSE,
#   edge.width=edgeweights,
#   edge.arrow.mode=0,
#   main="Correlations between bacterial species/viral/pathobiont load variables >= 0.3"
#   )


# edgew <- E(net.grph)$weight
# bad.vs <- V(net.grph)[degree(net.grph) == 0] 
# net.grph <- delete.vertices(net.grph, bad.vs)

# pdf("otitis_meta_correlation.pdf",height=15,width=15)

# Heatmap(cor.m,na_col = "grey", cluster_columns = T,cluster_rows = T,show_column_names = T,show_row_names = T,show_row_dend = F, show_column_dend = F,
#         row_names_gp = gpar(cex = .5), column_names_gp = gpar(cex = .5),
#         clustering_method_columns = "average",
#         clustering_method_rows = "average",
#         heatmap_legend_param = gpar(title = "Pearson correlation"),
#         # cell_fun = function(j, i, x, y, width, height, fill){
#         #   if(!is.na(cor.m[i, j])){
#         #     if(cor.m[i, j] > 0 & cor.m[i, j] != 1){
#         #       grid.text(sprintf("%.2f", cor.m[i, j]), x, y, gp = gpar(fontsize = 6))}          
#         #   }
#         # }
#         )
# dev.off()
