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

# Remove entries that don't vary across all samples
# dim(data.m)
# zv <- apply(data.m, 1, function(x) length(unique(x)) == 1)
# data.m <- data.m[!zv, ]#(suppose df is the name of your dataset)
# dim(data.m)

# temp <- calculate_feature_correlations(data.m, feature = "Gold_Star")
# plot_correlations(data.m, feature = "Gold_Star",top_n = 25)

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
  plot_correlations(data.m, feature = feature,top_n = 25, method = "pearson")
  dev.off()
}

# Calculate correlation matrix and p-values
correlation_results <- calculate_correlation_matrix(data.m, method = "pearson", adjust = "BH")
cor.m <- correlation_results$cor_matrix
cor_pval.m <- correlation_results$cor_pval_matrix

# Generate dataframe for graph generation
graph.df <- melt(cor.m)
graph.df$pvalue <- melt(cor_pval.m)$value

names(graph.df) <- c("Variable_1", "Variable_2", "Correlation", "P_value")
graph.df <- graph.df[graph.df$P_value <= 0.05,]
# graph.df <- graph.df[abs(graph.df$Correlation) >= .3,]
write.csv(graph.df,file = "culture_paper_analysis/results/correlation_data_network.csv", quote =F, row.names = F)


# browseVignettes("ggraph") 
graph.df <- as_tbl_graph(graph.df) %>%
  # Remove loops
  activate(edges) %>%
  filter(!edge_is_loop()) %>%
  # Remove isolated nodes
  activate(nodes) %>%
  filter(!node_is_isolated())


correlation_graph_plot <- ggraph(graph.df, layout = "kk") +
  geom_edge_link(aes(colour = Correlation), show.legend = T, width = 1, alpha = 1) +
  # scale_edge_color_continuous(low = "darkblue", high = "red", limits = c(-1,1)) +
  scale_edge_color_gradient2(low = "darkblue", high = "red", mid = "lightyellow", limits = c(-1,1)) +
  # scale_edge_color_gradient2(low = "darkblue", high = "red", mid = "white", limits = c(-1,1)) +
  # scale_edge_colour_manual(values = c("TRUE" = "red", "FALSE" = "royalblue"), name = "test") +
  # geom_point(aes(x = x, y = y), colour = "black", fill = "steelblue", pch = 21,size =6) 
  geom_node_point(colour = "black", fill = "darkolivegreen3", pch = 21,size =6) +
  geom_node_text(aes(label = name), size = 2, nudge_y = -.01,repel = T, fontface = "bold") +
  theme_graph()
# pdf(file = "culture_paper_analysis/results/correlation_graph_plot.pdf", width = 20,height = 20)
# correlation_graph_plot
# dev.off()

ggsave(filename = "culture_paper_analysis/results/correlation_graph_plot.jpeg", 
       plot = correlation_graph_plot, width = 20, height = 20,units = "cm",dpi = 300)

# data("highschool")
# graph <- graph_from_data_frame(highschool)

# Not specifying the layout - defaults to "auto"
# ggraph(graph,layout = "kk",maxiter  = 100) + 
#   geom_edge_link(aes(colour = factor(year))) + 
#   geom_node_point()





# calculate_feature_correlations(data.m, feature = "Mcat_>3rd_IQR")
# plot_correlations(data.m, feature = "Mcat_>3rd_IQR",top_n = 25)

# cor.m <- cor(t(data.m), method = "pearson",use="complete.obs") # generate correlation matrix
# cor_test.m <- corr.test(t(data.m),adjust = "BH",method = "pearson") # calculate adjusted p-values
# diag(cor.m) <- NA


# inds.select.df <- which(!is.na(correlation_results$cor_matrix) & 
#                           abs(correlation_results$cor_matrix) > 0.4 & 
#                           correlation_results$p_values < 0.05, arr.ind=TRUE)
# rnames.select = rownames(correlation_results$cor_matrix)[inds.select.df[,1]]
# cnames.select = colnames(correlation_results$cor_matrix)[inds.select.df[,2]]



net.grph <- graph.adjacency(correlation_results$cor_matrix, mode="undirected", weighted=TRUE, diag = F)
net.grph <- igraph::simplify(net.grph,remove.multiple=TRUE, remove.loops=TRUE)

temp <- melt(correlation_results$p_values)

E(net.grph)
# Colour negative correlation edges as blue
E(net.grph)[which(E(net.grph)$weight<0)]$color <- "royalblue"

# Colour positive correlation edges as red
E(net.grph)[which(E(net.grph)$weight>0)]$color <- "darkred"

# Convert edge weights to absolute values
E(net.grph)$weight <- abs(E(net.grph)$weight)

# Remove edges below absolute Pearson correlation 0.8
net.grph <- delete_edges(net.grph, E(net.grph)[which(E(net.grph)$weight<0.4)])


# Remove any vertices remaining that have no edges
net.grph <- delete.vertices(net.grph, degree(net.grph)==0)

# Assign names to the graph vertices (optional)
V(net.grph)$name <- V(net.grph)$name

# Change shape of graph vertices
V(net.grph)$shape <- "circle"

# Change colour of graph vertices
V(net.grph)$color <- "steelblue"

# Change colour of vertex frames
V(net.grph)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
# scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
# vSizes <- (scale01(apply(estrogenMainEffects, 1, mean)) + 1.0) * 10

# Amplify or decrease the width of the edges
edgeweights <- E(net.grph)$weight * 10.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
# mst <- mst(net.grph, algorithm="prim")

# Plot the tree object
plot(
  net.grph,
  layout=layout.fruchterman.reingold,
  edge.curved=F,
  # vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  vertex.label.cex=.8,
  # asp=FALSE,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="Correlations between bacterial species/viral/pathobiont load variables >= 0.3"
  )


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
