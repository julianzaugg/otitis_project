# Calculate the correlations in abundances between taxa

library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
library(igraph)
# install.packages("psych")
library(psych)


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


# ?? from david's code
robust_cor <- function(M, MARGIN=1) {
  robust_corA <- function(A) {
    apply(M, MARGIN = MARGIN, FUN = function(B) {
      if (length(which(A != 0 & B != 0)) > length(A)*0.2) {
        cor(A[!(A == 0 & B == 0)], B[!(A == 0 & B == 0)])
      } else {
        NA
      }
    })
  }
  apply(M, MARGIN=MARGIN, FUN = robust_corA)
}


# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv", sep =",", header = T)

genus.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts.csv", sep =",", header =T, row.names = 1))
genus_rel.m <-  as.matrix(read.table("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", sep =",", header =T, row.names = 1))

genus_decontaminated.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts_decontaminated.csv", sep =",", header =T, row.names = 1))

genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv", header = T)
genus_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata_decontaminated.csv", header = T)

# discrete_variables <- c("Remote_Community","Gold_Star","OM_6mo","Season","Nose","OM_Classification", "Remote_Community_Season",
#                         "Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae",
#                         "Remote_Community_OM_Classification")
# discrete_variables <- c("Nose")

# Generate summary
genus_decontaminated_taxa_summary.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus", group_by_columns = c("Nose"))
nose_taxa_summary.df <- generate_taxa_summary(mydata = genus_data_decontaminated.df, taxa_column = "taxonomy_genus", group_by_columns = c("Sample","Nose"))
# write.csv(genus_decontaminated_taxa_summary.df,file = "genus_nose_summary.csv", quote = F, row.names = F)
# write.csv(nose_taxa_summary.df,file = "genus_nose_summary_per_sample.csv", quote = F, row.names = F)

# colnames(genus_decontaminated.m)
# metadata_decontaminated.df
rownames(metadata.df) <- metadata.df$Sequence_file_ID_clean
genus_decontaminated_clr.m <-  clr(genus_decontaminated.m) # transform
genus_decontaminated_clr.m <- rbind(genus_decontaminated_clr.m, Nose = metadata.df[colnames(genus_decontaminated_clr.m),]$Nose)
genus_decontaminated.m <- rbind(genus_decontaminated.m, Nose = metadata.df[colnames(genus_decontaminated.m),]$Nose)
plot_correlations(genus_decontaminated_clr.m, feature = "Nose",top_n = 25)
calculate_feature_correlations(genus_decontaminated_clr.m, feature = "Nose")



# Filter to samples in groups of interest
genus_decontaminated.m <- genus_decontaminated.m[,metadata.df[metadata.df$Nose %in% c(0,1),]$Index]


# Filter by prevalence
genus_decontaminated_taxa_summary_filtered.df <- genus_decontaminated_taxa_summary.df %>% filter(Percent_group_samples > 40)

# Filter the matrix
dim(genus_decontaminated.m)
dim(genus_decontaminated.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),])
genus_decontaminated.m <- genus_decontaminated.m[unique(genus_decontaminated_taxa_summary_filtered.df$taxonomy_genus),]
# write.csv(genus_decontaminated.m, file = "genera_counts.csv",quote = F)
# temp <- metadata.df[,c("Index", "Nose"),drop= F]
# rownames(temp) <- temp$Index
# temp$Index <- NULL
# write.csv(temp[rownames(temp) %in% colnames(genus_decontaminated.m), ,drop = F],file = "genera_meta.csv",quote = F, row.names = T)

# NOTE - sparcc was provided counts in original study, not transformed
# however correlation networks were made from clr transformed values
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
sparcc_genus <- sparcc

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