# Calculate the correlations in abundances between taxa

library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
library(igraph)
# install.packages("psych")
library(psych)

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

#otu.m <- as.matrix(read.table("Result_tables/count_tables/OTU_counts.csv", sep =",", header =T, row.names = 1))
genus.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts.csv", sep =",", header =T, row.names = 1))
genus_rel.m <-  as.matrix(read.table("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", sep =",", header =T, row.names = 1))

#otu_decontaminated.m <- as.matrix(read.table("Result_tables/count_tables/OTU_counts_decontaminated.csv", sep =",", header =T, row.names = 1))
# genus_decontaminated.m <-  as.matrix(read.table("Result_tables/count_tables/Genus_counts_decontaminated.csv", sep =",", header =T, row.names = 1))
# genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv", header = T,)


discrete_variables <- c("Remote_Community","Gold_Star","OM_6mo","Season","Nose","OM_Classification", "Remote_Community_Season",
                        "Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae",
                        "Remote_Community_OM_Classification")
discrete_variables <- c("Nose")
genus.m <- genus.m[,metadata.df[metadata.df$Nose %in% c(0,1),]$Index]
genus_rel.m <- genus_rel.m[,metadata.df[metadata.df$Nose %in% c(0,1),]$Index]


# robust_cor(genus.m, MARGIN = 1)
Heatmap(cor(t(genus_rel.m),method = "pearson"),show_column_names = F,show_row_names = F)
corr.test(t(genus_rel.m),adjust = "BH")

spiec_genus <- spiec.easi(t(genus.m), method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50)) 
sparcc_genus <- sparcc 
dim(sparcc_genus$Cor)
Heatmap(sparcc_genus$Cor)
plot(density(sparcc_genus$Cor))

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