library(vegan)
library(ggplot2)
library(ggfortify)



############################################################
# Various colour palettes
my_colour_palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
my_colour_palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
my_colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
my_colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
my_colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
my_colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
my_colour_palette_32_distinct <- c("#ea7e00","#ca0074","#d1c69b","#474007","#bb00ad","#9c80ff","#be3300","#542e72","#00b9f5","#09436b","#8b0036","#9ac8e6","#ff1059","#959eff","#154a11","#0290f4","#ff7762","#7dbf00","#ff8194","#834c00","#006e73","#f9bb5d","#d6c943","#017229","#00d3a8","#732427","#36e191","#6a8200","#efb3ea","#3227bb","#ff90e1","#e92a12")
# lesion_palette_7 <- c("#8558d6","#6ee268","#d247ad","#c9d743","#d7453e","#59a237","#d78f2a")
# patient_palette_45 <- c("#d64530","#585fb1","#795d97","#9e4773","#3f6921","#71692c","#a2b93c","#d571cc","#9b3e97","#33947a","#98ad66","#448a4e","#869ae0","#5ce7af","#e085a3","#dfdc87","#d19be2","#5cb735","#e38269","#3db6c0","#50b565","#50902c","#a98a2c","#dde84a","#db3d76","#5fe485","#7c8329","#b3e791","#6fe965","#5ebce9","#3c86c1","#2a6a45","#65b688","#6651d1","#af4ed3","#df872f","#56e4db","#737cea","#ac464b","#dd37b5","#995b2b","#daac6f","#92e2be","#a2e24b","#e0be3a")
patient_palette_270 <- c("#456c00","#77bbff","#75fb3f","#273300","#f5006f","#ac008b","#125700","#ffef77","#00278e","#3d005e","#d84100","#015686","#01289f","#ff8c4c","#0070b5","#8015cd","#feffd6","#02d8aa","#019cf4","#4f2f00","#bbffe9","#c52900","#1b002c","#a3ac00","#5d9200","#f29fff","#231500","#934cfc","#988a00","#002cbb","#ffeb36","#ffa758","#f1f9ff","#000045","#b4004b","#602900","#390048","#e6c400","#00ce3c","#ff7bd0","#8cff56","#e60051","#b89aff","#00474b","#d5fbff","#ff79c2","#1d0016","#00635d","#ff8e33","#992300","#ff6e91","#ffa081","#534a00","#61002d","#ffe1c1","#8c0072","#00405d","#89ffc6","#607500","#64ff6f","#002e52","#9b97ff","#b1ccff","#02c5cd","#5dffba","#beff45","#00112b","#b8ff74","#7f0027","#0074cd","#005c6f","#3f00af","#dd7900","#cced00","#77ffd6","#ffc5b5","#99ffb1","#01ea72","#f0ff88","#007f63","#abff9d","#391200","#003a74","#114b00","#0a0100","#ff5fff","#ffccfb","#00d6b7","#c7ff93","#1efa5a","#005437","#f6af00","#a60024","#ffb7e6","#ea0043","#c7ffbc","#72ab00","#789300","#585500","#c3ff14","#00f398","#ab4a00","#9b7600","#85e5ff","#006235","#130020","#006825","#ff735c","#007a7f","#02a3a4","#4856ff","#bf52ff","#00edbc","#a31bd6","#009642","#e93bee","#e400ae","#ffbdd2","#00cfc7","#f1ffaa","#009b7a","#dd00c9","#ff697d","#004a14","#ff72ac","#ff3d1f","#fffaa3","#5d0400","#027ba4","#01c774","#002655","#00941f","#0a32d7","#82acff","#ff8af3","#ff4165","#001104","#ffd6f2","#efebff","#aebc00","#3e0030","#c5abff","#00402e","#ff4bae","#0275f1","#be89ff","#ffd899","#00c765","#01b208","#97ffd4","#7e9fff","#00fde1","#0050c9","#ff8eb5","#c800cd","#005173","#ff2b95","#76ff7a","#ea0074","#001d70","#009856","#f100a8","#ba6b00","#0293df","#00462d","#ff6862","#f6ff65","#02bbda","#2c2200","#01a876","#e35a00","#e3000f","#ff819e","#5a0039","#a558ff","#e2ffb2","#784800","#016def","#b400a2","#00143c","#00212d","#403d00","#ff75fe","#975300","#166c00","#260008","#917fff","#ff8d89","#01bf7a","#ffa6bf","#800086","#90a100","#cce4ff","#dad800","#52c900","#46a700","#0c0039","#0b0052","#79009d","#003c85","#bb0034","#59e7ff","#af0064","#64001e","#c0007e","#000897","#bd8400","#2b007f","#318400","#31f1ff","#7c8600","#807300","#ffc072","#6f005f","#770040","#e62c00","#2e0019","#005599","#6535e1","#5b0099","#006bd5","#0142a1","#baaf00","#00ab2d","#ffcc40","#edffec","#ef0031","#153100","#abe9ff","#6bbd00","#e5ff4e","#ffdb43","#ffa5ef","#01c4f3","#ffbd8f","#84d100","#bbff84","#9fcdff","#7b0059","#ffe897","#ff8711","#ffa869","#febaff","#20003a","#94002b","#5387ff","#756dff","#fff494","#a5c1ff","#e0ffcf","#002417","#530076","#ff8459","#ffe4ec","#00b650","#0119b7","#c963ff","#a2ff64","#9c6800","#03b6f8","#00a0c2","#00240b","#6297ff","#bd0010","#fff7af","#7d2d00","#cf7aff","#af5600","#322c00","#500028")
my_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")


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

# For each rowname (OTU), get the corresponding taxonomy_species
# Assumes "OTU.ID" and "taxonomy_species" columns in the provided map dataframe
assign_taxonomy_to_otu <- function(otutable, taxon_map){
  taxonomies <- c()
  for (otuid in rownames(otutable)){
    taxonomies <- c(taxonomies, as.character(taxon_map[taxon_map$OTU.ID == otuid,]$taxonomy_species))
  }
  return(taxonomies)
}

# Function that takes the metadata and a list of variables (column names) and returns those samples (rownames) with NA entries
get_samples_missing_data <- function(my_metadata, variables){
  samples_missing_data <- c()
  for (name in variables) {
    samples_missing_data <- c(samples_missing_data, rownames(my_metadata[is.na(my_metadata[[name]]),]))
  }
  return(unique(samples_missing_data))
}

bin_my_variable <- function(mydataframe, variable, breaks){
  temp <- cut(mydataframe[,variable], breaks = breaks,dig.lab = 4)
  temp <- gsub("\\(", "", temp)
  temp <- gsub("\\]", "", temp)
  temp <- gsub(",", "-", temp)
  return(temp)
}
############################################################

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project")
source("Code/helper_functions.R")

# Load count table at the OTU level. These are the counts for OTUs that were above our abundance thresholds
# otu_rare.df <- read.table("Result_tables/count_tables/OTU_counts_rarefied.csv", sep =",", header =T)
otu_decontaminated.df <- read.table("Result_tables/count_tables/OTU_counts_decontaminated.csv", sep =",", header =T)

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata_decontaminated.df) <- metadata_decontaminated.df$Index

# Since we likely removed samples from the count matrix
# in the main script, remove them from the metadata.df here
# samples_removed <- metadata.df$Index[!metadata.df$Index %in% names(otu_rare.df)]
# metadata.df <- metadata.df[! metadata.df$Index %in% samples_removed,]

# Remove samples from the OTU table that are not in the filtered metadata
# otu_rare.df <- otu_rare.df[,names(otu_rare.df) %in% c("OTU.ID", as.character(metadata.df$Index))]
otu_decontaminated.df <- otu_decontaminated.df[,names(otu_decontaminated.df) %in% c("OTU.ID", as.character(metadata_decontaminated.df$Index))]

# Remove samples from metadata that are not in the data
metadata_decontaminated.df <- metadata_decontaminated.df[metadata_decontaminated.df$Index %in% colnames(otu_decontaminated.df),]

# Order the metadata.df by the index value
metadata_decontaminated.df <- metadata_decontaminated.df[order(metadata_decontaminated.df$Index),]

# Create matrices
# otu_rare.m <- otu_rare.df
# rownames(otu_rare.m) <- otu_rare.df$OTU.ID
# otu_rare.m$OTU.ID <- NULL
# otu_rare.m <- as.matrix(otu_rare.m)
otu_decontaminated.m <- df2matrix(otu_decontaminated.df)

# Filter by reads per sample if you don't want to use the existing filtering
minimum_reads <- 0
# otu_rare.m <- otu_rare.m[,colSums(otu_rare.m) >= minimum_reads]
# otu_decontaminated.m <- otu_decontaminated.m[,colSums(otu_decontaminated.m) >= minimum_reads]
# metadata.df <- metadata.df[rownames(metadata.df) %in% colnames(otu_rare.m),]
metadata_decontaminated.df <- metadata_decontaminated.df[rownames(metadata_decontaminated.df) %in% colnames(otu_decontaminated.m),]

# Order the matrices and metadata to be the same order
metadata_decontaminated.df <- metadata_decontaminated.df[order(rownames(metadata_decontaminated.df)),]
# otu_rare.m <- otu_rare.m[,order(rownames(metadata.df))]
otu_decontaminated.m <- otu_decontaminated.m[,order(rownames(metadata_decontaminated.df))]


# Filter out OTUs that do not have at # reads in at least one sample
# dim(otu_rare.m)
# otu_rare_filtered.m <- otu_rare.m[apply(otu_rare.m,1,max) >= 50,]
# dim(otu_rare_filtered.m)

# dim(otu.m)
# otu_filtered.m <- otu_rare.m[apply(otu.m,1,max) >= 50,]
# dim(otu_filtered.m)

# CLR transform the otu matrix.
otu_clr_decontaminated.m <- clr(otu_decontaminated.m)
# otu_genus_clr.m <- clr(otu_genus.m)
# otu_class_clr.m <- clr(otu_class.m)

# If there are negative values, assign them a value of zero
# otu_rare_clr_filtered.m[which(otu_rare_clr_filtered.m < 0)] <- 0
# otu_clr_filtered.m[which(otu_clr_filtered.m < 0)] <- 0


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# Ordination analysis


discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season",
                        "Nose","Otitis_status_OM_6mo", "Remote_Community_Otitis_status", "OM_6mo_Type_OM","Remote_Community_Season")

# ------------------------
# pca_full_data <- rda(t(otu_rare_clr_filtered.m), data = metadata.df) # ~1 makes it unconstrained

# pca.scores <- scores(pca_full_data, choices=c(1,2,3),scaling = "symmetric")
# # Get component x,y coordinates
# pca_site_scores <- scores(pca_full_data, display = "sites")
# pca_specie_scores <- scores(pca_full_data, display = "species")
# pca_percentages <- (pca_full_data$CA$eig/sum(pca_full_data$CA$eig)) * 100
# 
# plot(pca_full_data,
#      type='n')
# grid(NULL,NULL, lty = 2, col = "grey80")
# points(pca.scores$sites)
# points(pca.scores$species, col = 'red',cex = .5, pch = "x")
# arrows(0,0,pca_specie_scores[,1],pca_specie_scores[,2],length =.5)
# 
# 
# pca.scores <- scores(pca_full_data)
# allOTUs <- pca.scores$species
# left_pca.v <- names(sort(pca.scores$species[,'PC2']))[1:3]
# right_pca.v <- names(sort(pca.scores$species[,'PC2'], decreasing = T))[1:3]
# left_rda.v <- names(sort(pca.scores$species[,'PC1']))[1:3]
# right_rda.v <- names(sort(pca.scores$species[,'PC1'], decreasing = T))[1:3]
# top_vars.v <- unique(c(left_pca.v, right_pca.v, left_rda.v, right_rda.v))
# # otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID %in% top_vars.v,]$taxonomy_species
# # pca_full_data$CA$v <- pca_full_data$CA$v[top_vars.v,]
# top_OTUs <- pca_full_data$CA$v[top_vars.v,]
# 
# lda.arrows <- function(x, myscale = 1, tex = 0.75, choices = c(1,2), ...){
#   ## adds `biplot` arrows to an lda using the discriminant function values
#   heads <- coef(x)
#   arrows(x0 = 0, y0 = 0, 
#          x1 = myscale * heads[,choices[1]], 
#          y1 = myscale * heads[,choices[2]], ...)
#   text(myscale * heads[,choices], labels = row.names(heads), 
#        cex = tex)
# }
# pca_full_data$CA
# 
# arrows(0,0,top_OTUs[,1],top_OTUs[,2], length = .5,angle =30,col = 'red')
# text(pca_full_data, display="bp", scaling=1)
# biplot(pca_full_data,)
# 
# plot(c(0:10),type="n")
# 
# arrows(1,0,2,1,length=0.2)
# arrows(1,1,2,2,length=0.1,angle=40,lwd=3)
# 
# 
# text(x = pca_specie_scores[,1],
#      y = pca_specie_scores[,2],
#      labels = rownames(pca_specie_scores),
#      cex = .5,
#      pos = 2)
# arrows(x = pca_specie_scores[,1],
#        y = pca_specie_scores[,2])
# arrows

# data.frame(pca_full_data$CA$v)
#biplot(pca_full_data)
# pca_full_data <- capscale(t(otu_rare.m)~1,distance="bray", data = metadata.df)
# --------------------------------------------------------------------------------------------------------
# https://stats.stackexchange.com/questions/276645/arrows-of-underlying-variables-in-pca-biplot-in-r
# X <- t(otu_rare_clr_filtered.m)
# 
# CEN = scale(X, center = T, scale = T) # Centering and scaling the data
# PCA = prcomp(CEN)
# 
# # EIGENVECTORS:
# (evecs.ei = eigen(cor(CEN))$vectors)       # Using eigen() method
# (evecs.svd = svd(CEN)$v)                   # PCA with SVD...
# (evecs = prcomp(CEN)$rotation)             # Confirming with prcomp()
# 
# # EIGENVALUES:
# (evals.ei = eigen(cor(CEN))$values)        # Using the eigen() method
# (evals.svd = svd(CEN)$d^2/(nrow(X) - 1))   # and SVD: sing.values^2/n - 1
# (evals = prcomp(CEN)$sdev^2)               # with prcomp() (needs squaring)
# 
# # SCORES:
# scr.svd = svd(CEN)$u %*% diag(svd(CEN)$d)  # with SVD
# scr = prcomp(CEN)$x                        # with prcomp()
# scr.mm = CEN %*% prcomp(CEN)$rotation      # "Manually" [data] [eigvecs]
# 
# # LOADINGS:
# 
# loaded = evecs %*% diag(prcomp(CEN)$sdev)  # [E-vectors] [sqrt(E-values)]
# 
# arrows(0, 0,
#        cor(X[,1], scr[,1]) * 0.8 * sqrt(nrow(X) - 1), 
#        cor(X[,1], scr[,2]) * 0.8 * sqrt(nrow(X) - 1), 
#        lwd = 1, angle = 30, length = 0.1, col = 4)
# arrows(0, 0,
#        cor(X[,1], scr[,1]) * 0.8 * sqrt(nrow(X) - 1), 
#        cor(X[,1], scr[,2]) * 0.8 * sqrt(nrow(X) - 1), 
#        lwd = 1, angle = 30, length = 0.1, col = 4)
# --------------------------------------------------------------------------------------------------------
otu_pca <- rda(t(otu_clr_decontaminated.m), data = metadata_decontaminated.df)

source("Code/helper_functions.R")
rownames(otu_pca$CA$u)[!rownames(otu_pca$CA$u) %in% metadata_decontaminated.df$Index]
"PGDCPos"


my_relabeller_function <- function(my_labels){
  unlist(lapply(my_labels, 
                function(x) {
                  phylostring <- unlist(strsplit(x, split = ";"))
                  # paste(phylostring[2],phylostring[3], phylostring[6], sep = ";")
                  paste(phylostring[3], phylostring[6], sep = ";")
                }))
}

generate_pca(otu_pca, mymetadata = metadata_decontaminated.df,
             plot_height = 5, plot_width = 5,
             legend_x = -9, legend_y = 3,
             # legend_x = -2, legend_y = 2,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
             legend_title = "Remote_Community",
             legend_cex = .5,
             plot_title = "",
             limits = c(-9,4,-4,3),
             # limits = c(-2,2,-2,2),
             plot_spiders = F,
             plot_ellipses = F,
             plot_hulls = F,
             use_shapes = T,
             ellipse_border_width = .5,
             include_legend = T,
             label_ellipse = F, ellipse_label_size = .3,
             colour_palette = my_colour_palette_15,
             variable_to_plot = "Remote_Community", legend_cols = 1,
             variable_colours_available = T,
             num_top_species = 3,
             plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 2,arrow_thickness = .5,
             label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
             specie_labeller_function = my_relabeller_function,arrow_label_offset = 0,)





for (myvar in discrete_variables){
  generate_pca(pca_full_data, mymetadata = metadata.df,
               plot_height = 5, plot_width =5,
               legend_x = -7, legend_y = 6,
               point_size = 1, point_line_thickness = .3,point_alpha =.9,
               legend_title = myvar,
               plot_title = "",
               limits = c(-7,7,-8,7),
               plot_spiders = F,
               plot_ellipses = F,
               use_shapes = T,
               ellipse_border_width = .5,
               label_ellipse = F, ellipse_label_size = .5,
               color_palette = my_colour_palette_10_distinct,
               variable_to_plot = myvar, legend_cols = 1,
               variable_colours_available = T,
               filename = paste0("Result_figures/pca_plots/",myvar, "_pca.pdf"))
}

# ----------------------------
# Split by community

community_0_meta.df <- metadata.df[metadata.df$Remote_Community == 0,]
community_1_meta.df <- metadata.df[metadata.df$Remote_Community == 1,]

otu_rare_clr_community_0.m <- otu_rare_clr_filtered.m[,colnames(otu_rare_clr_filtered.m) %in% community_0_meta.df$Index]
otu_rare_clr_community_1.m <- otu_rare_clr_filtered.m[,colnames(otu_rare_clr_filtered.m) %in% community_1_meta.df$Index]

pca_community_0 <- rda(t(otu_rare_clr_community_0.m), data = community_0_meta.df) # ~1 makes it unconstrained
pca_community_1 <- rda(t(otu_rare_clr_community_1.m), data = community_1_meta.df) # ~1 makes it unconstrained


for (myvar in discrete_variables){
  if (myvar == "Remote_Community") {next}
  generate_pca(pca_community_0, mymetadata = community_0_meta.df,
               plot_height = 5, plot_width =5,
               legend_x = -7, legend_y = 6,
               point_size = 1, point_line_thickness = .3,point_alpha =.9,
               legend_title = myvar,
               plot_title = "",
               limits = c(-7,7,-8,7),
               plot_spiders = F,
               plot_ellipses = F,
               use_shapes = T,
               ellipse_border_width = .5,
               label_ellipse = F, ellipse_label_size = .5,
               color_palette = my_colour_palette_10_distinct,
               variable_to_plot = myvar, legend_cols = 1,
               variable_colours_available = T,
               filename = paste0("Result_figures/pca_plots/community_0__",myvar, "_pca.pdf"))
}


for (myvar in discrete_variables){
  if (myvar == "Remote_Community") {next}
  generate_pca(pca_community_1, mymetadata = community_1_meta.df,
               plot_height = 5, plot_width =5,
               legend_x = -7, legend_y = 6,
               point_size = 1, point_line_thickness = .3,point_alpha =.9,
               legend_title = myvar,
               plot_title = "",
               limits = c(-7,7,-8,7),
               plot_spiders = F,
               plot_ellipses = F,
               use_shapes = T,
               ellipse_border_width = .5,
               label_ellipse = F, ellipse_label_size = .5,
               color_palette = my_colour_palette_10_distinct,
               variable_to_plot = myvar, legend_cols = 1,
               variable_colours_available = T,
               filename = paste0("Result_figures/pca_plots/community_1__",myvar, "_pca.pdf"))
}

# ----------------------------------------
# Now split by each group for each variable and colour by Remote Community

for (myvar in discrete_variables){
  if (myvar == "Remote_Community") {next}
  for (group in unique(metadata.df[,myvar])){
    if (is.na(group)) {next}
    meta_subset.df <- subset(metadata.df, get(myvar) == group)
    if(dim(meta_subset.df)[1] < 3){next}
    if (length(unique(meta_subset.df$Remote_Community)) == 1) {next}
    otu_rare_clr_subset.m <- otu_rare_clr_filtered.m[,colnames(otu_rare_clr_filtered.m) %in% meta_subset.df$Index]
    pca_subset <- rda(t(otu_rare_clr_subset.m), data = meta_subset.df) # ~1 makes it unconstrained
    generate_pca(pca_subset, mymetadata = meta_subset.df,
                 plot_height = 5, plot_width =5,
                 legend_x = -11, legend_y =10,
                 point_size = 1, point_line_thickness = .3,point_alpha =.9,
                 legend_title = "Remote Community",
                 plot_title = paste0(myvar,", Group ", group),
                 limits = c(-10,10,-10,10),
                 plot_spiders = F,
                 plot_ellipses = F,
                 use_shapes = T,
                 ellipse_border_width = .5,
                 label_ellipse = F, ellipse_label_size = .5,
                 color_palette = my_colour_palette_10_distinct,
                 variable_to_plot = "Remote_Community", legend_cols = 1,
                 variable_colours_available = T,
                 filename = paste0("Result_figures/pca_plots/", myvar, "_", group, "__Remote_Community_pca.pdf"))
  }
}



# [1] "Sequence_file_ID"       "Sequence_file_ID_clean" "Sample_ID"              "Remote_Community"      
# [5] "Otitis_status"          "Gold_Star"              "OM_6mo"                 "Type_OM"               
# [9] "Season"                 "Nose"                   "Index"                  "Sample_retained"   



# ------------------------------------------------------------------------
# PERMANOVA analysis

# Permutational Multivariate Analysis of Variance (PERMANOVA) can be used to 
# determine if the structure of the microbial communities is significantly different between
# environmental variables. This is done using the adonis function from Vegan with a distance metric, e.g. Bray-Curtis, and a
# specified number of permutations, e.g. 10,000.
# The analysis measures the degree each environmental variable affects the community composition and indicates 
# the significance of that effect on beta diversity (described by p-values and R2 values). 
# The R2 value corresponds to the proportion of variability observed in the dissimilarity.


# The function below will only calculate the the significance of individual variables
run_permanova <- function(my_community_data, my_metadata, my_variables){
  stat_sig_table <- NULL
  
  # Remove NA entries from the metadata
  for (var_name in my_variables) {
    result <- adonis(my_community_data~get(var_name),data = my_metadata, permu=999,method="euclidean")
    SumOfSqs <- round(result$aov.tab$SumsOfSqs[1], 3)
    meanSqs <- round(result$aov.tab$MeanSqs[1], 3)
    F.model <- round(result$aov.tab$F.Model[1], 3)
    R2 <- round(result$aov.tab$R2[1], 3)
    p_value <- round(result$aov.tab$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(var_name,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       R2,
                                                       p_value))
  }
  names(stat_sig_table) <- c("Variable", "SumOfSqs","MeanSqs","F.Model","R2","P-value")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"P-value"),]
  stat_sig_table
}
# : = interaction. Only include the variable interaction, not the variables themselves.
# * = crossing. Include variable interaction and the variables.
# + = add the variable, independent
# %in% = nesting. e.g. b %in% a
# /    = nesting with main effect term. e.g. a + b %in% a.

adonis(t(otu_rare_clr_filtered.m)~Remote_Community,data = metadata.df, permu=999,method="euclidean")
adonis(t(otu_rare_clr_filtered.m)~Nose,data = metadata.df, permu=999,method="euclidean")

metadata_permanova.df <- subset(metadata.df, Remote_Community == 0)
otu_rare_clr_filtered_permanova.m <- otu_rare_clr_filtered.m[,rownames(metadata_permanova.df)]
adonis(t(otu_rare_clr_filtered_permanova.m)~Otitis_status,data = metadata_permanova.df, permu=999,method="euclidean")

metadata_permanova.df <- subset(metadata.df, Remote_Community == 1)
otu_rare_clr_filtered_permanova.m <- otu_rare_clr_filtered.m[,rownames(metadata_permanova.df)]
adonis(t(otu_rare_clr_filtered_permanova.m)~Otitis_status,data = metadata_permanova.df, permu=999,method="euclidean")




dim(metadata_permanova.df)
metadata_permanova.df <- metadata_permanova.df[complete.cases(metadata_permanova.df),]
dim(metadata_permanova.df)
otu_rare_clr_filtered_permanova.m <- otu_rare_clr_filtered.m[,rownames(metadata_permanova.df)]
adonis(t(otu_rare_clr_filtered_permanova.m)~Sample_Type + DX_Groups +  Gender + BMI + AGE,data = metadata_permanova.df, permu=999,method="euclidean")
adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups + Sample_Type +  Gender + BMI + AGE,data = metadata_permanova.df, permu=999,method="euclidean")

# Nested DX_group in Sample_type
# metadata_permanova.df <- metadata.df[,c("DX_Groups", "Sample_Type", "Gender","PPI_Medications", "BMI","AGE")]
metadata_permanova.df <- metadata.df[,c("DX_Groups", "Sample_Type")]
metadata_permanova.df <- metadata_permanova.df[complete.cases(metadata_permanova.df),]
otu_rare_clr_filtered_permanova.m <- otu_rare_clr_filtered.m[,colnames(otu_rare_clr_filtered.m) %in% rownames(metadata_permanova.df)]
adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups %in% Sample_Type,data = metadata_permanova.df, permu=999,method="euclidean")

# Nested DX_group in Sample_type with main effect term
adonis(t(otu_rare_clr_filtered_permanova.m)~Sample_Type/DX_Groups,data = metadata_permanova.df, permu=999,method="euclidean")

# Interaction DX_group and Sample_type 
adonis(t(otu_rare_clr_filtered_permanova.m)~Sample_Type:DX_Groups,data = metadata_permanova.df, permu=999,method="euclidean")
# adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups:Sample_Type,data = metadata_permanova.df, permu=999,method="euclidean")

# Just sample type and disease state
adonis(t(otu_rare_clr_filtered_permanova.m)~Sample_Type + DX_Groups,data = metadata_permanova.df, permu=999,method="euclidean")
# adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups:Sample_Type,data = metadata_permanova.df, permu=999,method="euclidean")


# 
metadata_permanova.df <- metadata.df[,c("DX_Groups", "Sample_Type", "Gender","PPI_Medications", "BMI","AGE")]
metadata_permanova.df <- metadata_permanova.df[complete.cases(metadata_permanova.df),]
otu_rare_clr_filtered_permanova.m <- otu_rare_clr_filtered.m[,colnames(otu_rare_clr_filtered.m) %in% rownames(metadata_permanova.df)]
adonis(t(otu_rare_clr_filtered_permanova.m)~Sample_Type + DX_Groups + PPI_Medications + BMI + AGE + Gender,data = metadata_permanova.df, permu=999,method="euclidean")

# adonis(t(otu_rare_clr_filtered.m)~Sampletype,data = metadata.df, permu=999,method="euclidean")
# adonis(t(otu_rare_clr_filtered.m)~Sampletype,data = metadata.df, permu=999,method="bray")


# adonis(t(otu_rare_clr_filtered.m)~Sampletype,data = metadata.df, permu=999,method="euclidean")
# adonis(t(otu_rare_clr_filtered.m)~Sampletype_pooled,data = metadata.df, permu=999,method="euclidean")
# adonis(t(otu_rare_clr_filtered.m)~Patient+Sampletype+Sampletype_pooled+Patient:Sampletype+Patient:Sampletype_pooled,data = metadata.df, permu=999,method="euclidean")
# adonis(t(otu_rare_clr_filtered.m)~Patient+Sampletype_pooled+Sampletype+Patient:Sampletype+Patient:Sampletype_pooled,data = metadata.df, permu=999,method="euclidean")

# filter the metadata and OTU table to only those samples that have non-NA entries across the variables of interest
# metadata_permanova.df <- metadata.df[variables]
metadata_permanova.df <- metadata.df[,c("DX_Groups", "Sample_Type", "Gender","PPI_Medications", "BMI","AGE")]
metadata_permanova.df <- metadata_permanova.df[complete.cases(metadata_permanova.df),]
otu_rare_clr_filtered_permanova.m <- otu_rare_clr_filtered.m[,colnames(otu_rare_clr_filtered.m) %in% rownames(metadata_permanova.df)]
# dim(otu_rare_clr_filtered_permanova.m)
# dim(metadata_permanova.df)

# PERMANOVA adjusted for, for example, gender, age, BMI, PPI usage and Sample type


adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups * Sample_Type ,data = metadata_permanova.df, permu=999,method="euclidean")
adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups + Sample_Type + BMI + AGE + Gender + PPI_Medications + DX_Groups:Sample_Type,data = metadata_permanova.df, permu=999,method="euclidean")
adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups + Sample_Type + BMI + AGE + Gender + PPI_Medications,data = metadata_permanova.df, permu=999,method="euclidean")

adonis(t(otu_rare_clr_filtered_permanova.m)~DX_Groups+Gender+Smoking+Dysphagia,data = metadata_permanova.df, permu=999,method="euclidean")
adonis(t(otu_rare_clr_filtered_permanova.m)~Smoking,data = metadata_permanova.df, permu=999,method="euclidean")

permanova_results_OTU <- run_permanova(t(otu_rare_clr_filtered_permanova.m), metadata_permanova.df, variables)
write.csv(permanova_results_OTU,file="Result_tables/stats_various/PERMANOVA_otu_clr_rarified.csv",row.names = F)

