#************************************
# Build heatmaps at OTU/ASV and varying taxonomy levels
#************************************



library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
# library(vegan)
library(reshape2)
# library(gplots)
library(pheatmap)
library(grid)



####################################
# Define various colour palletes
# Various colour palletes
my_colour_pallete <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
my_colour_pallete_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
my_colour_pallete_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
my_colour_pallete_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
my_colour_pallete_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
my_colour_pallete_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
my_colour_pallete_32_distinct <- c("#ea7e00","#ca0074","#d1c69b","#474007","#bb00ad","#9c80ff","#be3300","#542e72","#00b9f5","#09436b","#8b0036","#9ac8e6","#ff1059","#959eff","#154a11","#0290f4","#ff7762","#7dbf00","#ff8194","#834c00","#006e73","#f9bb5d","#d6c943","#017229","#00d3a8","#732427","#36e191","#6a8200","#efb3ea","#3227bb","#ff90e1","#e92a12")
# lesion_pallete_7 <- c("#8558d6","#6ee268","#d247ad","#c9d743","#d7453e","#59a237","#d78f2a")
# patient_pallete_45 <- c("#d64530","#585fb1","#795d97","#9e4773","#3f6921","#71692c","#a2b93c","#d571cc","#9b3e97","#33947a","#98ad66","#448a4e","#869ae0","#5ce7af","#e085a3","#dfdc87","#d19be2","#5cb735","#e38269","#3db6c0","#50b565","#50902c","#a98a2c","#dde84a","#db3d76","#5fe485","#7c8329","#b3e791","#6fe965","#5ebce9","#3c86c1","#2a6a45","#65b688","#6651d1","#af4ed3","#df872f","#56e4db","#737cea","#ac464b","#dd37b5","#995b2b","#daac6f","#92e2be","#a2e24b","#e0be3a")
patient_pallete_270 <- c("#456c00","#77bbff","#75fb3f","#273300","#f5006f","#ac008b","#125700","#ffef77","#00278e","#3d005e","#d84100","#015686","#01289f","#ff8c4c","#0070b5","#8015cd","#feffd6","#02d8aa","#019cf4","#4f2f00","#bbffe9","#c52900","#1b002c","#a3ac00","#5d9200","#f29fff","#231500","#934cfc","#988a00","#002cbb","#ffeb36","#ffa758","#f1f9ff","#000045","#b4004b","#602900","#390048","#e6c400","#00ce3c","#ff7bd0","#8cff56","#e60051","#b89aff","#00474b","#d5fbff","#ff79c2","#1d0016","#00635d","#ff8e33","#992300","#ff6e91","#ffa081","#534a00","#61002d","#ffe1c1","#8c0072","#00405d","#89ffc6","#607500","#64ff6f","#002e52","#9b97ff","#b1ccff","#02c5cd","#5dffba","#beff45","#00112b","#b8ff74","#7f0027","#0074cd","#005c6f","#3f00af","#dd7900","#cced00","#77ffd6","#ffc5b5","#99ffb1","#01ea72","#f0ff88","#007f63","#abff9d","#391200","#003a74","#114b00","#0a0100","#ff5fff","#ffccfb","#00d6b7","#c7ff93","#1efa5a","#005437","#f6af00","#a60024","#ffb7e6","#ea0043","#c7ffbc","#72ab00","#789300","#585500","#c3ff14","#00f398","#ab4a00","#9b7600","#85e5ff","#006235","#130020","#006825","#ff735c","#007a7f","#02a3a4","#4856ff","#bf52ff","#00edbc","#a31bd6","#009642","#e93bee","#e400ae","#ffbdd2","#00cfc7","#f1ffaa","#009b7a","#dd00c9","#ff697d","#004a14","#ff72ac","#ff3d1f","#fffaa3","#5d0400","#027ba4","#01c774","#002655","#00941f","#0a32d7","#82acff","#ff8af3","#ff4165","#001104","#ffd6f2","#efebff","#aebc00","#3e0030","#c5abff","#00402e","#ff4bae","#0275f1","#be89ff","#ffd899","#00c765","#01b208","#97ffd4","#7e9fff","#00fde1","#0050c9","#ff8eb5","#c800cd","#005173","#ff2b95","#76ff7a","#ea0074","#001d70","#009856","#f100a8","#ba6b00","#0293df","#00462d","#ff6862","#f6ff65","#02bbda","#2c2200","#01a876","#e35a00","#e3000f","#ff819e","#5a0039","#a558ff","#e2ffb2","#784800","#016def","#b400a2","#00143c","#00212d","#403d00","#ff75fe","#975300","#166c00","#260008","#917fff","#ff8d89","#01bf7a","#ffa6bf","#800086","#90a100","#cce4ff","#dad800","#52c900","#46a700","#0c0039","#0b0052","#79009d","#003c85","#bb0034","#59e7ff","#af0064","#64001e","#c0007e","#000897","#bd8400","#2b007f","#318400","#31f1ff","#7c8600","#807300","#ffc072","#6f005f","#770040","#e62c00","#2e0019","#005599","#6535e1","#5b0099","#006bd5","#0142a1","#baaf00","#00ab2d","#ffcc40","#edffec","#ef0031","#153100","#abe9ff","#6bbd00","#e5ff4e","#ffdb43","#ffa5ef","#01c4f3","#ffbd8f","#84d100","#bbff84","#9fcdff","#7b0059","#ffe897","#ff8711","#ffa869","#febaff","#20003a","#94002b","#5387ff","#756dff","#fff494","#a5c1ff","#e0ffcf","#002417","#530076","#ff8459","#ffe4ec","#00b650","#0119b7","#c963ff","#a2ff64","#9c6800","#03b6f8","#00a0c2","#00240b","#6297ff","#bd0010","#fff7af","#7d2d00","#cf7aff","#af5600","#322c00","#500028")
my_colour_pallete_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
my_colour_pallete_10_soft <- c("#9E788F","#4C5B61","#678D58","#AD5233","#A0A083","#4D456A","#588578","#D0AC4C","#2A7BA0","#931621")
####################################
# Function
log_matrix <- function(mymat){
  out <- log(mymat, 10)
  out[is.infinite(out)] <- 0
  return(out)
}


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


filter_heatmap_matrix <- function(myheatmap, row_max = 0, prevalence = 0){
  internal_heatmap <- myheatmap
  internal_heatmap <- internal_heatmap[which(apply(internal_heatmap, 1, max) >= row_max), ]
  # keep only OTUs/taxa that are in more than this fraction of samples
  filter_fraction <- prevalence
  entry_prevalences <- apply(internal_heatmap, 1, function(x) {length(which(x > 0))})/dim(internal_heatmap)[2]
  entries_from_prevalences <- names(entry_prevalences)[entry_prevalences > filter_fraction]
  entries_from_prevalences <- entries_from_prevalences[!is.na(entries_from_prevalences)]
  return(internal_heatmap[entries_from_prevalences,])
}

####################################

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project/")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)


# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index


# Load count matrices and covert to log space
otu_rare_log.m <- log_matrix(as.matrix(read.table(file = "Result_tables/count_tables/OTU_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_genus_rare_log.m <- log_matrix(as.matrix(read.table(file = "Result_tables/count_tables/Genus_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_family_rare_log.m <- log_matrix(as.matrix(read.table(file = "Result_tables/count_tables/Family_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_order_rare_log.m <- log_matrix(as.matrix(read.table(file = "Result_tables/count_tables/Order_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_class_rare_log.m <- log_matrix(as.matrix(read.table(file = "Result_tables/count_tables/Class_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_phylum_rare_log.m <- log_matrix(as.matrix(read.table(file = "Result_tables/count_tables/Phylum_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))

# Load count matrices and apply CLR transformation
otu_rare_clr.m <- clr(as.matrix(read.table(file = "Result_tables/count_tables/OTU_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_rare_species_clr.m <- clr(as.matrix(read.table(file = "Result_tables/count_tables/Specie_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_rare_genus_clr.m <- clr(as.matrix(read.table(file = "Result_tables/count_tables/Genus_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_rare_family_clr.m <- clr(as.matrix(read.table(file = "Result_tables/count_tables/Family_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_rare_order_clr.m <- clr(as.matrix(read.table(file = "Result_tables/count_tables/Order_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_rare_class_clr.m <- clr(as.matrix(read.table(file = "Result_tables/count_tables/Class_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))
otu_rare_phylum_clr.m <- clr(as.matrix(read.table(file = "Result_tables/count_tables/Phylum_counts_rarefied.csv", sep = ",", header = T, row.names = 1)))

# Load relative abundance matrices
otu_genus_rare_rel.m <- as.matrix(read.table(file = "Result_tables/relative_abundance_tables/Genus_relative_abundances_rarefied.csv", sep = ",", header = T, row.names = 1))

# Cleanup column (sample) names - log matrices
colnames(otu_rare_log.m) <- gsub("_J.*", "", colnames(otu_rare_log.m))
colnames(otu_genus_rare_log.m) <- gsub("_J.*", "", colnames(otu_genus_rare_log.m))
colnames(otu_family_rare_log.m) <- gsub("_J.*", "", colnames(otu_family_rare_log.m))
colnames(otu_order_rare_log.m) <- gsub("_J.*", "", colnames(otu_order_rare_log.m))
colnames(otu_class_rare_log.m) <- gsub("_J.*", "", colnames(otu_class_rare_log.m))
colnames(otu_phylum_rare_log.m) <- gsub("_J.*", "", colnames(otu_phylum_rare_log.m))

# Cleanup column (sample) names - CLR matrices
colnames(otu_rare_clr.m) <- gsub("_J.*", "", colnames(otu_rare_clr.m))
colnames(otu_rare_species_clr.m) <- gsub("_J.*", "", colnames(otu_rare_species_clr.m))
colnames(otu_rare_genus_clr.m) <- gsub("_J.*", "", colnames(otu_rare_genus_clr.m))
colnames(otu_rare_family_clr.m) <- gsub("_J.*", "", colnames(otu_rare_family_clr.m))
colnames(otu_rare_order_clr.m) <- gsub("_J.*", "", colnames(otu_rare_order_clr.m))
colnames(otu_rare_class_clr.m) <- gsub("_J.*", "", colnames(otu_rare_class_clr.m))
colnames(otu_rare_phylum_clr.m) <- gsub("_J.*", "", colnames(otu_rare_phylum_clr.m))

# Cleanup column (sample) names - Relative abundance matrices
colnames(otu_genus_rare_rel.m) <- gsub("_J.*", "", colnames(otu_genus_rare_rel.m))

# and correct metadata
rownames(metadata.df) <- gsub("_J.*", "", rownames(metadata.df))
metadata.df$Index <- gsub("_J.*", "", metadata.df$Index)


# Since we likely removed samples from the count matrix
# in the main script, remove them from the metadata.df here
samples_removed <- metadata.df$Index[!metadata.df$Index %in% colnames(otu_genus_rare_log.m)]
metadata.df <- metadata.df[! metadata.df$Index %in% samples_removed,]

# Remove samples that are not in the metadata.
otu_rare_log.m <- otu_rare_log.m[,colnames(otu_rare_log.m) %in% metadata.df$Index]
otu_genus_rare_log.m <- otu_genus_rare_log.m[,colnames(otu_genus_rare_log.m) %in% metadata.df$Index]
otu_family_rare_log.m <- otu_family_rare_log.m[,colnames(otu_family_rare_log.m) %in% metadata.df$Index]
otu_order_rare_log.m <- otu_order_rare_log.m[,colnames(otu_order_rare_log.m) %in% metadata.df$Index]
otu_class_rare_log.m <- otu_class_rare_log.m[,colnames(otu_class_rare_log.m) %in% metadata.df$Index]
otu_phylum_rare_log.m <- otu_phylum_rare_log.m[,colnames(otu_phylum_rare_log.m) %in% metadata.df$Index]

# Remove samples that are not in the metadata.
otu_rare_clr.m <- otu_rare_clr.m[,colnames(otu_rare_clr.m) %in% metadata.df$Index]
otu_rare_species_clr.m <- otu_rare_species_clr.m[,colnames(otu_rare_species_clr.m) %in% metadata.df$Index]
otu_rare_genus_clr.m <- otu_rare_genus_clr.m[,colnames(otu_rare_genus_clr.m) %in% metadata.df$Index]
otu_rare_family_clr.m <- otu_rare_family_clr.m[,colnames(otu_rare_family_clr.m) %in% metadata.df$Index]
otu_rare_order_clr.m <- otu_rare_order_clr.m[,colnames(otu_rare_order_clr.m) %in% metadata.df$Index]
otu_rare_class_clr.m <- otu_rare_class_clr.m[,colnames(otu_rare_class_clr.m) %in% metadata.df$Index]
otu_rare_phylum_clr.m <- otu_rare_phylum_clr.m[,colnames(otu_rare_phylum_clr.m) %in% metadata.df$Index]

# Remove samples that are not in the metadata.
otu_genus_rare_rel.m <- otu_genus_rare_rel.m[,colnames(otu_genus_rare_rel.m) %in% metadata.df$Index,drop=F]


# If there are negative values in the CLR matrices, assign them a value of zero
otu_rare_clr.m[which(otu_rare_clr.m < 0)] <- 0
otu_rare_species_clr.m[which(otu_rare_species_clr.m < 0)] <- 0
otu_rare_genus_clr.m[which(otu_rare_genus_clr.m < 0)] <- 0
otu_rare_family_clr.m[which(otu_rare_family_clr.m < 0)] <- 0
otu_rare_order_clr.m[which(otu_rare_order_clr.m < 0)] <- 0
otu_rare_class_clr.m[which(otu_rare_class_clr.m < 0)] <- 0
otu_rare_phylum_clr.m[which(otu_rare_phylum_clr.m < 0)] <- 0


# Function to create heatmap
make_heatmap <- function(myheatmap_matrix,
                         mymetadata,
                         filename,
                         height,
                         width,
                         variables,
                         clust_cols = F,
                         clust_rows = F,
                         treeheight_col = 0,
                         treeheight_row = 0,
                         cutree_cols = 1, 
                         annotation_pallete,
                         my_breaks = NULL,
                         variable_colours_available = F,
                         ...){
  
  internal_heatmap_matrix.m <- myheatmap_matrix

  sample_annotations.df <- data.frame(matrix(nrow = length(colnames(internal_heatmap_matrix.m)), ncol = length(variables)))
  rownames(sample_annotations.df) <- colnames(internal_heatmap_matrix.m)
  colnames(sample_annotations.df) <- variables
  
  # Determine the values for each variable for each sample
  for (sample_id in rownames(sample_annotations.df)){
    for (var in variables){
      sample_variable_value <- paste(mymetadata[mymetadata$Index == sample_id,var])
      sample_annotations.df[sample_id, var] <- sample_variable_value
    }
  }
  for (var in variables){
    # Assign a distinct colour to each unique value for each variable
    matched_colours <- NULL
    i = 1
    for (var_val in unique(sample_annotations.df[,var])){
      if (is.na(var_val) | var_val == "" | var_val == "NA"){
        matched_colours <- rbind(matched_colours, data.frame(Group = var_val, Colour = "grey"))  
      } else{
        color_col_name <- paste0(var, "_colour")
        if (variable_colours_available == T & color_col_name %in% names(mymetadata)){
          metadata_subset.df <- subset(mymetadata, get(var) == var_val)
          meta_colour <- unique(as.character(metadata_subset.df[,color_col_name]))
          matched_colours <- rbind(matched_colours, data.frame(Group = var_val, Colour = meta_colour))
        }else{
          matched_colours <- rbind(matched_colours, data.frame(Group = var_val, Colour = annotation_pallete[i]))    
        }
        
      }
      
      i = i + 1
    }
    sample_annotations.df[,paste0(var, "_colour")] <- unlist(lapply(sample_annotations.df[,var], function(x) as.character(matched_colours[matched_colours$Group == x,]$Colour)))
  }
  
  # ********************************************************************
  
  # Colour pallete for heatmap
  if (is.null(my_breaks)){
    my_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))(13)  
  } else{
    my_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))(length(my_breaks))  
  }
  
  
  
  heatmap_ordered.m <- internal_heatmap_matrix.m[,do.call(order, sample_annotations.df[variables])]
  sample_annotations.df <- sample_annotations.df[do.call(order, sample_annotations.df[variables]),]
  
  # Create an annotation list for each unique group-colour
  annotation_colours <- list()
  for (var in variables){
    annotation_colours_var <- unique(sample_annotations.df[,paste0(var, "_colour")])
    names(annotation_colours_var) <- unique(sample_annotations.df[,var])
    annotation_colours[[var]] <- annotation_colours_var
  }
  
 
  
  pdf(filename,height=height,width=width)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(heatmap_ordered.m,
           treeheight_row = treeheight_row, # To hide row tree
           treeheight_col = treeheight_col, # To hide column tree
           cutree_cols = cutree_cols,
           # treeheight_col = 100,
           # gaps_col = annotation_breaks,
           # breaks = my_breaks,
           # scale = "column",
           annotation_names_col = T,
           color = my_palette,
           fontsize_row = 6,
           fontsize_col = 6,
           cluster_rows = clust_rows,
           cluster_cols = clust_cols,
           # legend_breaks = my_breaks,
           # legend_labels = my_breaks,
           annotation_col  = sample_annotations.df[,variables,  drop = FALSE],
           annotation_colors = annotation_colours,
           ...)
  setHook("grid.newpage", NULL, "replace")
  grid.text("Sample", y=-0.02, x =0.43, gp=gpar(fontsize=16))
  grid.text("Taxa", x=-0.02, y = 0.5, rot=90, gp=gpar(fontsize=16),)
  dev.off()
}

# Define the discrete variables
discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season","Nose",
                        "Otitis_status_OM_6mo","Remote_Community_Otitis_status","OM_6mo_Type_OM","Remote_Community_Season")

# discrete_variables <- c("Gold_Star","OM_6mo")
# discrete_variables <- c("Nose")
# Convert variables to factors
metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)

# log(50,10) = 1.69897, 0.1% assuming 50,000
# log(100,10) = 2, 0.2% assuming 50,000
# log(250,10) = 2.39794, 0.5% assuming 50,000
# log(500,10) = 2.69897, 1% assuming 50,000
# 10**2.69897 = 500

heatmap_genus.m <- filter_heatmap_matrix(otu_genus_rare_log.m, row_max = 2.39794, prevalence = 0.0)
heatmap_genus_clr.m <- filter_heatmap_matrix(otu_rare_genus_clr.m, row_max = 0, prevalence = 0.1)
heatmap_genus_rel.m <- filter_heatmap_matrix(otu_genus_rare_rel.m, row_max = 0.005, prevalence = 0.1)

apply(heatmap_genus_rel.m*100, 1, min)

seq(0,.01, .005)
seq(.01,1,.05)

my_breaks = c(c(0,0.5,1,5),seq(10,100,10))
my_breaks

make_heatmap(myheatmap_matrix = heatmap_genus_rel.m*100,
             mymetadata = metadata.df,
             height = 5,
             width = 20,
             filename = paste0("Result_figures/heatmaps/Remote_Community_genus_relative_abundance.pdf"), 
             variables = "Remote_Community",
             clustering_method = "average",
             # clustering_method = "complete",
             # my_breaks = my_breaks,
             breaks = my_breaks,
             legend_breaks = my_breaks,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 0, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct,variable_colours_available = T)
dev.off()


my_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))(10)
pheatmap(
  heatmap_genus_rel.m,
  clustering_method = "average",
  # clustering_method = "average",
  # clustering_distance = "euclidean",
  fontsize = 6,
  color = my_palette,
  breaks = seq(0,1,.1),
  legend_breaks = seq(0,1,.1),
  # legend_labels = 
)


# otu_rare_genus_clr.m
# make_heatmap(myheatmap_matrix = heatmap_genus.m,
#              mymetadata = metadata.df,
#              height = 5,
#              width = 20,
#              filename = "Result_figures/heatmaps/genus_log_read_count.pdf", 
#              variables = discrete_variables,
#              clust_cols = T,clust_rows = T, treeheight_col = 0, treeheight_row = 0, cutree_cols = 1,
#              annotation_pallete = my_colour_pallete_10_distinct,variable_colours_available = T)

for (var in discrete_variables){
  # make_heatmap(myheatmap_matrix = heatmap_genus.m,
  #              mymetadata = metadata.df,
  #              height = 5,
  #              width = 20,
  #              filename = paste0("Result_figures/heatmaps/",var,"_genus_log_read_count.pdf"), 
  #              variables = var,
  #              clust_cols = F,clust_rows = T, treeheight_col = 0, treeheight_row = 0, cutree_cols = 1,
  #              annotation_pallete = my_colour_pallete_10_distinct,variable_colours_available = T)
  # Genus
  # make_heatmap(myheatmap_matrix = heatmap_genus_clr.m,
  #              mymetadata = metadata.df,
  #              height = 5,
  #              width = 20,
  #              filename = paste0("Result_figures/heatmaps/",var,"_genus_read_count_CLR.pdf"), 
  #              variables = var,
  #              clust_cols = F,clust_rows = T, treeheight_col = 0, treeheight_row = 0, cutree_cols = 1,
  #              annotation_pallete = my_colour_pallete_10_distinct, variable_colours_available =T)
  make_heatmap(myheatmap_matrix = heatmap_genus_rel.m,
               mymetadata = metadata.df,
               height = 5,
               width = 20,
               filename = paste0("Result_figures/heatmaps/",var,"_genus_relative_abundance.pdf"), 
               variables = var,
               clust_cols = F,clust_rows = T, treeheight_col = 0, treeheight_row = 0, cutree_cols = 1,
               annotation_pallete = my_colour_pallete_10_distinct,variable_colours_available = T)
}



for (community in unique(metadata.df$Remote_Community)){
  for (var in discrete_variables){
    if (var == "Remote_Community") {next}
    make_heatmap(myheatmap_matrix = heatmap_genus.m,
                 mymetadata = metadata.df,
                 height = 5,
                 width = 20,
                 filename = paste0("Result_figures/heatmaps/Remote_Community_and_",var,"_genus_log_read_count.pdf"), 
                 variables = c(var, "Remote_Community"),
                 clust_cols = F,clust_rows = T, treeheight_col = 0, treeheight_row = 0, cutree_cols = 1,
                 annotation_pallete = my_colour_pallete_10_distinct,variable_colours_available = T)
  }
}


#### Filter the matrices to only the most abundant taxa and plot them
##  Log matrices
# Phylum
heatmap_phylum.m <- filter_heatmap_matrix(otu_phylum_rare_log.m, row_max = 1, prevalence = 0.1)
make_heatmap(myheatmap_matrix = heatmap_phylum.m,
             height = 5,
             width = 20,
             filename = "Result_figures/heatmaps/phylum_log_read_count.pdf", 
             variables = discrete_variables,
             clust_cols = F,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_soft,variable_colours_available = T)


# Class
heatmap_class.m <- filter_heatmap_matrix(otu_class_rare_log.m, row_max = 1, prevalence = 0.1)
make_heatmap(myheatmap_matrix = heatmap_class.m,
             height = 5,
             width = 95,
             filename = "Result_figures/heatmaps/class_log_read_count.pdf", 
             variables = discrete_variables,
             clust_cols = F,clust_rows = T, treeheight_col = 0, treeheight_row = 0, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_soft)


# Order
heatmap_order.m <- filter_heatmap_matrix(otu_order_rare_log.m, row_max =1, prevalence = 0.1)
make_heatmap(myheatmap_matrix = heatmap_order.m,
             height = 8,
             width = 95,
             filename = "Result_figures/heatmaps/order_log_read_count.pdf", 
             variables = discrete_variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)


# Family
heatmap_family.m <- filter_heatmap_matrix(otu_family_rare_log.m, row_max = 1, prevalence = 0.2)
make_heatmap(myheatmap_matrix = heatmap_family.m,
             height = 12,
             width = 75,
             filename = "Result_figures/heatmaps/family_log_read_count.pdf", 
             variables = discrete_variables,
             clust_cols = F,clust_rows = T, treeheight_col = 0, treeheight_row = 0, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_soft)


# Genus
heatmap_genus.m <- filter_heatmap_matrix(otu_genus_rare_log.m, row_max = 2, prevalence = 0.2)
make_heatmap(myheatmap_matrix = heatmap_genus.m,
             height = 20,
             width = 95,
             filename = "Result_figures/heatmaps/genus_log_read_count.pdf", 
             variables = discrete_variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)



# OTU
heatmap_otu.m <- filter_heatmap_matrix(otu_rare_log.m, row_max = 1, prevalence = 0.1)
make_heatmap(myheatmap_matrix = heatmap_otu.m,
             height = 25,
             width = 95,
             filename = "Result_figures/heatmaps/otu_log_read_count.pdf", 
             variables = discrete_variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)


#### CLR matrices

# Phylum
# metadata.df$Remote_Community_colour <- NULL
heatmap_phylum_clr.m <- filter_heatmap_matrix(otu_rare_phylum_clr.m, row_max = 0, prevalence = 0.1)
make_heatmap(myheatmap_matrix = heatmap_phylum_clr.m,
             height = 5,
             width = 95,
             filename = "Result_figures/heatmaps/phylum_clr.pdf", 
             variables = variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)

# Class
heatmap_class_clr.m <- filter_heatmap_matrix(otu_rare_class_clr.m, row_max = 0, prevalence = 0.1)
make_heatmap(myheatmap_matrix = heatmap_class_clr.m,
             height = 5,
             width = 95,
             filename = "Result_figures/heatmaps/class_clr.pdf", 
             variables = variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)

# Order
heatmap_order_clr.m <- filter_heatmap_matrix(otu_rare_order_clr.m, row_max = 0, 0.05)
make_heatmap(myheatmap_matrix = heatmap_order_clr.m,
             height = 8,
             width = 95,
             filename = "Result_figures/heatmaps/order_clr.pdf", 
             variables = variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)

# Family
heatmap_family_clr.m <- filter_heatmap_matrix(otu_rare_family_clr.m, row_max = 0, 0.1)
make_heatmap(myheatmap_matrix = heatmap_family_clr.m,
             height = 12,
             width = 85,
             filename = "Result_figures/heatmaps/family_clr.pdf", 
             variables = variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)

# Genus
heatmap_genus_clr.m <- filter_heatmap_matrix(otu_rare_genus_clr.m, row_max = 0, 0.2)
make_heatmap(myheatmap_matrix = heatmap_genus_clr.m,
             height = 20,
             width = 95,
             filename = "Result_figures/heatmaps/genus_clr.pdf", 
             variables = variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)

# OTU
heatmap_otu_clr.m <- filter_heatmap_matrix(otu_rare_clr.m, row_max = 0, 0.1)
make_heatmap(myheatmap_matrix = heatmap_otu_clr.m,
             height = 25,
             width = 95,
             filename = "Result_figures/heatmaps/otu_clr.pdf", 
             variables = variables,
             clust_cols = T,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)


# ------------------------------------------------------------
# Load significant OTU / Genus according to DESeq and make a unfiltered heatmap for them 

# TODO want to be able to provide custom axis labels
# Could provide a dataframe that has a specified mapping between the column (or row) name and the label.
# Column_name Column_label Row_name Row_label
# Currently just changing the rownames, this approach won't work with columns

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

wheat_otu_deseq.df <- read.csv("Result_tables/DESeq_results/OTU/Wheat_related_symptoms_and_group_combined.csv")
otus_of_interest.m <- otu_rare_log.m[unique(wheat_otu_deseq.df$OTU),]
rownames(otus_of_interest.m) <- as.character(lapply(rownames(otus_of_interest.m), function(x) paste0(x,"; ", otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$Genus)))
make_heatmap(myheatmap_matrix = otus_of_interest.m,
             height = 15,
             width = 95,
             filename = "Result_figures/heatmaps/DESeq/Wheat_related_symptoms_and_group_otu_log_read_count.pdf", 
             variables = "Wheat_related_symptoms_and_group",
             clust_cols = F,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)


DX_Groups_deseq.df <- read.csv("Result_tables/DESeq_results/OTU/DX_Groups_combined.csv")
otus_of_interest.m <- otu_rare_log.m[unique(DX_Groups_deseq.df$OTU),]
rownames(otus_of_interest.m) <- as.character(lapply(rownames(otus_of_interest.m), function(x) paste0(x,"; ", otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$Genus)))
make_heatmap(myheatmap_matrix = otus_of_interest.m,
             height = 25,
             width = 95,
             filename = "Result_figures/heatmaps/DESeq/DX_Groups_otu_log_read_count.pdf", 
             variables = "DX_Groups",
             clust_cols = F,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)


Sample_type_deseq.df <- read.csv("Result_tables/DESeq_results/OTU/Sample_type_combined.csv")
otus_of_interest.m <- otu_rare_log.m[unique(Sample_type_deseq.df$OTU),]
rownames(otus_of_interest.m) <- as.character(lapply(rownames(otus_of_interest.m), function(x) paste0(x,"; ", otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$Genus)))
make_heatmap(myheatmap_matrix = otus_of_interest.m,
             height = 25,
             width = 95,
             filename = "Result_figures/heatmaps/DESeq/Sample_type_otu_log_read_count.pdf", 
             variables = "Sample_Type",
             clust_cols = F,clust_rows = T, treeheight_col = 60, treeheight_row = 60, cutree_cols = 1,
             annotation_pallete = my_colour_pallete_10_distinct)

