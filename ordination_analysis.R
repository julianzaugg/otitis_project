detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

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

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index
rownames(metadata_decontaminated.df) <- metadata_decontaminated.df$Index

# Convert variables to factors
# discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season",
# "Nose","Otitis_status_OM_6mo", "Remote_Community_Otitis_status", "OM_6mo_Type_OM","Remote_Community_Season")
discrete_variables <- c("Remote_Community","Gold_Star","OM_6mo","Season","Nose","OM_Classification", "Remote_Community_Season",
                        "Streptococcus_pneumoniae", "Moraxella_catarrhalis", "Haemophilus_influenzae",
                        "Remote_Community_OM_Classification")

metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)
metadata_decontaminated.df[discrete_variables] <- lapply(metadata_decontaminated.df[discrete_variables], factor)

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the counts
otu.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts.csv", header =T, row.names = 1))
genus.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts.csv", header =T, row.names = 1))

otu_decontaminated.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts_decontaminated.csv", header =T, row.names = 1))
genus_decontaminated.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts_decontaminated.csv", header =T, row.names = 1))

colSums(otu.m)
# Remove samples from the OTU table that are not in the filtered metadata
otu.m <- otu.m[,colnames(otu.m) %in% c("OTU.ID", as.character(metadata.df$Index))]
genus.m <- genus.m[,colnames(genus.m) %in% c("OTU.ID", as.character(metadata.df$Index))]

otu_decontaminated.m <- otu_decontaminated.m[,colnames(otu_decontaminated.m) %in% c("OTU.ID", as.character(metadata_decontaminated.df$Index))]
genus_decontaminated.m <- genus_decontaminated.m[,colnames(genus_decontaminated.m) %in% c("OTU.ID", as.character(metadata_decontaminated.df$Index))]


# Remove samples from metadata that are not in the data
# metadata.df <- metadata.df[metadata.df$Index %in% colnames(otu.df),]
# metadata_decontaminated.df <- metadata_decontaminated.df[metadata_decontaminated.df$Index %in% colnames(otu_decontaminated.df),]

# Order the metadata.df by the index value
metadata.df <- metadata.df[order(metadata.df$Index),]
metadata_decontaminated.df <- metadata_decontaminated.df[order(metadata_decontaminated.df$Index),]

# Order the matrices and metadata to be the same order
metadata.df <- metadata.df[order(rownames(metadata.df)),]
metadata_decontaminated.df <- metadata_decontaminated.df[order(rownames(metadata_decontaminated.df)),]

otu.m <- otu.m[,order(rownames(metadata.df))]
genus.m <- genus.m[,order(rownames(metadata.df))]
otu_decontaminated.m <- otu_decontaminated.m[,order(rownames(metadata_decontaminated.df))]
genus_decontaminated.m <- genus_decontaminated.m[,order(rownames(metadata_decontaminated.df))]

# CLR transform the otu matrix.
otu_clr.m <- clr(otu.m)
genus_clr.m <- clr(genus.m)

otu_clr_decontaminated.m <- clr(otu_decontaminated.m)
genus_clr_decontaminated.m <- clr(genus_decontaminated.m)


# If there are negative values, assign them a value of zero
# otu_rare_clr_filtered.m[which(otu_rare_clr_filtered.m < 0)] <- 0
# otu_clr_filtered.m[which(otu_clr_filtered.m < 0)] <- 0


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# Ordination analysis

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
# otu_relabeller_function(rownames(otu.m))

genus_relabeller_function <- function(my_labels){
  unlist(lapply(my_labels, 
                function(x) {
                  phylostring <- unlist(strsplit(x, split = ";"))
                  # paste(phylostring[2],phylostring[3], phylostring[6], sep = ";")
                  paste(phylostring[3], phylostring[6], sep = ";")
                }))
}

# --------------------
# Generate PCA plots. If CLR transformed values with euclidean distances, these will be the same as
# the values calculated from betadisper...maybe not important
temp <- betadiver(t(otu_clr.m),method = "e")
# temp <- with(metadata.df, betadisper...maybe(vegdist(t(otu_clr.m),method = "euclidean"),group = Gold_Star))
# temp$eig["PCoA1"] / sum(temp$eig)
# temp$eig["PCoA2"] / sum(temp$eig)
# plot(temp)

# Generate ordination objects
otu_pca <- rda(t(otu_clr.m), data = metadata.df)
genus_pca <- rda(t(genus_clr.m), data = metadata.df)
otu_decontaminated_pca <- rda(t(otu_clr_decontaminated.m), data = metadata_decontaminated.df)
genus_decontaminated_pca <- rda(t(genus_clr_decontaminated.m), data = metadata_decontaminated.df)


source("code/helper_functions.R")
# calculate_PC_taxa_contributions(genus_pca)
otu_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata.csv",header = T)
otu_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata_decontaminated.csv",header = T)

genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv",header = T)
genus_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata_decontaminated.csv",header = T)

temp <- calculate_PC_abundance_correlations(genus_pca, mydata.df = genus_data.df,taxa_column = "taxonomy_genus",variables = discrete_variables)
temp <- calculate_PC_abundance_correlations(genus_decontaminated_pca, mydata.df = genus_data_decontaminated.df,taxa_column = "taxonomy_genus",variables = discrete_variables)
temp <- calculate_PC_abundance_correlations(otu_decontaminated_pca, mydata.df = otu_data_decontaminated.df,taxa_column = "OTU.ID",variables = discrete_variables)
# temp %>% filter(N_Samples > 5)
# temp <- rda(t(clr(t(rrarefy(t(genus_decontaminated.m[,colSums(genus_decontaminated.m) >= 3000]), 3000)))))
# temp_meta <- metadata_decontaminated.df[rownames(temp$CA$u),]
# par(mfrow=c(1, 2))
# generate_pca(temp,  mymetadata = temp_meta,variable_to_plot = "Nose", colour_palette = my_colour_palette_15)
# generate_pca(genus_decontaminated_pca,  mymetadata = metadata_decontaminated.df,variable_to_plot = "Nose", colour_palette = my_colour_palette_15)
generate_pca(otu_pca, mymetadata = metadata.df,
             plot_height = 5, plot_width = 5,
             legend_x = -8, legend_y = 4,
             # legend_x = -2, legend_y = 2,
             point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
             legend_title = myvar,
             legend_cex = .5,
             plot_title = myvar,
             limits = c(-8,5,-6,5),
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
             plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
             label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
             specie_labeller_function = otu_relabeller_function,arrow_label_offset = 0)


for (myvar in discrete_variables){
  
  # OTU, normal
  generate_pca(otu_pca, mymetadata = metadata.df,
               plot_height = 5, plot_width = 5,
               legend_x = -8, legend_y = 4,
               # legend_x = -2, legend_y = 2,
               point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
               legend_title = myvar,
               legend_cex = .5,
               plot_title = myvar,
               limits = c(-8,5,-6,5),
               plot_spiders = F,
               plot_ellipses = F,
               plot_hulls = F,
               use_shapes = T,
               ellipse_border_width = .5,
               include_legend = T,
               label_ellipse = F, ellipse_label_size = .3,
               colour_palette = my_colour_palette_15,
               variable_to_plot = myvar, legend_cols = 1,
               variable_colours_available = T,
               num_top_species = 3,
               plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
               label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
               specie_labeller_function = otu_relabeller_function,arrow_label_offset = 0,
               filename = paste0("Result_figures/pca_plots/otu/otu_",myvar, "_pca.pdf"))

  # Genus, normal
  generate_pca(genus_pca, mymetadata = metadata.df,
               plot_height = 5, plot_width = 5,
               legend_x = -4, legend_y = 3,
               # legend_x = -2, legend_y = 2,
               point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
               legend_title = myvar,
               legend_cex = .5,
               plot_title = myvar,
               limits = c(-4,6,-4,4),
               plot_spiders = F,
               plot_ellipses = F,
               plot_hulls = F,
               use_shapes = T,
               ellipse_border_width = .5,
               include_legend = T,
               label_ellipse = F, ellipse_label_size = .3,
               colour_palette = my_colour_palette_15,
               variable_to_plot = myvar, legend_cols = 1,
               variable_colours_available = T,
               num_top_species = 3,
               plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
               label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
               specie_labeller_function = genus_relabeller_function,
               arrow_label_offset = 0,
               filename = paste0("Result_figures/pca_plots/genus/genus_",myvar, "_pca.pdf"))
  # OTU, decontaminated
  generate_pca(otu_decontaminated_pca, mymetadata = metadata_decontaminated.df,
               plot_height = 5, plot_width = 5,
               legend_x = -8, legend_y = 4,
               # legend_x = -2, legend_y = 2,
               point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
               legend_title = myvar,
               legend_cex = .5,
               plot_title = myvar,
               limits = c(-8,4,-8,4),
               plot_spiders = F,
               plot_ellipses = F,
               plot_hulls = F,
               use_shapes = T,
               ellipse_border_width = .5,
               include_legend = T,
               label_ellipse = F, ellipse_label_size = .3,
               colour_palette = my_colour_palette_15,
               variable_to_plot = myvar, legend_cols = 1,
               variable_colours_available = T,
               num_top_species = 3,
               plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
               label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
               specie_labeller_function = otu_relabeller_function,
               arrow_label_offset = 0,
               filename = paste0("Result_figures/pca_plots/otu_decontaminated/otu_decontaminated_",myvar, "_pca.pdf"))
  # Genus, decontaminated
  generate_pca(genus_decontaminated_pca, mymetadata = metadata_decontaminated.df,
               plot_height = 5, plot_width = 5,
               legend_x = -4, legend_y = 4,
               # legend_x = -2, legend_y = 2,
               point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
               legend_title = myvar,
               legend_cex = .5,
               plot_title = myvar,
               limits = c(-4,6,-4,4),
               plot_spiders = F,
               plot_ellipses = F,
               plot_hulls = F,
               use_shapes = T,
               ellipse_border_width = .5,
               include_legend = T,
               label_ellipse = F, ellipse_label_size = .3,
               colour_palette = my_colour_palette_15,
               variable_to_plot = myvar, legend_cols = 1,
               variable_colours_available = T,
               num_top_species = 3,
               plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
               label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
               specie_labeller_function = genus_relabeller_function,
               arrow_label_offset = 0,
               filename = paste0("Result_figures/pca_plots/genus_decontaminated/genus_decontaminated_",myvar, "_pca.pdf"))
}

# TODO Perform PCA on only samples within the same community and generate plots
source("code/helper_functions.R")
for (community in unique(metadata.df$Remote_Community)){
  # Generate community specific metadata
  metadata_subset.df <- subset(metadata.df, Remote_Community == community)
  metadata_decontaminated_subset.df <- subset(metadata_decontaminated.df, Remote_Community == community)
  
  # Generate ordination objects for the community data
  otu_pca <- rda(t(otu_clr.m[,rownames(metadata_subset.df)]), data = metadata_subset.df)
  genus_pca <- rda(t(genus_clr.m[,rownames(metadata_subset.df)]), data = metadata_subset.df)
  otu_decontaminated_pca <- rda(t(otu_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]), data = metadata_decontaminated_subset.df)
  genus_decontaminated_pca <- rda(t(genus_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]), data = metadata_decontaminated_subset.df)
  
  # print(community)
  # pca_percentages <- (genus_pca$CA$eig/sum(genus_pca$CA$eig)) * 100
  # pca_percentages2 <- (genus_decontaminated_pca$CA$eig/sum(genus_decontaminated_pca$CA$eig)) * 100
  # print(pca_percentages)
  # print(pca_percentages2)

  for (myvar in discrete_variables){
    if (myvar == "Remote_Community") {next}
    # OTU, normal
    generate_pca(otu_pca, mymetadata = metadata_subset.df,
                 plot_height = 5, plot_width = 5,
                 # legend_x = -8, legend_y = 4,
                 # legend_x = -2, legend_y = 2,
                 point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
                 legend_title = myvar,
                 legend_cex = .5,
                 title_cex = .5,
                 plot_title = paste0(myvar,", Remote Community ", community),
                 # limits = c(-8,5,-6,5),
                 plot_spiders = F,
                 plot_ellipses = F,
                 plot_hulls = F,
                 use_shapes = T,
                 ellipse_border_width = .5,
                 include_legend = T,
                 label_ellipse = F, ellipse_label_size = .3,
                 colour_palette = my_colour_palette_15,
                 variable_to_plot = myvar, legend_cols = 1,
                 variable_colours_available = T,
                 num_top_species = 3,
                 plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
                 label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
                 specie_labeller_function = otu_relabeller_function,arrow_label_offset = 0,
                 filename = paste0("Result_figures/pca_plots/otu_within_community/otu_",myvar, "_Remote_Community_", community,"_pca.pdf"))
    
    # Genus, normal
    generate_pca(genus_pca, mymetadata = metadata_subset.df,
                 plot_height = 5, plot_width = 5,
                 # legend_x = -4, legend_y = 3,
                 # legend_x = -2, legend_y = 2,
                 point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
                 legend_title = myvar,
                 legend_cex = .5,
                 title_cex = .5,
                 plot_title = paste0(myvar,", Remote Community ", community),
                 # limits = c(-4,6,-4,4),
                 plot_spiders = F,
                 plot_ellipses = F,
                 plot_hulls = F,
                 use_shapes = T,
                 ellipse_border_width = .5,
                 include_legend = T,
                 label_ellipse = F, ellipse_label_size = .3,
                 colour_palette = my_colour_palette_15,
                 variable_to_plot = myvar, legend_cols = 1,
                 variable_colours_available = T,
                 num_top_species = 3,
                 plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
                 label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
                 specie_labeller_function = genus_relabeller_function,
                 arrow_label_offset = 0,
                 filename = paste0("Result_figures/pca_plots/genus_within_community/genus_",myvar, "_Remote_Community_", community,"_pca.pdf"))
    # OTU, decontaminated
    generate_pca(otu_decontaminated_pca, mymetadata = metadata_decontaminated_subset.df,
                 plot_height = 5, plot_width = 5,
                 # legend_x = -8, legend_y = 4,
                 # legend_x = -2, legend_y = 2,
                 point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
                 legend_title = myvar,
                 legend_cex = .5,
                 title_cex = .5,
                 plot_title = paste0(myvar,", Remote Community ", community),
                 # limits = c(-8,4,-8,4),
                 plot_spiders = F,
                 plot_ellipses = F,
                 plot_hulls = F,
                 use_shapes = T,
                 ellipse_border_width = .5,
                 include_legend = T,
                 label_ellipse = F, ellipse_label_size = .3,
                 colour_palette = my_colour_palette_15,
                 variable_to_plot = myvar, legend_cols = 1,
                 variable_colours_available = T,
                 num_top_species = 3,
                 plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
                 label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
                 specie_labeller_function = otu_relabeller_function,
                 arrow_label_offset = 0,
                 filename = paste0("Result_figures/pca_plots/otu_within_community_decontaminated/otu_decontaminated_",myvar, "_Remote_Community_", community,"_pca.pdf"))
    # Genus, decontaminated
    generate_pca(genus_decontaminated_pca, mymetadata = metadata_decontaminated_subset.df,
                 plot_height = 5, plot_width = 5,
                 # legend_x = -4, legend_y = 4,
                 # legend_x = -2, legend_y = 2,
                 point_size = .7, point_line_thickness = 0.3,point_alpha =.9,
                 legend_title = myvar,
                 legend_cex = .5,
                 title_cex = .5,
                 plot_title = paste0(myvar,", Remote Community ", community),
                 # limits = c(-4,6,-4,4),
                 plot_spiders = F,
                 plot_ellipses = F,
                 plot_hulls = F,
                 use_shapes = T,
                 ellipse_border_width = .5,
                 include_legend = T,
                 label_ellipse = F, ellipse_label_size = .3,
                 colour_palette = my_colour_palette_15,
                 variable_to_plot = myvar, legend_cols = 1,
                 variable_colours_available = T,
                 num_top_species = 3,
                 plot_arrows = T,arrow_alpha = .7, arrow_colour = "grey20",arrow_scalar = 1,arrow_thickness = .7,
                 label_arrows = T, arrow_label_size = .25, arrow_label_colour = "black", arrow_label_font_type = 1,
                 specie_labeller_function = genus_relabeller_function,
                 arrow_label_offset = 0,
                 filename = paste0("Result_figures/pca_plots/genus_within_community_decontaminated/genus_decontaminated_",myvar, "_Remote_Community_", community,"_pca.pdf"))    
  }
}



# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# PERMANOVA tests whether distance differ between groups.

# Genus, clr euclidean
print("Centred-log ratio transformed counts - Euclidean distance")
# otu_clr.m <- clr(otu.m)
# genus_clr.m <- clr(genus.m)
# otu_clr_decontaminated.m <- clr(otu_decontaminated.m)
# genus_clr_decontaminated.m <- clr(genus_decontaminated.m)
otu_permanova_results <- data.frame()
genus_permanova_results <- data.frame()
otu_decontaminated_permanova_results <- data.frame()
genus_decontaminated_permanova_results <- data.frame()

otu_within_community_permanova_results <- data.frame()
genus_within_community_permanova_results <- data.frame()
otu_within_community_decontaminated_permanova_results <- data.frame()
genus_within_community_decontaminated_permanova_results <- data.frame()

for (myvar in discrete_variables){
  print(myvar)
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  metadata_decontaminated_subset.df <- metadata_decontaminated.df[!is.na(metadata_decontaminated.df[,myvar]),]
  
  otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
  genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
  otu_clr_decontaminated_subset.m <- otu_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
  genus_clr_decontaminated_subset.m <- genus_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
  
  otu_permanova_results <- rbind(otu_permanova_results,run_permanova_custom(my_metadata = metadata_subset.df, 
                                                      my_formula = as.formula(paste0("t(otu_clr_subset.m)~", myvar)),
                                                      my_method = "euclidean",label = "CLR",permutations = 9999))
  
  genus_permanova_results <- rbind(genus_permanova_results,run_permanova_custom(my_metadata = metadata_subset.df, 
                                                      my_formula = as.formula(paste0("t(genus_clr_subset.m)~", myvar)),
                                                      my_method = "euclidean",label = "CLR",permutations = 9999))
  
  otu_decontaminated_permanova_results <- rbind(otu_decontaminated_permanova_results, run_permanova_custom(my_metadata = metadata_decontaminated_subset.df, 
                                                      my_formula = as.formula(paste0("t(otu_clr_decontaminated_subset.m)~", myvar)),
                                                      my_method = "euclidean",label = "CLR",permutations = 9999))
  
  genus_decontaminated_permanova_results <- rbind(genus_decontaminated_permanova_results, run_permanova_custom(my_metadata = metadata_decontaminated_subset.df, 
                                                        my_formula = as.formula(paste0("t(genus_clr_decontaminated_subset.m)~", myvar)),
                                                        my_method = "euclidean",label = "CLR",permutations = 9999))
  if (myvar == "Remote_Community") {next}
  for (community in unique(metadata.df$Remote_Community)){
    
    metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
    metadata_subset.df <- subset(metadata_subset.df, Remote_Community = community)
    metadata_decontaminated_subset.df <- metadata_decontaminated.df[!is.na(metadata_decontaminated.df[,myvar]),]
    metadata_decontaminated_subset.df <- subset(metadata_decontaminated_subset.df, Remote_Community = community)
    
    otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
    genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
    otu_clr_decontaminated_subset.m <- otu_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
    genus_clr_decontaminated_subset.m <- genus_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
    
    temp <- run_permanova_custom(my_metadata = metadata.df, 
                                 my_formula = as.formula(paste0("t(otu_clr_subset.m)~", myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 9999)
    temp$Remote_Community <- community
    otu_within_community_permanova_results <- rbind(otu_within_community_permanova_results, temp)
    
    temp <- run_permanova_custom(my_metadata = metadata.df, 
                         my_formula = as.formula(paste0("t(genus_clr_subset.m)~", myvar)),
                         my_method = "euclidean",label = "CLR",permutations = 9999)
    temp$Remote_Community <- community
    genus_within_community_permanova_results <- rbind(genus_within_community_permanova_results,temp)
    
    temp <- run_permanova_custom(my_metadata = metadata_decontaminated_subset.df, 
                                 my_formula = as.formula(paste0("t(otu_clr_decontaminated_subset.m)~", myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 9999)
    temp$Remote_Community <- community
    otu_within_community_decontaminated_permanova_results <- rbind(otu_within_community_decontaminated_permanova_results, temp)
    
    temp <- run_permanova_custom(my_metadata = metadata_decontaminated_subset.df, 
                                 my_formula = as.formula(paste0("t(genus_clr_decontaminated_subset.m)~", myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 9999)
    temp$Remote_Community <- community
    genus_within_community_decontaminated_permanova_results <- rbind(genus_within_community_decontaminated_permanova_results, temp)
  }
}

write.csv(otu_permanova_results, file = "Result_tables/stats_various/otu_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_permanova_results, file = "Result_tables/stats_various/genus_PERMANOVA.csv", row.names = F, quote = F)
write.csv(otu_decontaminated_permanova_results, file = "Result_tables/stats_various/otu_decontaminated_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_decontaminated_permanova_results, file = "Result_tables/stats_various/genus_decontaminated_PERMANOVA.csv", row.names = F, quote = F)

write.csv(otu_within_community_permanova_results, file = "Result_tables/stats_various/otu_within_community_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_within_community_permanova_results, file = "Result_tables/stats_various/genus_within_community_PERMANOVA.csv", row.names = F, quote = F)
write.csv(otu_decontaminated_permanova_results, file = "Result_tables/stats_various/otu_within_community_decontaminated_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_decontaminated_permanova_results, file = "Result_tables/stats_various/genus_within_community_decontaminated_PERMANOVA.csv", row.names = F, quote = F)

# run_permanova_custom(my_metadata = metadata.df, my_formula = as.formula(t(genus_clr.m)~Remote_Community),
# my_method = "euclidean",label = "CLR")

# myvar <- "Remote_Community"
# run_permanova_custom(my_metadata = metadata.df, my_formula = as.formula(t(otu_clr.m)~get(myvar)),
#                      my_method = "euclidean",label = "CLR")
# 
# run_permanova_custom(my_metadata = metadata.df, my_formula = as.formula(t(otu_clr.m)~OM_6mo),
#                      my_method = "euclidean",label = "CLR")
# 
# run_permanova_custom(my_metadata = metadata.df, my_formula = as.formula(paste0("t(otu_clr.m)~", myvar)),
#                      my_method = "euclidean",label = "CLR")
# write.csv(permanova_results, file = "Result_tables/stats_various/PERMANOVA.csv", row.names = F, quote = F)



# ---------------------------------------------
# PERMDISP (betadisper)
# See: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
# "Test the homogeneity of within-group multivariate dispersions on the basis of any resemblance measure."

# temp <- with(metadata.df, betadisper(vegdist(t(genus_clr.m), method = "euclidean"), group = Remote_Community))

# temp <- permutest(temp, permutations = 999, parallel = 2)
# plot(temp, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
# boxplot(temp, main = "", xlab = "")
# 
# temp <- run_permdisp_custom(metadata.df, 
#                     my_data = genus_clr.m,
#                     my_group = "Remote_Community",
#                     my_method = "euclidean",
#                     permutations = 999, label = NULL)

otu_permdisp_results <- data.frame()
genus_permdisp_results <- data.frame()
otu_decontaminated_permdisp_results <- data.frame()
genus_decontaminated_permdisp_results <- data.frame()

otu_within_community_permdisp_results <- data.frame()
genus_within_community_permdisp_results <- data.frame()
otu_within_community_decontaminated_permdisp_results <- data.frame()
genus_within_community_decontaminated_permdisp_results <- data.frame()


# source("code/helper_functions.R")

for (myvar in discrete_variables){
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  metadata_decontaminated_subset.df <- metadata_decontaminated.df[!is.na(metadata_decontaminated.df[,myvar]),]
  
  otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
  genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
  otu_clr_decontaminated_subset.m <- otu_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
  genus_clr_decontaminated_subset.m <- genus_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
  
  otu_permdisp_results <- rbind(otu_permdisp_results, run_permdisp_custom(my_metadata = metadata_subset.df, 
                                                                          my_data = otu_clr_subset.m,
                                                                          my_group = myvar,
                                                                          my_method = "euclidean",
                                                                          permutations = 9999,
                                                                          label = "CLR"))
  
  genus_permdisp_results <- rbind(genus_permdisp_results, run_permdisp_custom(my_metadata = metadata_subset.df, 
                                                                              my_data = genus_clr_subset.m,
                                                                              my_group = myvar,
                                                                              my_method = "euclidean",
                                                                              permutations = 9999,
                                                                              label = "CLR"))
  
  otu_decontaminated_permdisp_results <- rbind(otu_decontaminated_permdisp_results, run_permdisp_custom(my_metadata = metadata_decontaminated_subset.df, 
                                                                                                        my_data = otu_clr_decontaminated_subset.m,
                                                                                                        my_group = myvar,
                                                                                                        my_method = "euclidean",
                                                                                                        permutations = 9999,
                                                                                                        label = "CLR"))
  
  genus_decontaminated_permdisp_results <- rbind(genus_decontaminated_permdisp_results, run_permdisp_custom(my_metadata = metadata_decontaminated_subset.df, 
                                                                                                            my_data = genus_clr_decontaminated_subset.m,
                                                                                                            my_group = myvar,
                                                                                                            my_method = "euclidean",
                                                                                                            permutations = 9999,
                                                                                                            label = "CLR"))
  if (myvar == "Remote_Community") {next}
  for (community in unique(metadata.df$Remote_Community)){
    
    metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
    metadata_subset.df <- subset(metadata_subset.df, Remote_Community == community)
    metadata_decontaminated_subset.df <- metadata_decontaminated.df[!is.na(metadata_decontaminated.df[,myvar]),]
    metadata_decontaminated_subset.df <- subset(metadata_decontaminated_subset.df, Remote_Community = community)
    
    otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
    genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
    otu_clr_decontaminated_subset.m <- otu_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
    genus_clr_decontaminated_subset.m <- genus_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
    
    temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                                my_data = otu_clr_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 9999,
                                label = "CLR")
    temp$Remote_Community <- community
    otu_within_community_permdisp_results <- rbind(otu_within_community_permdisp_results, temp)
    
    temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                                my_data = genus_clr_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 9999,
                                label = "CLR")
    temp$Remote_Community <- community
    genus_within_community_permdisp_results <- rbind(genus_within_community_permdisp_results, temp)
    
    temp <- run_permdisp_custom(my_metadata = metadata_decontaminated_subset.df, 
                                my_data = otu_clr_decontaminated_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 9999,
                                label = "CLR")
    temp$Remote_Community <- community
    otu_within_community_decontaminated_permdisp_results <- rbind(otu_within_community_decontaminated_permdisp_results, temp)
    
    temp <- run_permdisp_custom(my_metadata = metadata_decontaminated_subset.df, 
                                my_data = genus_clr_decontaminated_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 9999,
                                label = "CLR")
    temp$Remote_Community <- community
    genus_within_community_decontaminated_permdisp_results <- rbind(genus_within_community_decontaminated_permdisp_results, temp)
  }
  
}
  

write.csv(otu_permdisp_results, file = "Result_tables/stats_various/otu_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_permdisp_results, file = "Result_tables/stats_various/genus_PERMDISP.csv", row.names = F, quote = F)
write.csv(otu_decontaminated_permdisp_results, file = "Result_tables/stats_various/otu_decontaminated_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_decontaminated_permdisp_results, file = "Result_tables/stats_various/genus_decontaminated_PERMDISP.csv", row.names = F, quote = F)

write.csv(otu_within_community_permdisp_results, file = "Result_tables/stats_various/otu_within_community_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_within_community_permdisp_results, file = "Result_tables/stats_various/genus_within_community_PERMDISP.csv", row.names = F, quote = F)
write.csv(otu_decontaminated_permdisp_results, file = "Result_tables/stats_various/otu_within_community_decontaminated_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_decontaminated_permdisp_results, file = "Result_tables/stats_various/genus_within_community_decontaminated_PERMDISP.csv", row.names = F, quote = F)




# ---------------------------------------------
# Takes awhile to calculate, uncomment to run
# ANOSIM tests whether distances between groups are greater than within groups.
# "Nonparametric procedure for testing the hypothesis of no difference between two or more groups of entities
# based on permutation test of among- and within-group similarities"
# R = 1 when all pairs of samples within groups are more similar than to any pair of samples from different groups
# R = 0 expected value under the null model that among-and within- group dissimilarities are the same on average.

# "If you have very different group sizes, you may consider analysis of similarities (ANOSIM) instead of PERMANOVA. 
# This test does not assume equal group variances."

otu_anosim_results <- data.frame()
genus_anosim_results <- data.frame()
otu_decontaminated_anosim_results <- data.frame()
genus_decontaminated_anosim_results <- data.frame()

otu_within_community_anosim_results <- data.frame()
genus_within_community_anosim_results <- data.frame()
otu_within_community_decontaminated_anosim_results <- data.frame()
genus_within_community_decontaminated_anosim_results <- data.frame()

for (myvar in discrete_variables){
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  metadata_decontaminated_subset.df <- metadata_decontaminated.df[!is.na(metadata_decontaminated.df[,myvar]),]
  
  otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
  genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
  otu_clr_decontaminated_subset.m <- otu_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
  genus_clr_decontaminated_subset.m <- genus_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
  
  otu_anosim_results <- rbind(otu_anosim_results, run_anosim_custom(my_metadata = metadata_subset.df, 
                                                                    my_data = otu_clr_subset.m,
                                                                    my_group = myvar,
                                                                    my_method = "euclidean",
                                                                    permutations = 9999,
                                                                    label = "CLR"))
  
  genus_anosim_results <- rbind(genus_anosim_results, run_anosim_custom(my_metadata = metadata_subset.df, 
                                                                        my_data = genus_clr_subset.m,
                                                                        my_group = myvar,
                                                                        my_method = "euclidean",
                                                                        permutations = 9999,
                                                                        label = "CLR"))
  
  
  otu_decontaminated_anosim_results <- rbind(otu_decontaminated_anosim_results, run_anosim_custom(my_metadata = metadata_decontaminated_subset.df, 
                                                                                                  my_data = otu_clr_decontaminated_subset.m,
                                                                                                  my_group = myvar,
                                                                                                  my_method = "euclidean",
                                                                                                  permutations = 9999,
                                                                                                  label = "CLR"))
  
  genus_decontaminated_anosim_results <- rbind(genus_decontaminated_anosim_results, run_anosim_custom(my_metadata = metadata_decontaminated_subset.df, 
                                                                                                      my_data = genus_clr_decontaminated_subset.m,
                                                                                                      my_group = myvar,
                                                                                                      my_method = "euclidean",
                                                                                                      permutations = 9999,
                                                                                                      label = "CLR"))
  
  if (myvar == "Remote_Community") {next}
  for (community in unique(metadata.df$Remote_Community)){
    metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
    metadata_subset.df <- subset(metadata_subset.df, Remote_Community == community)
    metadata_decontaminated_subset.df <- metadata_decontaminated.df[!is.na(metadata_decontaminated.df[,myvar]),]
    metadata_decontaminated_subset.df <- subset(metadata_decontaminated_subset.df, Remote_Community = community)
    
    otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
    genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
    otu_clr_decontaminated_subset.m <- otu_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
    genus_clr_decontaminated_subset.m <- genus_clr_decontaminated.m[,rownames(metadata_decontaminated_subset.df)]
    
    temp <- run_anosim_custom(my_metadata = metadata_subset.df, 
                              my_data = otu_clr_subset.m,
                              my_group = myvar,
                              my_method = "euclidean",
                              permutations = 9999,
                              label = "CLR")
    temp$Remote_Community <- community
    otu_within_community_anosim_results <- rbind(otu_within_community_anosim_results, temp)
    
    temp <- run_anosim_custom(my_metadata = metadata_subset.df, 
                              my_data = genus_clr_subset.m,
                              my_group = myvar,
                              my_method = "euclidean",
                              permutations = 9999,
                              label = "CLR")
    temp$Remote_Community <- community
    genus_within_community_anosim_results <- rbind(genus_within_community_anosim_results, temp)
    
    temp <- run_anosim_custom(my_metadata = metadata_decontaminated_subset.df, 
                              my_data = otu_clr_decontaminated_subset.m,
                              my_group = myvar,
                              my_method = "euclidean",
                              permutations = 9999,
                              label = "CLR")
    temp$Remote_Community <- community
    otu_within_community_decontaminated_anosim_results <- rbind(otu_within_community_decontaminated_anosim_results, temp)
    
    temp <- run_anosim_custom(my_metadata = metadata_decontaminated_subset.df, 
                              my_data = genus_clr_decontaminated_subset.m,
                              my_group = myvar,
                              my_method = "euclidean",
                              permutations = 9999,
                              label = "CLR")
    temp$Remote_Community <- community
    genus_within_community_decontaminated_anosim_results <- rbind(genus_within_community_decontaminated_anosim_results, temp)
  }
}


write.csv(otu_anosim_results, file = "Result_tables/stats_various/otu_ANOSIM.csv", row.names = F, quote = F)
write.csv(genus_anosim_results, file = "Result_tables/stats_various/genus_ANOSIM.csv", row.names = F, quote = F)
write.csv(otu_decontaminated_anosim_results, file = "Result_tables/stats_various/otu_decontaminated_ANOSIM.csv", row.names = F, quote = F)
write.csv(genus_decontaminated_anosim_results, file = "Result_tables/stats_various/genus_decontaminated_ANOSIM.csv", row.names = F, quote = F)

write.csv(otu_within_community_anosim_results, file = "Result_tables/stats_various/otu_within_community_ANOSIM.csv", row.names = F, quote = F)
write.csv(genus_within_community_anosim_results, file = "Result_tables/stats_various/genus_within_community_ANOSIM.csv", row.names = F, quote = F)
write.csv(otu_decontaminated_anosim_results, file = "Result_tables/stats_various/otu_within_community_decontaminated_ANOSIM.csv", row.names = F, quote = F)
write.csv(genus_decontaminated_anosim_results, file = "Result_tables/stats_various/genus_within_community_decontaminated_ANOSIM.csv", row.names = F, quote = F)



