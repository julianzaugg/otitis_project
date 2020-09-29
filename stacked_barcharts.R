library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")



# Load metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv")

# Ensure rownames are the Index
rownames(metadata.df) <- metadata.df$Index

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load abundance data
otu_rel.df <- read.csv("Result_tables/relative_abundance_tables/OTU_relative_abundances.csv", header = T)
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)

# First create a global palette of the top taxa (more than we need though not everything)
palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
palette_20b <- c("#d05c2a","#a9b635","#b3a95e","#3fbcc1","#397f50","#e19073","#68752b","#d790ce","#7663cd","#d99c34","#97518d","#cf4248","#63c34e","#5fbf85","#d94983","#59993a","#658acc","#ad5262","#c653bd","#95642d")
palette_20c <- c("#a84b54","#a2b432","#c057c5","#6f64cf","#5aae36","#527c2e","#826729","#6690cf","#7fb875","#d68c62","#a8a352","#e484a0","#d2408e","#d49731","#ca5729","#3a8862","#49c1b5","#ab6daa","#dc3f53","#4bc06d")
palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
palette_100_distinct <- c("#b2e54d","#acb5ff","#7ce5e2","#ab3d00","#a4e75a","#699164","#00b9e2","#5a3e6b","#dbcf00","#26511d","#00b979","#f7d152","#008885","#0072b9","#493d85","#6a67fe","#01cdbb","#f8a7ff","#ff54ac","#a8b800","#e675ff","#3080ff","#edd37c","#15ca42","#d1db7b","#014cb1","#01d179","#800489","#ffa27c","#b70094","#ffcc6a","#c1db03","#00ac42","#ffb971","#e00097","#5aef87","#fac0ff","#ff74b9","#f1053b","#1f5a00","#898ac0","#7fbf00","#016cc8","#514804","#00578b","#abc994","#4ba300","#ffad4d","#9720ba","#cd27bf","#ff68d7","#01b0b1","#bd6700","#87a1ff","#a20019","#842708","#ea0090","#006e04","#facc9c","#8ae6c3","#b28400","#4d8cff","#d86400","#84223c","#918b55","#234198","#ff8242","#6d3260","#e31db5","#7c7f00","#3a8500","#70ec9f","#4a5700","#a065fd","#2f41d0","#f0a300","#bc0042","#c184ff","#93e973","#6d4400","#58699b","#ff8d1c","#00997f","#ff6655","#88066c","#b6e190","#ff2899","#ce2a00","#ff2b91","#a20062","#7e2d23","#ff7c88","#783403","#009f49","#a29600","#017a67","#ffaeb6","#4e329e","#00722f","#5e23a7")
top_genus <- abundance_summary(genus_rel.df,20)$taxonomy_genus

genus_palette <- setNames(c(palette_20c[1:length(top_genus)], "grey"), c(top_genus, "Other"))

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# All samples (grouped by gold star)

# Get list of samples from each group
gold_star_healthy_samples.v <- metadata.df[metadata.df$Gold_Star == "Healthy",]$Index
gold_star_not_healthy_samples.v <- metadata.df[metadata.df$Gold_Star == "Hx/Current OM/URTI",]$Index

# Generate summaries
gold_star_healthy_samples_summary.df <- abundance_summary(genus_rel.df[,c("taxonomy_genus", gold_star_healthy_samples.v)],8)
gold_star_not_healthy_samples_summary.df <- abundance_summary(genus_rel.df[,c("taxonomy_genus", gold_star_not_healthy_samples.v)],8)

# Calculate top genus
gold_star_top_genus <- 
  unique(rbind(gold_star_healthy_samples_summary.df,gold_star_not_healthy_samples_summary.df)$taxonomy_genus)

# Collapse abundance dataframe
gold_star_genus.df <- collapse_abundance_dataframe(genus_rel.df, gold_star_top_genus)

# Combine with metadata
gold_star_genus.df <- left_join(gold_star_genus.df, metadata.df[,c("Index", "Gold_Star")], by = c("Sample"= "Index"))

# Get order of top taxa by mean abundance
temp <- gold_star_genus.df[,c("taxonomy_genus","Sample", "Relative_abundance")]
temp <- df2matrix(dcast(temp,taxonomy_genus~Sample))
ordered_taxa <- names(sort(apply(temp,1,mean),decreasing = T))
ordered_taxa <- c(ordered_taxa[!ordered_taxa == "Other"],"Other")

# Order dataframe by abundance
gold_star_genus.df <- arrange(gold_star_genus.df, dplyr::desc(Relative_abundance),taxonomy_genus)
# gold_star_genus.df <- arrange(gold_star_genus.df, taxonomy_genus, dplyr::desc(Relative_abundance))

# Factorise sample
gold_star_genus.df$Sample <- factor(gold_star_genus.df$Sample, 
                                    levels = unique(as.character(gold_star_genus.df$Sample)))
# Factorise taxonomy
gold_star_genus.df$taxonomy_genus <- factor(gold_star_genus.df$taxonomy_genus, levels = ordered_taxa)


gold_star_genus.df$Label <- unlist(lapply(gold_star_genus.df$taxonomy_genus,
                                          function(x) ifelse(grepl(";", x),
                                                             first_resolved_taxonomy(as.character(x)),
                                                             as.character(x)))
                                   )

labels <- unlist(lapply(as.character(unique(gold_star_genus.df$taxonomy_genus)),function(x) ifelse(grepl(";", x),
                                                                                         first_resolved_taxonomy(as.character(x)),
                                                                                         as.character(x))))

label_list <- setNames(labels,as.character(unique(gold_star_genus.df$taxonomy_genus)))
# my_labels <- unique(setNames(gold_star_genus.df$Label,as.character(gold_star_genus.df$taxonomy_genus)))

gold_star_genus_all_samples_plot <- 
  ggplot(gold_star_genus.df, aes(x = Sample, 
                                 y = Relative_abundance*100, 
                                 fill = taxonomy_genus))+
  coord_flip() +
  geom_bar(stat = "identity", colour = "black", lwd = 0.1, width = 0.75) +
  xlab("Sample") + 
  ylab("Relative abundance %") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = genus_palette, 
                    labels = label_list,
                    name = "Taxonomy") +
  facet_wrap(~Gold_Star, scales = "free_y") +
  theme(
    # axis.text.x = element_text(angle = 90, size = 6,hjust = 0, vjust = 0.5),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 10,face = "bold"),
    )
ggsave(
  filename = "Result_figures/abundance_analysis_plots/gold_star_genus_all_samples_stacked_barchart.pdf",
  plot = gold_star_genus_all_samples_plot,
  height = 30,
  width = 20,
  units = "cm")

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# All samples where Otitis Status is either Never OM or HxOM


# Get list of samples from each group
otitis_status_never_om_samples.v <- metadata.df[metadata.df$Otitis_Status == "Never OM",]$Index
otitis_status_hxom_samples.v <- metadata.df[metadata.df$Otitis_Status == "HxOM",]$Index

# Generate summaries
otitis_never_om_samples_summary.df <- abundance_summary(genus_rel.df[,c("taxonomy_genus", otitis_status_never_om_samples.v)],8)
otitis_status_hxom_samples_summary.df <- abundance_summary(genus_rel.df[,c("taxonomy_genus", otitis_status_hxom_samples.v)],8)

# Calculate top genus
otitis_status_top_genus <- 
  unique(rbind(otitis_never_om_samples_summary.df,otitis_status_hxom_samples_summary.df)$taxonomy_genus)

# Collapse abundance dataframe
otitis_status_genus.df <- collapse_abundance_dataframe(genus_rel.df, otitis_status_top_genus)

# Combine with metadata
otitis_status_genus.df <- left_join(otitis_status_genus.df, metadata.df[,c("Index", "Otitis_Status")], by = c("Sample"= "Index"))
otitis_status_genus.df <- otitis_status_genus.df[otitis_status_genus.df$Sample %in% c(otitis_status_never_om_samples.v, otitis_status_hxom_samples.v),]

# Get order of top taxa by mean abundance
temp <- otitis_status_genus.df[,c("taxonomy_genus","Sample", "Relative_abundance")]
temp <- df2matrix(dcast(temp,taxonomy_genus~Sample))
ordered_taxa <- names(sort(apply(temp,1,mean),decreasing = T))
ordered_taxa <- c(ordered_taxa[!ordered_taxa == "Other"],"Other")

# Order dataframe by abundance
otitis_status_genus.df <- arrange(otitis_status_genus.df, dplyr::desc(Relative_abundance),taxonomy_genus)

# Factorise sample
otitis_status_genus.df$Sample <- factor(otitis_status_genus.df$Sample, 
                                    levels = unique(as.character(otitis_status_genus.df$Sample)))
# Factorise taxonomy
otitis_status_genus.df$taxonomy_genus <- factor(otitis_status_genus.df$taxonomy_genus, levels = ordered_taxa)


otitis_status_genus.df$Label <- unlist(lapply(otitis_status_genus.df$taxonomy_genus,
                                          function(x) ifelse(grepl(";", x),
                                                             first_resolved_taxonomy(as.character(x)),
                                                             as.character(x)))
)

labels <- unlist(lapply(as.character(unique(otitis_status_genus.df$taxonomy_genus)),function(x) ifelse(grepl(";", x),
                                                                                                   first_resolved_taxonomy(as.character(x)),
                                                                                                   as.character(x))))

label_list <- setNames(labels,as.character(unique(otitis_status_genus.df$taxonomy_genus)))
otitis_status_genus.df$Otitis_Status <- factor(otitis_status_genus.df$Otitis_Status, levels = c("Never OM", "HxOM"))
otitis_status_genus_all_samples_plot <- 
  ggplot(otitis_status_genus.df, aes(x = Sample, 
                                 y = Relative_abundance*100, 
                                 fill = taxonomy_genus))+
  coord_flip() +
  geom_bar(stat = "identity", colour = "black", lwd = 0.1, width = 0.75) +
  xlab("Sample") + 
  ylab("Relative abundance %") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = genus_palette, 
                    labels = label_list,
                    name = "Taxonomy") +
  facet_wrap(~Otitis_Status, scales = "free_y") +
  theme(
    # axis.text.x = element_text(angle = 90, size = 8,hjust = 0, vjust = 0.5),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 10,face = "bold"),
  )
ggsave(
  filename = "Result_figures/abundance_analysis_plots/otitis_status_genus_all_samples_plot_stacked_barchart.pdf",
  plot = otitis_status_genus_all_samples_plot,
  height = 30,
  width = 20,
  units = "cm")

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
