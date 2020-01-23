# Analysis specific to the cultured species


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


setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project")
source("code/helper_functions.R")


# Load metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv")
metadata_decontaminated.df <- read.csv("Result_tables/other/processed_metadata_decontaminated.csv")


# Load Genus level data
# genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv")
# genus_data_decontaminated.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata_decontaminated.csv")


# Plot and calculate correlation between presence/absence of Haemophilus_influenzae, Moraxella_catarrhalis and Streptococcus_pneumoniae
# against abundances and/or presence/absence of all other taxa
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv")


haemophilus.df <- genus_rel.df[grep("g__Haemophilus", genus_rel.df$taxonomy_genus),]
haemophilus.df <- melt(haemophilus.df, id.vars = c("taxonomy_genus"), variable.name = "Sample", value.name = "Relative_abundance")
haemophilus.df <- left_join(haemophilus.df, metadata.df,  by = c("Sample" = "Index"))
haemophilus.df <- haemophilus.df[!is.na(haemophilus.df$Haemophilus_influenzae),]

moraxella.df <- genus_rel.df[grep("g__Moraxella", genus_rel.df$taxonomy_genus),]
moraxella.df <- melt(moraxella.df, id.vars = c("taxonomy_genus"), variable.name = "Sample", value.name = "Relative_abundance")
moraxella.df <- left_join(moraxella.df, metadata.df,  by = c("Sample" = "Index"))
moraxella.df <- moraxella.df[!is.na(moraxella.df$Moraxella_catarrhalis),]

streptococcus.df <- genus_rel.df[grep("g__Streptococcus", genus_rel.df$taxonomy_genus),]
streptococcus.df <- melt(streptococcus.df, id.vars = c("taxonomy_genus"), variable.name = "Sample", value.name = "Relative_abundance")
streptococcus.df <- left_join(streptococcus.df, metadata.df,  by = c("Sample" = "Index"))
streptococcus.df <- streptococcus.df[!is.na(streptococcus.df$Streptococcus_pneumoniae),]

# colour_palette <- with(unique(temp[,c("Haemophilus_influenzae", "Haemophilus_influenzae_colour")]), setNames(as.character(Haemophilus_influenzae_colour), Haemophilus_influenzae))

ggplot(haemophilus.df, aes(x = factor(Haemophilus_influenzae), y = Relative_abundance, fill = Haemophilus_influenzae_colour)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = .5) + 
  # xlab("Haemophilus influenzae present / absent (culture)") +
  ylab("Relative abundance") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_fill_identity() +
  common_theme
  
ggplot(moraxella.df, aes(x = factor(Moraxella_catarrhalis), y = Relative_abundance, fill = Moraxella_catarrhalis_colour)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = .5) + 
  # xlab("Haemophilus influenzae present / absent (culture)") +
  ylab("Relative abundance") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_fill_identity() +
  common_theme

ggplot(streptococcus.df, aes(x = factor(Streptococcus_pneumoniae), y = Relative_abundance, fill = Streptococcus_pneumoniae_colour)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = .5) + 
  # xlab("Haemophilus influenzae present / absent (culture)") +
  ylab("Relative abundance") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = .1)) +
  scale_fill_identity() +
  common_theme

