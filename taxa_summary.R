library(reshape2)
library(dplyr)


setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_16S_project/")
source("code/helper_functions.R")


metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", header = T)
<<<<<<< HEAD

genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)


genus_rel.df <- melt(genus_rel.df, variable.name = "Sample", value.name = "Relative_abundance")
=======
# metadata.df <- read.table("data/metadata.tsv", header = T, sep = "\t",row.names = F)
# metadata.df$Otitis_Status
# metadata.df$Index <- metadata.df$Sequence_file_ID_clean
# metadata.df <- metadata.df[metadata.df$Index %in% colnames(genus_rel.df),]

genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)

temp <- df2matrix(genus_rel.df)["d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Erwiniaceae;g__Pantoea",
             subset(metadata.df, Nose == "Serous")$Index]
mean(temp)
median(temp)
max(temp)
# abundance_summary(genus_rel.df)
genus_rel.df <- melt(genus_rel.df, variable.name = "Sample", value.name = "Relative_abundance")

>>>>>>> ac8d1af8b0bcc3ccad540a152432b91af0083eb7
genus_rel.df <- left_join(genus_rel.df, metadata.df, by = c("Sample" = "Index"))


discrete_variables <- c("Community", "Nose", "Otitis_Status", "Season", "No_peop_res_discrete")

out <- data.frame()
for (variable in discrete_variables){
  temp <- genus_rel.df %>% 
    dplyr::group_by(get(variable), taxonomy_genus) %>%
    dplyr::mutate(Samples_total_in_group = n_distinct(Sample)) %>% # number of unique samples/index
    # dplyr::group_by_at(1,.drop = F) %>% 
    dplyr::mutate(Number_of_samples_taxa_present = sum(Relative_abundance > 0, na.rm = TRUE)) %>%
    dplyr::summarise(
      Number_of_samples_taxa_present = max(Number_of_samples_taxa_present),
      Samples_total_in_group = max(Samples_total_in_group),
      # Mean_relative_abundance2 = round(mean(sum(Relative_abundance)/max(Samples_total)), 5),
      Min_relative_abundance = round(min(Relative_abundance),5),
      Max_relative_abundance = round(max(Relative_abundance),5),
      Mean_relative_abundance = round(mean(Relative_abundance), 5),
      Stdev_relative_abundance = round(sd(Relative_abundance), 5),
      Median_relative_abundance = round(median(Relative_abundance), 5),
    ) %>% 
    as.data.frame()
  colnames(temp)[1] <- "Group"
  temp$Variable <-  variable
  temp <- 
    temp %>%
    select("Variable", everything())
  temp <- temp %>% dplyr::arrange(Variable,Group, desc(Mean_relative_abundance))
  out <- rbind(out,temp)
}

write.csv(x = out, file = "Result_tables/abundance_analysis_tables/genus_variable_summary.csv", quote = F, row.names = F)

<<<<<<< HEAD
=======
genus_rel.df

>>>>>>> ac8d1af8b0bcc3ccad540a152432b91af0083eb7

