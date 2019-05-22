library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)



# install.packages("lemon")
library(lemon)



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


############################################################
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
my_colour_pallete_12_soft <-c("#9E788F","#4C5B61","#678D58","#AD5233","#A0A083","#4D456A","#588578","#D0AC4C","#2A7BA0","#931621", "#c75a93", "#7c7731")
############################################################
# Set the working directory
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/otitis_project/")



plot_genus_deseq <- function(deseq_result_table, facet_plot = T, limit_to = NULL, title = NULL, pallete = my_colour_pallete_30_distinct,
                             include_grid = F,point_size =1 ){
  internal_table <- deseq_result_table
  # Separate into different taxa levels
  internal_table <- tidyr::separate(internal_table,col = "Taxonomy", into = c("Domain", "Phylum", "Class", "Order","Family", "Genus"), sep = ";")
  # Combine the groups into single variable
  internal_table$Group_1_Group_2 <- with(internal_table, paste0(Group_1, "_vs_", Group_2))
  
  if (!is.null(limit_to)){
    internal_table <- internal_table[grep(limit_to, internal_table$Group_1_Group_2),]
  }
  
  # Order by the log fold change
  internal_table <- internal_table[order(internal_table$log2FoldChange),]
  # if (facet_plot == T){
  #   internal_table <- internal_table[order(internal_table$Group_1_Group_2),]
  # }
  
  # Get the maximum log fold value for each genus
  x <- with(internal_table, tapply(log2FoldChange, Genus, function(x) max(x)))
  x <- sort(x, T)
  
  # Order the levels of the genus by the value
  internal_table$Genus <- factor(internal_table$Genus, levels = names(x))
  
  myplot <- ggplot(internal_table, aes(x= Genus, y= log2FoldChange))
  
  if (include_grid == T){
    myplot <- myplot + background_grid(major = "x", colour.major = "grey", size.major = .1)
                                       # minor = "x", colour.minor = "red", size.minor = .1)
    # only_minor
  }
  myplot <- myplot + 
    geom_hline(yintercept = 0, linetype = "dashed", size =.3, colour = "grey") +
    geom_point(aes(fill = Class), shape = 21,stroke = .1,size = point_size, color = "black", alpha= .8) +
    # geom_point(aes(fill = Class, color = Sample_Type), shape = 21,stroke = .1,size = 2, alpha= .8) +
    scale_fill_manual(values = pallete) +
    theme(axis.text.x = element_text(angle = 90,size = 6,vjust = .5),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 6),
          plot.title = element_text(size = 8),
          strip.text = element_text(size= 6),
          strip.background = element_blank(),
          legend.title = element_text(size = 6,face = "bold"),
          legend.text = element_text(size = 6),
          legend.title.align = 0.5,
          legend.margin = margin(c(2,2,2,2)),
          legend.key.height=unit(.4,"cm"),
          axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          axis.ticks = element_line(size = 0.3)
          )
  
  if (!is.null(title)){
    myplot <- myplot + ggtitle(title)
  }
  
  if (facet_plot == T){
    myplot <- myplot + 
      # facet_wrap(~Group_1_Group_2, ncol = 1) 
      facet_rep_wrap(~Group_1_Group_2, ncol = 1, repeat.tick.labels = F) # from 'lemon'. To keep ticks in facets.
      # geom_hline(yintercept=max(internal_table$log2FoldCchange))
      # scale_x_discrete(breaks = levels(internal_table$Genus))
  }
  return(myplot)
}

# Load the deseq data. Assume that there is a Group_1 and Group_2 column

# --------------------------------------------------------------------------------

# Define the discrete variables
discrete_variables <- c("Remote_Community","Otitis_status","Gold_Star","OM_6mo","Type_OM","Season",
                        "Nose","Otitis_status_OM_6mo", "Remote_Community_Otitis_status", "OM_6mo_Type_OM","Remote_Community_Season")

otu_groups <- read.csv("Result_tables/DESeq_results/OTU_variable_groups.csv", header = T)
otu_groups <- otu_groups[abs(otu_groups$log2FoldChange) > 1,]

# Otitis_status
data_subset <- subset(otu_groups, Variable == "Otitis_status")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Otitis_status", 
                   pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Otitis_status.pdf"),plot = myplot, width = 8, height = 8, units = "cm")   

# OM_6mo
data_subset <- subset(otu_groups, Variable == "OM_6mo")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "OM_6mo", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/OM_6mo.pdf"),plot = myplot, width = 8, height = 8, units = "cm")   

# Otitis_status_OM_6mo
data_subset <- subset(otu_groups, Variable == "Otitis_status_OM_6mo")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Otitis_status_OM_6mo", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Otitis_status_OM_6mo.pdf"),plot = myplot, width = 8, height = 16, units = "cm")   

# Remote_Community
data_subset <- subset(otu_groups, Variable == "Remote_Community")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Remote_Community", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Remote_Community.pdf"),plot = myplot, width = 9, height = 8, units = "cm")   

# Remote_Community_Otitis_status
data_subset <- subset(otu_groups, Variable == "Remote_Community_Otitis_status")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Remote_Community_Otitis_status", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Remote_Community_Otitis_status.pdf"),plot = myplot, width = 9, height = 18, units = "cm")   

# Gold_star
data_subset <- subset(otu_groups, Variable == "Gold_Star")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Gold_Star", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Gold_Star.pdf"),plot = myplot, width = 9, height = 8, units = "cm")   


# Type_OM
data_subset <- subset(otu_groups, Variable == "Type_OM")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = "(0|1|2)_vs_(0|1|2)", title = "Type_OM", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Type_OM.pdf"),plot = myplot, width = 9, height = 12, units = "cm")   


# OM_6mo_Type_OM
data_subset <- subset(otu_groups, Variable == "OM_6mo_Type_OM")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = "1_(0|2)_vs_1_(0|2)", title = "OM_6mo_Type_OM", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/OM_6mo_Type_OM.pdf"),plot = myplot, width = 9, height = 8, units = "cm")   

# Season 
data_subset <- subset(otu_groups, Variable == "Season")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Season", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Season.pdf"),plot = myplot, width = 9, height = 14, units = "cm")   

# Remote_Community_Season
data_subset <- subset(otu_groups, Variable == "Remote_Community_Season")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Remote_Community_Season", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Remote_Community_Season.pdf"),plot = myplot, width = 9, height = 30, units = "cm")   

# Nose
data_subset <- subset(otu_groups, Variable == "Nose")
myplot <- plot_genus_deseq(data_subset,facet_plot = T,limit_to = ".", title = "Nose", 
                           pallete = my_colour_pallete_12_soft,include_grid = T,point_size = 2)
myplot
ggsave(filename = paste0("Result_figures/DESeq_plots/Nose.pdf"),plot = myplot, width = 9, height = 22, units = "cm")   





data_subset <- subset(otu_groups, Variable == "OM_6mo_Type_OM")
myplot <- 
  plot_genus_deseq(data_subset,facet_plot = T,
                   limit_to = ".", 
                   title = "OM_6mo_Type_OM", 
                   pallete = my_colour_pallete_12_soft,
                   include_grid = T) 
myplot






for (myvar in discrete_variables){
  data_subset <- subset(otu_groups, Variable == myvar)
  myplot <- 
    plot_genus_deseq(data_subset,facet_plot = T,
                     limit_to = ".", 
                     title = myvar, 
                     pallete = my_colour_pallete_12_soft,
                     include_grid = T) 
  
  if(length(unique(with(data_subset, paste0(Group_1,"_", Group_2)))) == 3){
    ggsave(filename = paste0("Result_figures/DESeq_plots/",myvar, ".pdf"),
           plot = myplot,
           width = 10,
           height = 12,
           units = "cm")   
  } else if(length(unique(with(data_subset, paste0(Group_1,"_", Group_2)))) > 3){
    ggsave(filename = paste0("Result_figures/DESeq_plots/",myvar, ".pdf"),
           plot = myplot,
           width = 10,
           height = 22,
           units = "cm")   
    
  } else{
    ggsave(filename = paste0("Result_figures/DESeq_plots/",myvar, ".pdf"),
           plot = myplot,
           width = 12,
           height = 8,
           units = "cm")    
  }
  
}



otu_groups <- read.csv("Result_tables/DESeq_results/OTU_variable_groups_within_community.csv", header = T)
otu_groups <- otu_groups[abs(otu_groups$log2FoldChange) > 1,]


for (community in unique(otu_groups$Remote_Community)){
  
  for (myvar in discrete_variables){
    data_subset <- subset(otu_groups, Variable == myvar & Remote_Community == community)
    if (grepl("Remote_Community", myvar)){next}
    print(myvar)
    myplot <- 
      plot_genus_deseq(data_subset,facet_plot = T,
                       limit_to = ".", 
                       title = paste0("community ", community,"; ", myvar), 
                       pallete = my_colour_pallete_12_soft,
                       include_grid = T) 
    
    if(length(unique(with(data_subset, paste0(Group_1,"_", Group_2)))) == 3){
      ggsave(filename = paste0("Result_figures/DESeq_plots/community_",community, "__", myvar, ".pdf"),
             plot = myplot,
             width = 10,
             height = 12,
             units = "cm")   
    } else if(length(unique(with(data_subset, paste0(Group_1,"_", Group_2)))) > 3){
      ggsave(filename = paste0("Result_figures/DESeq_plots/community_",community, "__", myvar, ".pdf"),
             plot = myplot,
             width = 10,
             height = 22,
             units = "cm")   
      
    } else{ # 1
      ggsave(filename = paste0("Result_figures/DESeq_plots/community_",community, "__", myvar, ".pdf"),
             plot = myplot,
             width = 8,
             height = 6,
             units = "cm")    
    }
  }
}


myplot




for (st in unique(otu_dx_groups$Sample_Type)){
  temp <- subset(otu_dx_groups, Sample_Type == st)
  myplot <- plot_genus_deseq(temp,facet_plot = T,
                   limit_to = ".", title = st, pallete = my_colour_pallete_12_soft,include_grid = T)
  if(length(unique(with(temp, paste0(Group_1,"_", Group_2)))) == 3){
    ggsave(filename = paste0("Result_figures/DESeq_plots/paper_UC_CD_CONTROL/", st, "__DX_Groups_DESeq.pdf"),
           plot = myplot,
           width = 10,
           height = 12,
           units = "cm")    
  } else{
    ggsave(filename = paste0("Result_figures/DESeq_plots/paper_UC_CD_CONTROL/", st, "__DX_Groups_DESeq.pdf"),
           plot = myplot,
           width = 12,
           height = 8,
           units = "cm")    
  }

}




# --------------------------------------------------------------------------------

