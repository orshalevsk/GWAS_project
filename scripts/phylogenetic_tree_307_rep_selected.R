library(ggtree)
library(ape)
library(ggplot2)
library(phyloseq)
library(ggjoy)
library(dplyr)
library(phytools)

# read GWAS data
setwd("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS/results/Green_pixels")

dataset_GWAS <- read.csv2("full_dataset_updated_20190415.csv")
dataset_GWAS <- dataset_GWAS[!dataset_GWAS$treatment%in%c("flexible",""),]
metadata <- read.csv2("../../FINAL_Pseudomonas_Phyletic_Pattern_all_geneClusters.Clusterd_by_Jaccard_cd_hit_like.0.99.rep_strains.NonOTU5.csv")
metadata$GWAS_index <- as.character(metadata$GWAS_index)
metadata$GWAS_index[!metadata$GWAS_index%in%c("pathogen_only","control")] <- paste(metadata$GWAS_index[!metadata$GWAS_index%in%c("pathogen_only","control")],"+p",sep = "")
for (treat in unique(dataset_GWAS$treatment)){
  dataset_GWAS$strain[dataset_GWAS$treatment==treat] <- as.character(metadata$strain[metadata$GWAS_index==treat])
}




#read general phylogeny data about 307 rep strains
metadata <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/all_Pseudomonas_genomes_meta_data.csv")

metadata_307 <- metadata[metadata$species %in% tree_307_rep$tip.label,]


tree_307_rep <- read.tree("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Pseudomonas_Phyletic_Pattern_all_geneClusters.Clusterd_by_Jaccard_cd_hit_like.0.99.rep_strains.Core_OG_frac0.9.MSA.AA.LG.FastTree.tre")

OTU5_remove <- as.character(metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU5"])
metadata_307 <- metadata_307[!metadata_307$species %in% OTU5_remove,]

write.csv(x = metadata_307,file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/metadata_127_commensals.csv")


tree_307_rep <- drop.tip(phy = tree_307_rep, tip = OTU5_remove )
#! write.tree(phy = tree_307_rep, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/127_rep_commensals.tree")

cls <- list(ATUE2=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU2"],
            ATUE3=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU3"],
            ATUE4=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU4"],
            #ATUE5=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU5"],
            ATUE6=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU6"],
            ATUE7=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU7"],
            ATUE8=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU8"],
            ATUE9=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU9"],
            ATUE10=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU10"],
            ATUE11=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU11"],
            ATUE_NEW=metadata_307$species[metadata_307$OTU_assignment_8_2017=="OTU_NEW"])

tree_gg <- ggtree(tree_307_rep,layout = "circular")
tree_gg_OTU <- groupOTU(tree_gg, cls)
tree_gg_OTU <- tree_gg_OTU  %<+% metadata_307
library("colorspace")
tree_gg_OTU$data$OTU_assignment_8_2017 <- as.character(tree_gg_OTU$data$OTU_assignment_8_2017)
tree_gg_OTU$data$OTU_assignment_8_2017[is.na(tree_gg_OTU$data$OTU_assignment_8_2017)] <- 0
tree_gg_OTU$data$OTU_assignment_8_2017 <- factor(tree_gg_OTU$data$OTU_assignment_8_2017, levels = c("0","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10","OTU11","OTU_NEW"))



# add grouping by chosen strains
strains <- unique(dataset_GWAS$strain)
strains <- strains[!is.na(strains)]
strains_df <- data.frame("strain"=strains, "treatment"="tested_strain")
strains_df$treatment <- as.character(strains_df$treatment)
strains_df$treatment[strains_df$strain=="p4.C9"] <- "pathogen"
strains_df$treatment[strains_df$strain=="p5.F2"] <- "protective strain"

# add the subset of strains i used to the tree dataframe
tree_gg_OTU$data$chosen_strains<- strains_df$treatment[match(tree_gg_OTU$data$label, strains_df$strain)]

colors <- c("#6300ff",	"#ff00f2",	"#00fff1","#ffff00",	"#0080ff",	"#685e5e",	"#d33b00",	"#ff8000", "#ff0000",	"#4a8220",	"#3a8403", "white", "darkred")
#colors=c("white","darkred")


tree_gg_OTU <- scaleClade(tree_gg_OTU, node=352, scale=4) 
tree_gg_OTU  %>% scale_x_ggtree() + 
  geom_tiplab(size=2, align=TRUE, linesize=.2, aes(color=OTU_assignment_8_2017),offset = 0.03,hjust = 1,geom = 'text',angle=90) + 
  geom_tippoint(aes(color=chosen_strains),shape=15, size=1.5) + scale_color_manual(values = colors ) 

OTU_307 <- data.frame(metadata_307$OTU_assignment_8_2017)
rownames(OTU_307) <- metadata_307$species

p_h=gheatmap(tree_gg_OTU, OTU_307, offset=0.003, width =0.8, font.size=7, 
             colnames_angle=90, hjust=1,color = "black",colnames = T) 


#tree_gg_OTU + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

