library(ggtree)
library(ape)
library(ggplot2)
library(phyloseq)
library(ggjoy)
library(dplyr)
library(rstanarm)
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(rstanarm)
library(RColorBrewer)
library(gganimate)


setwd("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS/results/Green_pixels")
dataset_GWAS <- read.csv2("full_dataset_updated_20190415.csv")
dataset_GWAS <- dataset_GWAS[!dataset_GWAS$treatment%in%c("flexible",""),]
metadata <- read.csv2("../../FINAL_Pseudomonas_Phyletic_Pattern_all_geneClusters.Clusterd_by_Jaccard_cd_hit_like.0.99.rep_strains.NonOTU5.csv")
metadata$GWAS_index <- as.character(metadata$GWAS_index)
metadata$GWAS_index[!metadata$GWAS_index%in%c("pathogen_only","control")] <- paste(metadata$GWAS_index[!metadata$GWAS_index%in%c("pathogen_only","control")],"+p",sep = "")
for (treat in unique(dataset_GWAS$treatment)){
  dataset_GWAS$strain[dataset_GWAS$treatment==treat] <- as.character(metadata$strain[metadata$GWAS_index==treat])
}

#re-arrange the pxls data to fit R growth models
dataset_GWAS_nopxls <- dataset_GWAS[,c(1:7,length(dataset_GWAS))]
dataset_GWAS_nopxls_7_no <- cbind(dataset_GWAS_nopxls,"pxls"=dataset_GWAS$X7dpi_nolid)
dataset_GWAS_nopxls_7_no$dpi <- "7_nolid"

dataset_GWAS_nopxls_7 <- cbind(dataset_GWAS_nopxls,"pxls"=dataset_GWAS$X7dpi)
dataset_GWAS_nopxls_7$dpi <- "7"

dataset_GWAS_nopxls_6 <- cbind(dataset_GWAS_nopxls,"pxls"=dataset_GWAS$X6dpi)
dataset_GWAS_nopxls_6$dpi <- "6"

dataset_GWAS_nopxls_5 <- cbind(dataset_GWAS_nopxls,"pxls"=dataset_GWAS$X5dpi)
dataset_GWAS_nopxls_5$dpi <- "5"

dataset_GWAS_nopxls_4 <- cbind(dataset_GWAS_nopxls,"pxls"=dataset_GWAS$X4dpi)
dataset_GWAS_nopxls_4$dpi <- "4"

dataset_GWAS_nopxls_1 <- cbind(dataset_GWAS_nopxls,"pxls"=dataset_GWAS$X1dpi)
dataset_GWAS_nopxls_1$dpi <- "1"

dataset_GWAS_nopxls_0 <- cbind(dataset_GWAS_nopxls,"pxls"=dataset_GWAS$X0dpi)
dataset_GWAS_nopxls_0$dpi <- "0"

dataset_GWAS <- rbind(dataset_GWAS_nopxls_7_no,dataset_GWAS_nopxls_7,dataset_GWAS_nopxls_6,dataset_GWAS_nopxls_5,
                      dataset_GWAS_nopxls_4,dataset_GWAS_nopxls_1,dataset_GWAS_nopxls_0)

# creating 7-0 pxls dataset, to bin by OTU
dataset_GWAS$pxls_7dpi_minus_0dpi[dataset_GWAS$dpi==7] <- dataset_GWAS$pxls[dataset_GWAS$dpi==7]-dataset_GWAS$pxls[dataset_GWAS$dpi==0]
dataset_GWAS_delta_pixels <- dataset_GWAS[dataset_GWAS$dpi==7,]
dataset_GWAS_delta_pixels <- dataset_GWAS_delta_pixels[!dataset_GWAS_delta_pixels$treatment%in% c("boiled_strain_17+p","control"),] # remove boiled protector and control

for (strain in unique(dataset_GWAS_delta_pixels$strain)){
  dataset_GWAS_delta_pixels$median_delta_pxls[dataset_GWAS_delta_pixels$strain==strain] <- (median(dataset_GWAS_delta_pixels$pxls_7dpi_minus_0dpi[dataset_GWAS_delta_pixels$strain==strain]))
}


# read phyletic pattern (a bit long)
phy_pattern <- read.csv("/Volumes/small_projects/hashkenazy/Pseudomonas/PhyleticPatterns/Pseudomonas_Phyletic_Pattern_all_geneClusters.csv")
phy_pattern_my_strains <- phy_pattern[phy_pattern$Species%in%unique(dataset_GWAS$strain[!dataset_GWAS$treatment%in%c("boiled_strain_17+p","control")]),] # subset to my strains only

phy_pattern_my_strains$median_delta_pxls <- dataset_GWAS_delta_pixels$median_delta_pxls[match(phy_pattern_my_strains$Species, dataset_GWAS_delta_pixels$strain)] # add the phenotype (median pxls 7dpi minus 0 dpi) to phyletic pattern variable

#### chose OG of interest
#subset_siginifact_OGs <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/OTU2_significant_OGs.csv", header = F) #subset to best OTU2 hits (signal from all strains)
subset_siginifact_OGs <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/OTU2_significant_OGs_updated.csv", header = F) #subset to best OTU2 hits (signal from all strains)
#subset_siginifact_OGs <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/preliminary_data/GWAS_allstrains_hit_list.csv", header = F) #subset to best OTU2 hits (signal from all strains)



# Add OTU for subsequent subsetting (e.g. heatmap of only OTU2 strains)
dataset_GWAS_delta_pixels$OTU <- metadata$O_T_U[match(dataset_GWAS_delta_pixels$strain,metadata$strain)]
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)

dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="control"] <- "Bacteria free" #adding "OTU" level to control for subsequent plotting
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$OTU==""] <- "OTU5" # adding OTU5






#####heatmapping part
#strains_by_OTU <- unique(dataset_GWAS_delta_pixels$strain[dataset_GWAS_delta_pixels$OTU %in% "OTU2"])
strains_by_OTU <- unique(dataset_GWAS_delta_pixels$strain)

strains_by_OTU <- strains_by_OTU[!is.na(strains_by_OTU)]
phy_pattern_my_strains_by_OTU <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_by_OTU,]

OTU <- metadata$O_T_U[metadata$strain %in% strains_by_OTU & !(metadata$GWAS_index%in%c("control","boiled_strain_17+p"))]
strain <- metadata$strain[metadata$strain %in% strains_by_OTU & !(metadata$GWAS_index%in%c("control","boiled_strain_17+p"))]
OTU <- as.character(OTU)
OTU[OTU==""] <- "OTU5"
OTU_strain <- data.frame("strain"=strain,"OTU"=OTU)


best_hits_phy_pattern <- phy_pattern_my_strains_by_OTU[,as.character(subset_siginifact_OGs$V1)]
row.names(best_hits_phy_pattern) <- phy_pattern_my_strains_by_OTU$Species

##### stats of strains with and without the significant hits
strains_with_sig_hits <- rownames(best_hits_phy_pattern[rowSums(best_hits_phy_pattern) > 0,])
strains_without_sig_hits <- rownames(best_hits_phy_pattern[rowSums(best_hits_phy_pattern) == 0,])
sig_hits_pheno_df <- phy_pattern_my_strains_by_OTU[,c("Species","median_delta_pxls")]
sig_hits_pheno_df$group[sig_hits_pheno_df$Species %in% strains_with_sig_hits] <- "significant_strain"
sig_hits_pheno_df$group[sig_hits_pheno_df$Species %in% strains_without_sig_hits] <- "not_significant_strain"
sig_hits_pheno_df$Species
sig_hits_pheno_df$sig_hits_num <- rowSums(best_hits_phy_pattern)
sig_hits_pheno_df_final <- merge(x = sig_hits_pheno_df, y = OTU_strain, by.x = "Species", by.y = "strain")

sig_hits_pheno_df_final$group_and_OTU[sig_hits_pheno_df_final$Species %in% strains_with_sig_hits] <- "sig. orthologs>1"
sig_hits_pheno_df_final$group_and_OTU[(sig_hits_pheno_df_final$Species %in% strains_without_sig_hits) & (sig_hits_pheno_df_final$OTU=="OTU2")] <- "sig. orthologs<1"
sig_hits_pheno_df_final$group_and_OTU[(sig_hits_pheno_df_final$Species %in% strains_without_sig_hits) & (sig_hits_pheno_df_final$OTU!="OTU2")] <- "sig. orthologs<1\nnon-ATUE2"

sig_hits_pheno_df_final$group_and_OTU <- factor(sig_hits_pheno_df_final$group_and_OTU, levels = c("sig. orthologs<1\nnon-ATUE2", "sig. orthologs<1", "sig. orthologs>1")) 



ggplot(data = sig_hits_pheno_df_final, aes(x = group_and_OTU, y=median_delta_pxls))+
  geom_boxplot() +
  geom_point() +
  theme_classic()

ggplot(data = sig_hits_pheno_df_final[sig_hits_pheno_df_final$OTU=="OTU2",], aes(x = group, y=median_delta_pxls))+
  geom_boxplot() +
  geom_point()

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/lm_number_of_sig_hits_OTU2_only.pdf", useDingbats = F)
ggplot(data = sig_hits_pheno_df_final[sig_hits_pheno_df_final$OTU=="OTU2",], aes(x = sig_hits_num, y=median_delta_pxls))+
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9))
#! dev.off()

fit_sig_hits <- stan_glm(formula = median_delta_pxls ~ group_and_OTU, data = sig_hits_pheno_df_final)


grid = sig_hits_pheno_df_final %>%
  data_grid(group_and_OTU)

fits = grid %>%
  add_fitted_draws(fit_sig_hits)

preds = grid %>%
  add_predicted_draws(fit_sig_hits)



#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/sig_hits_strains_effect.pdf", useDingbats = F)
sig_hits_pheno_df_final %>%
  ggplot(aes(y = group_and_OTU, x = median_delta_pxls)) +
  stat_interval(aes(x = .prediction), data = preds) +
  stat_pointinterval(aes(x = .value), data = fits, .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()

quantile(x = fits$.value[fits$group_and_OTU=="sig. orthologs<1\nnon-ATUE2"],probs = c(0.025,0.5,0.975))
quantile(x = fits$.value[fits$group_and_OTU=="sig. orthologs<1"],probs = c(0.025,0.5,0.975))
quantile(x = fits$.value[fits$group_and_OTU=="sig. orthologs>1"],probs = c(0.025,0.5,0.975))



cor.test(sig_hits_pheno_df_final$median_delta_pxls[sig_hits_pheno_df_final$OTU=="OTU2"],sig_hits_pheno_df_final$sig_hits_num[sig_hits_pheno_df_final$OTU=="OTU2"])
plot(sig_hits_pheno_df_final$median_delta_pxls[sig_hits_pheno_df_final$OTU=="OTU2"],sig_hits_pheno_df_final$sig_hits_num[sig_hits_pheno_df_final$OTU=="OTU2"])

fit_sig_hits_within_OTU2 <- stan_glm(formula = median_delta_pxls ~ sig_hits_num, data = sig_hits_pheno_df_final[sig_hits_pheno_df_final$OTU=="OTU2",])
fit_sig_hits_ <- stan_glm(formula = median_delta_pxls ~ sig_hits_num, data = sig_hits_pheno_df_final)


fit_sig_hits$stan_summary
fit_sig_hits_within_OTU2$stan_summary

plot(fit_sig_hits)
plot(fit_sig_hits_within_OTU2)

#create a phenotype df
phenotype_continious <- "median_delta_pxls"
my_pheno <- data.frame("nodes"=phy_pattern_my_strains_by_OTU$Species,"phenotype"=phy_pattern_my_strains_by_OTU[,phenotype_continious])

tree_phyl_all_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_my_strains_aligned_no_p4.C9.fasta.treefile"
tree_phyl_OTU2_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_OTU2_strains_aligned.fasta.treefile"
OTU_of_interest="all"

if (OTU_of_interest=="all"){
  tree=  tree_phyl_all_strains
  OTU <- paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
} else if (OTU_of_interest=="OTU2"){
  tree=  tree_phyl_OTU2_strains
  OTU <- "OTU2"
}

cls <- list(other_strains=OTU_strain$strain[OTU_strain$OTU!="OTU2"], OTU2=OTU_strain$strain[OTU_strain$OTU=="OTU2"])

tree_phyl <- read.tree(tree)
tree_phyl_OTU <- groupOTU(tree_phyl, cls)

p<- ggtree(tree_phyl_OTU)

p <- p %<+% my_pheno

library("colorspace")

#p <- ggtree(tr = tree_phyl_OTU)
  #+geom_tippoint(aes(color=group))

best_hits_phy_pattern_phenotype <- merge(x = best_hits_phy_pattern, y = my_pheno, by.x=0, by.y = "nodes")
rownames(best_hits_phy_pattern_phenotype) <- best_hits_phy_pattern_phenotype$Row.names
best_hits_phy_pattern_phenotype <- best_hits_phy_pattern_phenotype[,!(names(best_hits_phy_pattern_phenotype) %in% "Row.names")]
best_hits_phy_pattern_phenotype$phenotype <- best_hits_phy_pattern_phenotype$phenotype/max(best_hits_phy_pattern_phenotype$phenotype)
best_hits_phy_pattern_phenotype$phenotype[best_hits_phy_pattern_phenotype$phenotype<0] <- 0
best_hits_phy_pattern_phenotype$phenotype[best_hits_phy_pattern_phenotype$phenotype==1] <- 0.99



p_h=gheatmap(p, best_hits_phy_pattern_phenotype, offset=0.05, width =0.8, font.size=10, 
             colnames_angle=90, hjust=1,low = "white", high = "red",color = "black",colnames = T) 



color_OTU <- OTU_strain$OTU!="OTU2"
color_OTU[color_OTU==T] <- "red"
color_OTU[color_OTU==F] <- "black"

#tt <- p_h %>% scale_x_ggtree() + geom_tiplab(size=2, align=TRUE, linesize=.5, aes(color=group))
tt <- p_h  + geom_tiplab(size=2, align=TRUE, linesize=.5, aes(color=group))

tt <- tt + scale_color_manual(values = c("black","#6300ff")) + scale_fill_gradient(low = "white",high = "darkgreen") # add differential color for OTU2

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/heatmap_sig_OG_all_strains_and_phenotype.pdf", width = 14, height = 15,useDingbats = F)
#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/heatmap_ALL_sig_OG_all_strains_and_phenotype.pdf", width = 14, height = 15,useDingbats = F)
print(tt)
#! dev.off()


##### now producing the same, but coloring ATUE2, ATUE3, and ATUE4 seperately, 
##### and presenting all significant orthologs with SC >0

subset_siginifact_OGs <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/TableS2/positive_sig_OGs.csv", header = T) #subset to best OTU2 hits (signal from all strains)


# Add OTU for subsequent subsetting (e.g. heatmap of only OTU2 strains)
dataset_GWAS_delta_pixels$OTU <- metadata$O_T_U[match(dataset_GWAS_delta_pixels$strain,metadata$strain)]
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)

dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="control"] <- "Bacteria free" #adding "OTU" level to control for subsequent plotting
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$OTU==""] <- "OTU5" # adding OTU5






#####heatmapping part
#strains_by_OTU <- unique(dataset_GWAS_delta_pixels$strain[dataset_GWAS_delta_pixels$OTU %in% "OTU2"])
strains_by_OTU <- unique(dataset_GWAS_delta_pixels$strain)

strains_by_OTU <- strains_by_OTU[!is.na(strains_by_OTU)]
phy_pattern_my_strains_by_OTU <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_by_OTU,]

OTU <- metadata$O_T_U[metadata$strain %in% strains_by_OTU & !(metadata$GWAS_index%in%c("control","boiled_strain_17+p"))]
strain <- metadata$strain[metadata$strain %in% strains_by_OTU & !(metadata$GWAS_index%in%c("control","boiled_strain_17+p"))]
OTU <- as.character(OTU)
OTU[OTU==""] <- "OTU5"
OTU_strain <- data.frame("strain"=strain,"OTU"=OTU)


best_hits_phy_pattern <- phy_pattern_my_strains_by_OTU[,as.character(subset_siginifact_OGs$OG)]
row.names(best_hits_phy_pattern) <- phy_pattern_my_strains_by_OTU$Species


cls <- list(other_strains=OTU_strain$strain[!(OTU_strain$OTU %in% c("OTU2", "OTU3", "OTU4"))], OTU2=OTU_strain$strain[OTU_strain$OTU=="OTU2"], 
            OTU3=OTU_strain$strain[OTU_strain$OTU=="OTU3"], OTU4=OTU_strain$strain[OTU_strain$OTU=="OTU4"])

tree_phyl <- read.tree(tree)
tree_phyl_OTU <- groupOTU(tree_phyl, cls)

p<- ggtree(tree_phyl_OTU)

p <- p %<+% my_pheno

library("colorspace")

#p <- ggtree(tr = tree_phyl_OTU)
#+geom_tippoint(aes(color=group))

best_hits_phy_pattern_phenotype <- merge(x = best_hits_phy_pattern, y = my_pheno, by.x=0, by.y = "nodes")
rownames(best_hits_phy_pattern_phenotype) <- best_hits_phy_pattern_phenotype$Row.names
best_hits_phy_pattern_phenotype <- best_hits_phy_pattern_phenotype[,!(names(best_hits_phy_pattern_phenotype) %in% "Row.names")]
best_hits_phy_pattern_phenotype$phenotype <- best_hits_phy_pattern_phenotype$phenotype/max(best_hits_phy_pattern_phenotype$phenotype)
best_hits_phy_pattern_phenotype$phenotype[best_hits_phy_pattern_phenotype$phenotype<0] <- 0
best_hits_phy_pattern_phenotype$phenotype[best_hits_phy_pattern_phenotype$phenotype==1] <- 0.99



p_h=gheatmap(p, best_hits_phy_pattern_phenotype, offset=0.05, width =0.8, font.size=10, 
             colnames_angle=90, hjust=1,low = "white", high = "red",color = "black",colnames = T) 



color_OTU <- OTU_strain$OTU=="OTU2"
color_OTU[color_OTU==T] <- "#6300ff"
color_OTU[color_OTU==F] <- "black"
color_OTU[OTU_strain$OTU=="OTU3"] <- "#b2df8a"
color_OTU[OTU_strain$OTU=="OTU4"] <- "#1f78b4"


#tt <- p_h %>% scale_x_ggtree() + geom_tiplab(size=2, align=TRUE, linesize=.5, aes(color=group))
tt <- p_h  + geom_tiplab(size=2, align=TRUE, linesize=.5, aes(color=group))

tt <- tt + scale_color_manual(values = c("black","#6300ff","#b2df8a","#1f78b4")) + scale_fill_gradient(low = "white",high = "darkgreen") # add differential color for OTU2

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/heatmap_sig_OG_all_strains_and_phenotype.pdf", width = 14, height = 15,useDingbats = F)
#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/heatmap_ALL_sig_OG_all_strains_and_phenotype.pdf", width = 14, height = 15,useDingbats = F)
print(tt)
#! dev.off()




###now the same, focusing only on OTU2 group


OTU_of_interest="OTU2"

if (OTU_of_interest=="all"){
  tree=  tree_phyl_all_strains
  OTU <- paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
} else if (OTU_of_interest=="OTU2"){
  tree=  tree_phyl_OTU2_strains
  OTU <- "OTU2"
}

tree_phyl <- read.tree(tree)

p<- ggtree(tree_phyl)

p <- p %<+% my_pheno


p <- ggtree(tr = tree_phyl)
#+geom_tippoint(aes(color=group))

best_hits_phy_pattern_phenotype <- merge(x = best_hits_phy_pattern, y = my_pheno, by.x=0, by.y = "nodes")
rownames(best_hits_phy_pattern_phenotype) <- best_hits_phy_pattern_phenotype$Row.names
best_hits_phy_pattern_phenotype <- best_hits_phy_pattern_phenotype[,!(names(best_hits_phy_pattern_phenotype) %in% "Row.names")]
best_hits_phy_pattern_phenotype$phenotype <- best_hits_phy_pattern_phenotype$phenotype/max(best_hits_phy_pattern_phenotype$phenotype)
best_hits_phy_pattern_phenotype$phenotype[best_hits_phy_pattern_phenotype$phenotype<0] <- 0
best_hits_phy_pattern_phenotype$phenotype[best_hits_phy_pattern_phenotype$phenotype==1] <- 0.99



p_h=gheatmap(p, best_hits_phy_pattern_phenotype, offset=0.003, width =0.8, font.size=7, 
             colnames_angle=90, hjust=1,low = "white", high = "red",color = "black",colnames = T) 


tt <- p_h %>% scale_x_ggtree() + geom_tiplab(size=2, align=TRUE, linesize=.5)
tt <- tt + scale_fill_gradient(low = "white",high = "darkgreen")

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/heatmap_sig_OG_OTU2_strains_and_phenotype.pdf", width = 15, height = 15, useDingbats = F)
print(tt)
#! dev.off()



# with barplot
#phenotype <- my_pheno$phenotype[match(p$data$label,my_pheno$nodes)]
#d2 <- data.frame(id=p$data$label, value = phenotype)
#d2$value[is.na(d2$id)] <- NA
#facet_plot(tt, panel='bar', data=d2, geom=geom_segment, aes(x=0, xend=value, y=y, yend=y), size=3, color='steelblue') + theme_tree2()


######growth scatter plot by number of significant OGs




# now lets scatter plot this, using mean or median
MinMeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

#now lets create the dataframe that can show mean/median + error bars
scatter_weight_df <- data.frame(t(rep(NA,5)))
full_dataset <- dataset_GWAS
for (dpi in unique(full_dataset$dpi)){
  for (treat in unique(full_dataset$treatment)){
    mean_per_treat_per_gen <- full_dataset$pxls[full_dataset$dpi==dpi & full_dataset$treatment==treat]
    mean_pxls <- median(mean_per_treat_per_gen[!is.na(mean_per_treat_per_gen)])
    sd_pxls <- sd(mean_per_treat_per_gen[!is.na(mean_per_treat_per_gen)])
    sem_pxls <- sd_pxls/sqrt(length(mean_per_treat_per_gen[!is.na(mean_per_treat_per_gen)]))
    scatter_weight_df <- rbind(scatter_weight_df,c(mean_pxls,sd_pxls,sem_pxls,dpi,treat))
  }
}
colnames(scatter_weight_df) <- c("pxls_average","sd","sem","dpi","treatment")
scatter_weight_df$pxls_average <- as.numeric(scatter_weight_df$pxls_average)
scatter_weight_df$sd <- as.numeric(scatter_weight_df$sd)
scatter_weight_df$sem <- as.numeric(scatter_weight_df$sem )
scatter_weight_df$dpi <- as.factor(scatter_weight_df$dpi)
scatter_weight_df$treatment <- as.factor(scatter_weight_df$treatment)
scatter_weight_df_subset <- scatter_weight_df[scatter_weight_df$dpi%in%c("0","1","4","5","7","6") ,]

scatter_weight_df_subset$strain <- metadata$strain[match(scatter_weight_df_subset$treatment,metadata$GWAS_index)]
scatter_weight_df_subset <- scatter_weight_df_subset[!scatter_weight_df_subset$treatment %in% c("boiled_strain_17+p","pathogen_only", "control"),]


#strating with OTU2
chosen_OTU <- "OTU2"
phy_pattern_best_hits <- phy_pattern_my_strains[,c(as.character(subset_siginifact_OGs[[1]]),"Species")]
rownames(phy_pattern_my_strains) <- phy_pattern_my_strains$Species
phy_best_hits_sums <- data.frame("strain"=phy_pattern_best_hits$Species,"number_of_genes"= as.numeric(rowSums(phy_pattern_best_hits[,subset_siginifact_OGs[[1]]])))

scatter_weight_df_subset$number_of_genes <- phy_best_hits_sums$number_of_genes[match(scatter_weight_df_subset$strain,phy_best_hits_sums$strain)] 
scatter_weight_df_subset$number_of_genes[is.na(scatter_weight_df_subset$number_of_genes)] <- 0

scatter_weight_df_subset$OTU <- metadata$O_T_U[match(scatter_weight_df_subset$strain,metadata$strain)]

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/OTU2_sig_OG_growth_scatter_plot.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes), 
       data = scatter_weight_df_subset[scatter_weight_df_subset$OTU==chosen_OTU,]) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  scale_colour_gradient(low = "grey87", high = "#756bb1")+
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()


#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/OTU2_sig_OG_growth_scatter_plot_all_strains.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes), 
       data = scatter_weight_df_subset) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  scale_colour_gradient(low = "grey87", high = "#756bb1")+
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()




