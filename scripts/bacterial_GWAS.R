'''
This one is to for bacterial GWAs, done with 103 strains vs. one pathogen (p4.C9). 
will use different phenotypes: auc_l, median_day7-day0, 50% of credible intervals, stan model
'''
library(ggplot2)
library(randomForest)
library(growthcurver)
library(rstanarm)
library(treeWAS)
library(gplots)
library(tiger)
library(RColorBrewer)

#loading data
source("~/ownCloud/My papers/bacterial_GWAS_paper/scripts/Database_return.R")
dataset_GWAS <-database_return()
levels(dataset_GWAS$O_T_U)[1] <- "OTU5"

##scatter plot the data

source("~/ownCloud/My papers/bacterial_GWAS_paper/scripts/data_to_scatter_plot.R")
scatter_weight_df <- scatter_plot_data(dataset_GWAS)
scatter_weight_df <- scatter_weight_df[!is.na(scatter_weight_df$treatment),]
#adding controls in OTU, to avoid NA
scatter_weight_df$OTU <- as.character(scatter_weight_df$OTU)
scatter_weight_df$OTU[scatter_weight_df$treatment %in% c("control","boiled_strain_17+p","pathogen_only")] <- "control"
scatter_weight_df$OTU <- as.factor(scatter_weight_df$OTU)

# coloring by controls
colors_scatter_by_controls <- rep("gray87",length(scatter_weight_df$treatment))
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df$treatment))=="control",6)] <- "black"
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df$treatment))=="boiled_strain_17+p",6)] <- "orange"
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df$treatment))=="pathogen_only",6)] <- "red"
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df$treatment))=="30+p",6)] <- "green"

#plotting general results, colors by control / treatment

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS_few_pathogens_by_protective/figures/1st_rep_p4.C9.pdf")
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), 
            data = scatter_weight_df[scatter_weight_df$treatment%in%c("pathogen_only","control","30+p","boiled_strain_17+p"),]) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), 
             data = scatter_weight_df[scatter_weight_df$treatment%in%c("pathogen_only","control","30+p","boiled_strain_17+p"),])+
  scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()




#using the "growthcurver package" - calculating k - capacity,r - growth rate, and auc_l - taking both r and k
growth_results <- data.frame(t(c(NA,NA,NA,NA)))
colnames(growth_results) <- c("treatment","r","k","auc_l")
for (treat in unique(dataset_GWAS$treatment)){
  gc_fit <- SummarizeGrowth(as.numeric(dataset_GWAS$dpi[dataset_GWAS$treatment==treat]), dataset_GWAS$pxls[dataset_GWAS$treatment==treat])
  growth_results <- rbind(growth_results,c(treat,gc_fit$vals[[7]],gc_fit$vals[[1]],gc_fit$vals[[14]]))
}

growth_results <- growth_results[!is.na(growth_results$r),]
growth_results <- growth_results[order(growth_results$k,decreasing = T),]

growth_results$k <- as.numeric(growth_results$k)
growth_results$r <- as.numeric(growth_results$r)
growth_results$auc_l <- as.numeric(growth_results$auc_l)

growth_results$strain <- dataset_GWAS$strain[match(growth_results$treatment,dataset_GWAS$treatment)]
growth_results$strain <- factor(growth_results$strain,levels = growth_results$strain[order(growth_results$auc_l,decreasing = T)])

#sanity check - how does it looks like, ordered
ggplot(aes(x =strain,y=auc_l),data = growth_results) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# merging the growth curver data with the original data
dataset_GWAS <- merge(dataset_GWAS, growth_results[,c("strain","r","k","auc_l")], by = "strain")

# creating now median 7days-0days phenotype
for (treat in unique(dataset_GWAS$treatment)){
  dpi7<- dataset_GWAS$pxls[dataset_GWAS$treatment==treat & dataset_GWAS$dpi==7]
  dpi0 <- dataset_GWAS$pxls[dataset_GWAS$treatment==treat & dataset_GWAS$dpi==0]
  dpi7_0 <- dpi7-dpi0
  dataset_GWAS$dpi7_0median[dataset_GWAS$treatment==treat & dataset_GWAS$dpi==7] <- median(dpi7_0)
  dataset_GWAS$dpi7_0[dataset_GWAS$treatment==treat & dataset_GWAS$dpi==7] <- dpi7_0
}

dataset_GWAS$strain <- factor(dataset_GWAS$strain, levels = unique(dataset_GWAS$strain[order(dataset_GWAS$dpi7_0median,decreasing = T)]))

#sanity check - how does it looks like, ordered
ggplot(aes(x = strain, y=dpi7_0), data = dataset_GWAS[dataset_GWAS$dpi==7,])+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#now i will add the last phenotype - the stan model, credible intervals, for 7 dpi - 0 dpi, per strain.
strain_pxls_7_0_stan_mod <- stan_glm(dpi7_0 ~ strain, data = dataset_GWAS, iter = 8000, seed=12345)
strain_pxls_7_0_stan_mod_df <- as.data.frame(strain_pxls_7_0_stan_mod$stan_summary)
strain_pxls_7_0_stan_mod_df$strain <- rownames(strain_pxls_7_0_stan_mod_df)
strain_pxls_7_0_stan_mod_df <- strain_pxls_7_0_stan_mod_df[!strain_pxls_7_0_stan_mod_df$strain %in% c("sigma","mean_PPD","log-posterior"),]
strain_pxls_7_0_stan_mod_df$strain[strain_pxls_7_0_stan_mod_df$strain=="(Intercept)"] <- paste("strain",levels(dataset_GWAS$strain)[1],sep = "")
strain_pxls_7_0_stan_mod_df$strain <- vapply(strsplit(strain_pxls_7_0_stan_mod_df$strain, split = "strain"), `[`, 2, FUN.VALUE=character(1))

strain_pxls_7_0_stan_mod_df$strain <- factor(strain_pxls_7_0_stan_mod_df$strain, levels = strain_pxls_7_0_stan_mod_df$strain[order(strain_pxls_7_0_stan_mod_df$`50%`,decreasing = T)])

# sanity check again. How does it looks - ordered
ggplot(aes(x = strain, y=`50%`, ymin = `2.5%`, ymax = `97.5%`), data = strain_pxls_7_0_stan_mod_df)+
  geom_pointrange() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# merging the stan model results with the original data
dataset_GWAS <- merge(dataset_GWAS, strain_pxls_7_0_stan_mod_df, by="strain")

#agreement between the ranking (ordered) of the three different phenotypes.
match(levels(growth_results$strain),levels(dataset_GWAS$strain))
match(levels(strain_pxls_7_0_stan_mod_df$strain),levels(dataset_GWAS$strain))


#Now starts the GWAS part
#load phyletic pattern and filter genes with at least 2 occurences (filter unique and 0)
phy_pattern <- read.csv("/Volumes/small_projects/hashkenazy/Pseudomonas/PhyleticPatterns/Pseudomonas_Phyletic_Pattern_all_geneClusters.csv")
phy_pattern_my_strains <- phy_pattern[phy_pattern$Species%in%unique(dataset_GWAS$strain[!dataset_GWAS$treatment%in%c("boiled_strain_17+p","control")]),]

#merging the phenotypes, with the phyletic pattern
dataset_GWAS_7dpi <- dataset_GWAS[dataset_GWAS$dpi==7,]
phy_pattern_my_strains$cdl50 <- dataset_GWAS_7dpi$`50%`[match(phy_pattern_my_strains$Species,dataset_GWAS_7dpi$strain)]
phy_pattern_my_strains$dpi7_0median <- dataset_GWAS_7dpi$dpi7_0median[match(phy_pattern_my_strains$Species,dataset_GWAS_7dpi$strain)]
phy_pattern_my_strains$auc_l <- dataset_GWAS_7dpi$auc_l[match(phy_pattern_my_strains$Species,dataset_GWAS_7dpi$strain)]

#adding one more phenotype - cdl50, making it binray: while everything below 0 pxls, get the same value - 0, and everything above = 1
phy_pattern_my_strains$cdl50_binary[phy_pattern_my_strains$cdl50>0] <- 1
phy_pattern_my_strains$cdl50_binary[phy_pattern_my_strains$cdl50<0] <- 0


#limiting only to OGs with at least 2 occurences, and with at least one strain who miss it
number_of_occurences_per_gene<- colSums(phy_pattern_my_strains[,!names(phy_pattern_my_strains) %in% c("cdl50","auc_l","dpi7_0median","Species")])
at_least=names(number_of_occurences_per_gene[number_of_occurences_per_gene>1 & number_of_occurences_per_gene<length(phy_pattern_my_strains$Species)])
phy_pattern_my_strains <- phy_pattern_my_strains[,c("Species",at_least,"cdl50","auc_l","dpi7_0median")]
phy_pattern_my_strains_OGs <- phy_pattern_my_strains[,!names(phy_pattern_my_strains) %in% c("cdl50","auc_l","dpi7_0median","Species")]



#first - spearman corr, for all strains, cdl50 as phenotype
cor_spearman_cdl50 <- cor(y = phy_pattern_my_strains$cdl50, 
                                x= phy_pattern_my_strains[,grepl(x = names(phy_pattern_my_strains),pattern = "GC")],
                                method = "spearman")
cor_spearman_cdl50 <- cor_spearman_cdl50[order(cor_spearman_cdl50, decreasing = T), , drop = FALSE] # sort by correlation
cor_spearman_cdl50 <- data.frame("OG"=row.names(cor_spearman_cdl50),"cor_all"=cor_spearman_cdl50)


##first - spearman corr, OTU-specific, cdl50 as phenotype
#OTU2
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(2),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]

cor_spearman_cdl50_OTU2 <- cor(y = phy_pattern_my_strains_subset$cdl50, 
                          x= phy_pattern_my_strains_subset[,grepl(x = names(phy_pattern_my_strains_subset),pattern = "GC")],
                          method = "spearman")
cor_spearman_cdl50_OTU2 <- cor_spearman_cdl50_OTU2[order(cor_spearman_cdl50_OTU2, decreasing = T), , drop = FALSE] # sort by correlation
cor_spearman_cdl50_OTU2 <- data.frame("OG"=row.names(cor_spearman_cdl50_OTU2),"cor_OTU2"=cor_spearman_cdl50_OTU2)

#OTU3
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(3),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]

cor_spearman_cdl50_OTU3 <- cor(y = phy_pattern_my_strains_subset$cdl50, 
                               x= phy_pattern_my_strains_subset[,grepl(x = names(phy_pattern_my_strains_subset),pattern = "GC")],
                               method = "spearman")
cor_spearman_cdl50_OTU3 <- cor_spearman_cdl50_OTU3[order(cor_spearman_cdl50_OTU3, decreasing = T), , drop = FALSE] # sort by correlation
cor_spearman_cdl50_OTU3 <- data.frame("OG"=row.names(cor_spearman_cdl50_OTU3),"cor_OTU3"=cor_spearman_cdl50_OTU3)

#OTU4
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(4),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]

cor_spearman_cdl50_OTU4 <- cor(y = phy_pattern_my_strains_subset$cdl50, 
                               x= phy_pattern_my_strains_subset[,grepl(x = names(phy_pattern_my_strains_subset),pattern = "GC")],
                               method = "spearman")
cor_spearman_cdl50_OTU4 <- cor_spearman_cdl50_OTU4[order(cor_spearman_cdl50_OTU4, decreasing = T), , drop = FALSE] # sort by correlation
cor_spearman_cdl50_OTU4 <- data.frame("OG"=row.names(cor_spearman_cdl50_OTU4),"cor_OTU4"=cor_spearman_cdl50_OTU4)

cor_spearman_cdl50_merged <- merge(merge(merge(cor_spearman_cdl50,cor_spearman_cdl50_OTU2, by="OG"),cor_spearman_cdl50_OTU3, by="OG"), cor_spearman_cdl50_OTU4, by="OG")
#now i will rename the "OG" column to "Var1". This part will be just thechnically important for later - to merge with treeWAS dataset.
colnames(cor_spearman_cdl50_merged)[colnames(cor_spearman_cdl50_merged)=="OG"] <- "Var1"



####################treeWAS part
##now i will run treeWAS for different subset of strains (OTU-wise), or with different phenotyoe / correction parameters

#first -loading the small function i wrote for running treeWAS
source("~/ownCloud/My papers/bacterial_GWAS_paper/scripts/GWAS_for_OTU_function.R")

#tree_path. Should be the same for all OTUs (later i will prune the tree, in the  function itself)
tree_phyl_all_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_my_strains_aligned.fasta.treefile"


###all strains treeWAS
#strating with all strains, phenotype = cdl50 (median of credible intervals, after the stan model)
#subset to the strains/OTU if needed
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl_50_all_strains_2nd_replicate.pdf"
#run the treeWAS function
cdl50_all_strains <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50", output_file = outputfile)

#now cdl50 binary as phenotype
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl50_binary_all_strains_2nd_replicate.pdf"
#run the treeWAS function
cdl50_binary_all_strains <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50_binary", output_file = outputfile)



#all strains, phenotype = median 7dpi-0dpi (median of credible intervals, after the stan model)
#subset to the strains/OTU if needed
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_7dpi-0dpi/treewas_dpi7_0median_all_strains_2nd_replicate.pdf"
#run the treeWAS function
dpi7_all_strains <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "dpi7_0median", output_file = outputfile)


#all strains, phenotype = auc_l (growthcurver )
#subset to the strains/OTU if needed
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_auc_l/treewas_auc_l_all_strains_2nd_replicate.pdf"
#run the treeWAS function
auc_l_all_strains <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "auc_l", output_file = outputfile)



###OTU4 treeWAS
#phenotype = cdl50 (median of credible intervals, after the stan model)
#subset to OTU4
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(4),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl_50_OTU4_2nd_replicate.pdf"
#run the treeWAS function
cdl50_OTU4 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50", output_file = outputfile)


#now cdl50 binary as phenotype
#subset to OTU4
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(4),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl50_binary_OTU4_2nd_replicate.pdf"
#run the treeWAS function
cdl50_binary_OTU4 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50_binary", output_file = outputfile)



#phenotype = median 7dpi-0dpi (median of credible intervals, after the stan model)
#subset to OTU4
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(4),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_7dpi-0dpi/treewas_dpi7_0median_OTU4_2nd_replicate.pdf"
#run the treeWAS function
dpi7_OTU4 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "dpi7_0median", output_file = outputfile)


#phenotype = auc_l (growthcurver )
#subset to OTU4
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(4),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_auc_l/treewas_auc_l_OTU4_2nd_replicate.pdf"
#run the treeWAS function
auc_l_OTU4 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "auc_l", output_file = outputfile)




###OTU3 treeWAS
#phenotype = cdl50 (median of credible intervals, after the stan model)
#subset to 3
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(3),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl_50_OTU3_2nd_replicate.pdf"
#run the treeWAS function
cdl50_OTU3 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50", output_file = outputfile)


#now cdl50 binary as phenotype
#subset to OTU3
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(3),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl50_binary_OTU3_2nd_replicate.pdf"
#run the treeWAS function
cdl50_binary_OTU3 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50_binary", output_file = outputfile)



#phenotype = median 7dpi-0dpi (median of credible intervals, after the stan model)
#subset to OTU3
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(3),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_7dpi-0dpi/treewas_dpi7_0median_OTU3_2nd_replicate.pdf"
#run the treeWAS function
dpi7_OTU3 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "dpi7_0median", output_file = outputfile)


#phenotype = auc_l (growthcurver )
#subset to OTU3
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(3),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_auc_l/treewas_auc_l_OTU3_2nd_replicate.pdf"
#run the treeWAS function
auc_l_OTU3 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "auc_l", output_file = outputfile)





###OTU2 treeWAS
#phenotype = cdl50 (median of credible intervals, after the stan model)
#subset to OTU2
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(2),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl_50_OTU2_2nd_replicate.pdf"
#run the treeWAS function
cdl50_OTU2 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50", output_file = outputfile)


#now cdl50 binary as phenotype
#subset to OTU2
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(2),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_credible_interval_median/treewas_cdl50_binary_OTU2_2nd_replicate.pdf"
#run the treeWAS function
cdl50_binary_OTU2 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "cdl50_binary", output_file = outputfile)



#phenotype = median 7dpi-0dpi (median of credible intervals, after the stan model)
#subset to OTU2
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(2),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_7dpi-0dpi/treewas_dpi7_0median_OTU2_2nd_replicate.pdf"
#run the treeWAS function
dpi7_OTU2 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "dpi7_0median", output_file = outputfile)


#phenotype = auc_l (growthcurver )
#subset to OTU2
OTU_desired <- dataset_GWAS$O_T_U %in% paste("OTU",c(2),sep = "")
strains_desired <- as.character(phy_pattern_my_strains$Species[match(unique(dataset_GWAS$strain[OTU_desired]),phy_pattern_my_strains$Species)])
phy_pattern_my_strains_subset <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_desired,]
#choose output file
outputfile <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/phenotype_auc_l/treewas_auc_l_OTU2_2nd_replicate.pdf"
#run the treeWAS function
auc_l_OTU2 <- treeWAS_by_OTU(tree_path = tree_phyl_all_strains,phyletic_pattern = phy_pattern_my_strains_subset,phenotype = "auc_l", output_file = outputfile)




#now i will summerize all results in a dataframe.
#strating from all strains, the 3 phenotypes
cdl50_all_strains_sum <- as.data.frame(table(unlist(cdl50_all_strains$treeWAS.combined$treeWAS)))
cdl50_all_strains_sum$pheno <- "cdl50"

cdl50_binary_all_strains_sum <- as.data.frame(table(unlist(cdl50_binary_all_strains$treeWAS.combined$treeWAS)))
cdl50_binary_all_strains_sum$pheno <- "cdl50_binary"

dpi7_all_strains_sum <- as.data.frame(table(unlist(dpi7_all_strains$treeWAS.combined$treeWAS)))
dpi7_all_strains_sum$pheno <- "dpi7"

auc_l_all_strains_sum <- as.data.frame(table(unlist(auc_l_all_strains$treeWAS.combined$treeWAS)))
auc_l_all_strains_sum$pheno <- "auc_l"


all_strains_sum <- rbind(cdl50_all_strains_sum,dpi7_all_strains_sum,auc_l_all_strains_sum,cdl50_binary_all_strains_sum)

#OTU2
cdl50_OTU2_sum <- as.data.frame(table(unlist(cdl50_OTU2$treeWAS.combined$treeWAS)))
cdl50_OTU2_sum$pheno <- "cdl50"

cdl50_binary_OTU2_sum <- as.data.frame(table(unlist(cdl50_binary_OTU2$treeWAS.combined$treeWAS)))
cdl50_binary_OTU2_sum$pheno <- "cdl50_binary"

dpi7_OTU2_sum <- as.data.frame(table(unlist(dpi7_OTU2$treeWAS.combined$treeWAS)))
dpi7_OTU2_sum$pheno <- "dpi7"

auc_l_OTU2_sum <- as.data.frame(table(unlist(auc_l_OTU2$treeWAS.combined$treeWAS)))
auc_l_OTU2_sum$pheno <- "auc_l"

OTU2_sum <- rbind(cdl50_OTU2_sum,dpi7_OTU2_sum,auc_l_OTU2_sum,cdl50_binary_OTU2_sum)

#OTU3
cdl50_OTU3_sum <- as.data.frame(table(unlist(cdl50_OTU3$treeWAS.combined$treeWAS)))
cdl50_OTU3_sum$pheno <- "cdl50"

cdl50_binary_OTU3_sum <- as.data.frame(table(unlist(cdl50_binary_OTU3$treeWAS.combined$treeWAS)))
cdl50_binary_OTU3_sum$pheno <- "cdl50_binary"

dpi7_OTU3_sum <- as.data.frame(table(unlist(dpi7_OTU3$treeWAS.combined$treeWAS)))
dpi7_OTU3_sum$pheno <- "dpi7"

auc_l_OTU3_sum <- as.data.frame(table(unlist(auc_l_OTU3$treeWAS.combined$treeWAS)))
auc_l_OTU3_sum$pheno <- "auc_l"

OTU3_sum <- rbind(cdl50_OTU3_sum,dpi7_OTU3_sum,auc_l_OTU3_sum,cdl50_binary_OTU3_sum)


#OTU4
cdl50_OTU4_sum <- as.data.frame(table(unlist(cdl50_OTU4$treeWAS.combined$treeWAS)))
cdl50_OTU4_sum$pheno <- "cdl50"

cdl50_binary_OTU4_sum <- as.data.frame(table(unlist(cdl50_binary_OTU4$treeWAS.combined$treeWAS)))
cdl50_binary_OTU4_sum$pheno <- "cdl50_binary"

dpi7_OTU4_sum <- as.data.frame(table(unlist(dpi7_OTU4$treeWAS.combined$treeWAS)))
dpi7_OTU4_sum$pheno <- "dpi7"

auc_l_OTU4_sum <- as.data.frame(table(unlist(auc_l_OTU4$treeWAS.combined$treeWAS)))
auc_l_OTU4_sum$pheno <- "auc_l"

OTU4_sum <- rbind(cdl50_OTU4_sum,dpi7_OTU4_sum,auc_l_OTU4_sum,cdl50_binary_OTU4_sum)

#now sum it all...
OTU2_sum$OTU <- "OTU2"
OTU3_sum$OTU <- "OTU3"
OTU4_sum$OTU <- "OTU4"
all_strains_sum$OTU <- "all"

all_treeWAS_results <- rbind(OTU2_sum,OTU3_sum,OTU4_sum,all_strains_sum)
all_treeWAS_results <- merge(all_treeWAS_results, cor_spearman_cdl50_merged, by="Var1")

#! write.csv2(x = all_treeWAS_results, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/all_treeWAS_results_2nd_replicate.csv")

##now i will look for OGs that are shared between the test of all strains, and OTU-specific test
OGs_all_strains <- as.character(unique(all_treeWAS_results$Var1[all_treeWAS_results$OTU=="all"]))
OGs_OTU_specific <- as.character(unique(all_treeWAS_results$Var1[all_treeWAS_results$OTU %in% c("OTU2","OTU3","OTU4")]))
OGs_both_OTU_specific_and_all_strains <- all_treeWAS_results[all_treeWAS_results$Var1 %in% intersect(OGs_all_strains,OGs_OTU_specific),]

#! write.csv2(x = OGs_both_OTU_specific_and_all_strains, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/all_treeWAS_results_OGs_shared_between_OTUs_2nd_replicate.csv")


##now i will look for OGs that are shared between different phenotypes
cdl50_dpi7 <- intersect(all_treeWAS_results$Var1[all_treeWAS_results$pheno=="cdl50"], all_treeWAS_results$Var1[all_treeWAS_results$pheno=="dpi7"])
cdl50_aucl <- intersect(all_treeWAS_results$Var1[all_treeWAS_results$pheno=="cdl50"], all_treeWAS_results$Var1[all_treeWAS_results$pheno=="auc_l"])
dpi7_aucl <- intersect(all_treeWAS_results$Var1[all_treeWAS_results$pheno=="dpi7"], all_treeWAS_results$Var1[all_treeWAS_results$pheno=="auc_l"])
cdl50_dpi7_aucl <- intersect(intersect(x = cdl50_aucl,dpi7_aucl),cdl50_aucl)

all_treeWAS_results_shared_pheno <- rbind(data.frame("OG" = cdl50_dpi7_aucl, "intersect"="cdl50_dpi7_aucl"),data.frame("OG" = cdl50_dpi7, "intersect"="cdl50_dpi7"),
      data.frame("OG" = cdl50_aucl, "intersect"="cdl50_aucl"),data.frame("OG" = dpi7_aucl, "intersect"="dpi7_aucl"))

#! write.csv2(x = all_treeWAS_results_shared_pheno, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/all_treeWAS_results_OGs_shared_between_phenotypes_2nd_replicate.csv")



##now i will look for OGs got a significant hit in more than one treeWAS test (treeWAS use three different tests)
treeWAS_2_test <- all_treeWAS_results$Var1[all_treeWAS_results$Freq>1]
all_treeWAS_results_at_least_two_tests <- all_treeWAS_results[all_treeWAS_results$Var1 %in% treeWAS_2_test,]

#! write.csv2(x = all_treeWAS_results_at_least_two_tests, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/all_treeWAS_results_OGs_more_than_one_test_2nd_replicate.csv")


##intersecting shared phenotypes and more than one test OGs
shared_pheno_and_two_test <- intersect(all_treeWAS_results_shared_pheno$OG,all_treeWAS_results_at_least_two_tests$Var1)
all_treewas_results_shared_pheno_and_two_test <- all_treeWAS_results[all_treeWAS_results$Var1 %in% shared_pheno_and_two_test,]
#! write.csv2(x = all_treewas_results_shared_pheno_and_two_test, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/intersection_both_shared_phenotypes_and_more_than_one_test_2nd_replicate.csv")

#intersecting everything
shared_pheno_two_tests_both_OTU_specific_and_all_strains <- intersect(shared_pheno_and_two_test, OGs_both_OTU_specific_and_all_strains$Var1)
all_intersections <- all_treeWAS_results[all_treeWAS_results$Var1 %in% shared_pheno_two_tests_both_OTU_specific_and_all_strains,]

#! write.csv2(x = all_intersections,file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/all_intersections_2nd_replicate.csv")


# all OGs with either shared phenotypes, more than one test, or all vs OTu-specific test.
# these are the OGs that are strong enough to think about knocking them out
strong_OGs <- unique(c(as.character(all_treeWAS_results_shared_pheno$OG),as.character(OGs_both_OTU_specific_and_all_strains$Var1),as.character(all_treeWAS_results_at_least_two_tests$Var1)))

#! write.csv2(x = strong_OGs, file = "OGs_with_either_shared_pheno_OR_OTU_all_and_specific_OR_more_one_test_2nd_replicate.csv")

#now taking into account the spearman correlation, and filtering only for OGs with >0 correlation (only protective)
mask_all_strains <- (all_treeWAS_results$cor_all>0)
mask_OTU2 <- (all_treeWAS_results$cor_OTU2>0 & !is.na(all_treeWAS_results$cor_OTU2))
mask_OTU3 <- (all_treeWAS_results$cor_OTU3>0 & !is.na(all_treeWAS_results$cor_OTU3))
mask_OTU4 <- (all_treeWAS_results$cor_OTU4>0 & !is.na(all_treeWAS_results$cor_OTU4))
all_positive_cor_treeWAS_results <- all_treeWAS_results[mask_all_strains | mask_OTU2 | mask_OTU3 | mask_OTU4,]

# spliting into OTUs, and reordering by spearman correlation
positive_cor_treeWAS_all_strains <- all_positive_cor_treeWAS_results[all_positive_cor_treeWAS_results$OTU=="all" & all_positive_cor_treeWAS_results$cor_all>0,]
positive_cor_treeWAS_OTU2 <- all_positive_cor_treeWAS_results[all_positive_cor_treeWAS_results$OTU=="OTU2"& all_positive_cor_treeWAS_results$cor_OTU2>0,]
positive_cor_treeWAS_OTU3 <- all_positive_cor_treeWAS_results[all_positive_cor_treeWAS_results$OTU=="OTU3"& all_positive_cor_treeWAS_results$cor_OTU3>0,]
positive_cor_treeWAS_OTU4 <- all_positive_cor_treeWAS_results[all_positive_cor_treeWAS_results$OTU=="OTU4"& all_positive_cor_treeWAS_results$cor_OTU4>0,]


positive_cor_treeWAS_all_strains<- positive_cor_treeWAS_all_strains[order(positive_cor_treeWAS_all_strains$cor_all,decreasing = T),]
positive_cor_treeWAS_OTU2<- positive_cor_treeWAS_OTU2[order(positive_cor_treeWAS_OTU2$cor_all,decreasing = T),]
positive_cor_treeWAS_OTU3<- positive_cor_treeWAS_OTU3[order(positive_cor_treeWAS_OTU3$cor_all,decreasing = T),]
positive_cor_treeWAS_OTU4<- positive_cor_treeWAS_OTU4[order(positive_cor_treeWAS_OTU4$cor_all,decreasing = T),]


#! write.csv2(x = positive_cor_treeWAS_all_strains, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_all_strains_2nd_replicate.csv")
#! write.csv2(x = positive_cor_treeWAS_OTU2, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_OTU2_2nd_replicate.csv")
#! write.csv2(x = positive_cor_treeWAS_OTU3, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_OTU3_2nd_replicate.csv")
#! write.csv2(x = positive_cor_treeWAS_OTU4, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_OTU4_2nd_replicate.csv")

###################### now comes the annotation part
OG_annotation <- read.csv("/Volumes/small_projects/hashkenazy/Pseudomonas/List_of_all_Pseudomonas_geneClusters_aa_aln.info.csv")

Majority_annotation_all_results <- OG_annotation$Majority_annotation[match(all_treeWAS_results$Var1,as.character(OG_annotation$Cluster_name))]
all_treeWAS_results$Majority_annotation <- Majority_annotation_all_results
#! write.csv2(x = all_treeWAS_results, "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/all_treeWAS_results_2nd_replicate_annotated.csv")

all_strains_best_hits <- as.character(unique(positive_cor_treeWAS_all_strains$Var1))
OTU2_best_hits <- as.character(unique(positive_cor_treeWAS_OTU2$Var1))
OTU3_best_hits <- as.character(unique(positive_cor_treeWAS_OTU3$Var1))
OTU4_best_hits <- as.character(unique(positive_cor_treeWAS_OTU4$Var1))

Majority_annotation_all_strains <- OG_annotation$Majority_annotation[match(all_strains_best_hits,as.character(OG_annotation$Cluster_name))]
Majority_annotation_OTU2 <- OG_annotation$Majority_annotation[match(OTU2_best_hits,as.character(OG_annotation$Cluster_name))]
Majority_annotation_OTU3 <- OG_annotation$Majority_annotation[match(OTU3_best_hits,as.character(OG_annotation$Cluster_name))]
Majority_annotation_OTU4 <- OG_annotation$Majority_annotation[match(OTU4_best_hits,as.character(OG_annotation$Cluster_name))]

Majority_annotation_all_strains_sum <- data.frame("OG"=all_strains_best_hits,"Majority_annotation"=as.character(Majority_annotation_all_strains))
Majority_annotation_OTU2_sum <- data.frame("OG"=OTU2_best_hits,"Majority_annotation"=as.character(Majority_annotation_OTU2))
Majority_annotation_OTU3_sum <- data.frame("OG"=OTU3_best_hits,"Majority_annotation"=as.character(Majority_annotation_OTU3))
Majority_annotation_OTU4_sum <- data.frame("OG"=OTU4_best_hits,"Majority_annotation"=as.character(Majority_annotation_OTU4))

#! write.csv2(x = Majority_annotation_all_strains_sum, "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_all_strains_Majority_annotation_2nd_replicate.csv")
#! write.csv2(x = Majority_annotation_OTU2_sum, "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_OTU2_Majority_annotation_2nd_replicate.csv")
#! write.csv2(x = Majority_annotation_OTU3_sum, "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_OTU3_Majority_annotation_2nd_replicate.csv")
#! write.csv2(x = Majority_annotation_OTU4_sum, "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/positive_cor_treeWAS_OTU4_Majority_annotation_2nd_replicate.csv")

##################### now plotting the significant candidates with heatmap
source("~/ownCloud/My papers/bacterial_GWAS_paper/scripts/heatmap_significant_hits.R")

#first all strains
heatmap_sig_OGs(OTU_of_interest = "all",sig_hits = all_strains_best_hits, phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/all_strains_positive_cor_treeWAS_heatmap_2nd_replicate.pdf")

#OTU2 - no significant hits. Thus, i used the hits from all strains - which are mainly OTU2
heatmap_sig_OGs(OTU_of_interest = "OTU2",sig_hits = all_strains_best_hits, phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/OTU2_positive_cor_treeWAS_heatmap_2nd_replicate.pdf")

#OTU3
heatmap_sig_OGs(OTU_of_interest = "OTU3",sig_hits = OTU3_best_hits, phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/OTU3_positive_cor_treeWAS_heatmap_2nd_replicate.pdf")
#OTU3 sig hits, among all strains (to see if it shows up in other strains)
heatmap_sig_OGs(OTU_of_interest = "all",sig_hits = OTU3_best_hits, phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/OTU3_positive_cor_treeWAS_heatmap_among_all_strains_2nd_replicate.pdf")


#OTU4
heatmap_sig_OGs(OTU_of_interest = "OTU4",sig_hits = OTU4_best_hits, phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/OTU4_positive_cor_treeWAS_heatmap_2nd_replicate.pdf")
#OTU4 sig hits, among all strains (to see if it shows up in other strains)
heatmap_sig_OGs(OTU_of_interest = "all",sig_hits = OTU4_best_hits, phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/OTU4_positive_cor_treeWAS_heatmap_among_all_strains_2nd_replicate.pdf")

## only best hits - after my manual filtering by spearman correlation per OTU / all strains, and by proximity. Took only the best hits
best_hits_afer_filtering <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/Final_sig_hits_K.O._for_heatmap.csv")

# best OTU2 hits
heatmap_sig_OGs(OTU_of_interest = "all",sig_hits = as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U=="OTU2"]), phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/final_sig_hits_best_OTU2_among_all_strains.pdf")
heatmap_sig_OGs(OTU_of_interest = "OTU2",sig_hits = as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U=="OTU2"]), phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/final_sig_hits_best_OTU2.pdf")



# best OTU3 hits
heatmap_sig_OGs(OTU_of_interest = "all",sig_hits = as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U=="OTU3"]), phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/final_sig_hits_best_OTU3_among_all_strains.pdf")
heatmap_sig_OGs(OTU_of_interest = "OTU3",sig_hits = as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U=="OTU3"]), phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/final_sig_hits_best_OTU3.pdf")

# best OTU4 hits
heatmap_sig_OGs(OTU_of_interest = "all",sig_hits = as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U=="OTU4"]), phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/final_sig_hits_best_OTU4_among_all_strains.pdf")
heatmap_sig_OGs(OTU_of_interest = "OTU4",sig_hits = as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U=="OTU4"]), phy_pattern_my_strains=phy_pattern_my_strains,
                phenotype_continious = "cdl50",output_file_path = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/final_sig_hits_best_OTU4.pdf")



#plotting scatter plot of strains holding significant OGs (results from the treeWAS in this script. This is a manual list that can always be changed)
best_hits_afer_filtering <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/final_filters_and_heatmaps/Final_sig_hits_K.O._for_heatmap.csv")
filter_by_estimation <- c("V","good")
best_hits_afer_filtering <- best_hits_afer_filtering[best_hits_afer_filtering$hits_estimation %in% filter_by_estimation,]


#strating with OTU2
chosen_OTU <- "OTU2"
OG_best_hits <- as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U==chosen_OTU])
phy_pattern_best_hits <- phy_pattern_my_strains[,c(OG_best_hits,"Species")]
rownames(phy_pattern_my_strains) <- phy_pattern_my_strains$Species
phy_best_hits_sums <- data.frame("strain"=phy_pattern_best_hits$Species,"number_of_genes"= as.numeric(rowSums(phy_pattern_best_hits[,OG_best_hits])))
scatter_weight_df$number_of_genes <- phy_best_hits_sums$number_of_genes[match(scatter_weight_df_subset$strain,phy_best_hits_sums$strain)] 
scatter_weight_df$number_of_genes[is.na(scatter_weight_df$number_of_genes)] <- 0

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/scatter_plot_sig_hits/OTU2_sig_hits_estimation_good_V.pdf")
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes), 
       data = scatter_weight_df[scatter_weight_df$OTU==chosen_OTU,]) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  scale_colour_gradient(low = "grey", high = "green")+
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/scatter_plot_sig_hits/OTU2_sig_hits_estimation_good_V_among_all_strains.pdf")
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes), 
       data = scatter_weight_df) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  scale_colour_gradient(low = "grey", high = "green")+
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()



#Now OTU4
chosen_OTU <- "OTU4"
OG_best_hits <- as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U==chosen_OTU])
phy_pattern_best_hits <- phy_pattern_my_strains[,c(OG_best_hits,"Species")]
rownames(phy_pattern_my_strains) <- phy_pattern_my_strains$Species
phy_best_hits_sums <- data.frame("strain"=phy_pattern_best_hits$Species,"number_of_genes"= as.numeric(rowSums(phy_pattern_best_hits[,OG_best_hits])))
scatter_weight_df$number_of_genes <- phy_best_hits_sums$number_of_genes[match(scatter_weight_df_subset$strain,phy_best_hits_sums$strain)] 
scatter_weight_df$number_of_genes[is.na(scatter_weight_df$number_of_genes)] <- 0

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/scatter_plot_sig_hits/OTU4_sig_hits_estimation_good_V.pdf")
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes), 
       data = scatter_weight_df[scatter_weight_df$OTU==chosen_OTU,]) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  scale_colour_gradient(low = "gray", high = "green")+
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/treewas/Summary/scatter_plot_sig_hits/OTU4_sig_hits_estimation_good_V_among_all_strains.pdf")
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes), 
       data = scatter_weight_df) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  scale_colour_gradient(low = "grey", high = "green")+
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()


#Now OTU3
chosen_OTU <- "OTU3"
OG_best_hits <- as.character(best_hits_afer_filtering$OG[best_hits_afer_filtering$O_T_U==chosen_OTU])
phy_pattern_best_hits <- phy_pattern_my_strains[,c(OG_best_hits,"Species")]
rownames(phy_pattern_my_strains) <- phy_pattern_my_strains$Species
phy_best_hits_sums <- data.frame("strain"=phy_pattern_best_hits$Species,"number_of_genes"= as.numeric(rowSums(phy_pattern_best_hits[,OG_best_hits])))
scatter_weight_df$number_of_genes <- phy_best_hits_sums$number_of_genes[match(scatter_weight_df_subset$strain,phy_best_hits_sums$strain)] 
scatter_weight_df$number_of_genes[is.na(scatter_weight_df$number_of_genes)] <- 0

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS_few_pathogens_by_protective/figures/1st_rep_p4.C9.pdf")
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes), 
       data = scatter_weight_df[scatter_weight_df$OTU==chosen_OTU,]) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=number_of_genes)) +
  scale_colour_gradient(low = "grey", high = "green")+
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()


