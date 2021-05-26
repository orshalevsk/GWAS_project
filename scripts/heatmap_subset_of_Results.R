'''
Presenting significant results in heatmap, clustered by the collection of strains of interest.
'''
library("devtools")
library("lme4")
library("lmtest")
library(ggplot2)
library(MASS)
library(FSA)
library(agricolae)
require(vegan)
library(multcomp)
require(nlme)
library(randomForest)
library(tidyr)
library(treeWAS)
library(ape)
library(adegenet)
#library(emma)
library(gplots)
library(plyr)
library(RColorBrewer)
library(tiger)

# read experimental results and metadata

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
subset_siginifact_OGs <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/OTU2_significant_OGs.csv", header = F) #subset to best OTU2 hits (signal from all strains)


# Add OTU for subsequent subsetting (e.g. heatmap of only OTU2 strains)
dataset_GWAS_delta_pixels$OTU <- metadata$O_T_U[match(dataset_GWAS_delta_pixels$strain,metadata$strain)]
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)

dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="control"] <- "Bacteria free" #adding "OTU" level to control for subsequent plotting
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$OTU==""] <- "OTU5" # adding OTU5



#########plotting heatmap, using Haims script. a function so i can do it repeatedly

heatmap_sig_OGs <- function(OTU_of_interest ="all", sig_hits, phenotype_continious, output_file_path, phy_pattern_my_strains){
  ScriptPath="~/ownCloud/My papers/bacterial_GWAS_paper/scripts/"
  source (paste (ScriptPath,"plot_phyletic_and_trait.R",sep=""))
  
  #store all possible trees path
  tree_phyl_all_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_my_strains_aligned_no_p4.C9.fasta.treefile"
  tree_phyl_non_OTU2_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_non_OTU2_strains_aligned.fasta.treefile"
  tree_phyl_OTU2_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_OTU2_strains_aligned.fasta.treefile"
  tree_phyl_OTU4_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_OTU4_strains_aligned.fasta.treefile"
  tree_phyl_OTU3_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_OTU3_strains_aligned.fasta.treefile"
  
  
  #choose OTU to work with
  if (OTU_of_interest=="all"){
    tree=  tree_phyl_all_strains
    OTU <- paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
  } else if (OTU_of_interest=="OTU2"){
    tree=  tree_phyl_OTU2_strains
    OTU <- "OTU2"
    
  } else if (OTU_of_interest=="OTU3"){
    tree=  tree_phyl_OTU3_strains
    OTU <- "OTU3"
    
  } else if (OTU_of_interest=="OTU4"){
    tree=  tree_phyl_OTU4_strains
    OTU <- "OTU4"
    
  }
  
  strains_by_OTU <- unique(dataset_GWAS_delta_pixels$strain[dataset_GWAS_delta_pixels$OTU %in% OTU])
  strains_by_OTU <- strains_by_OTU[!is.na(strains_by_OTU)]
  phy_pattern_my_strains_by_OTU <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_by_OTU,]
  
  OTU <- metadata$O_T_U[metadata$strain %in% strains_by_OTU & !(metadata$GWAS_index%in%c("control","boiled_strain_17+p"))]
  strain <- metadata$strain[metadata$strain %in% strains_by_OTU & !(metadata$GWAS_index%in%c("control","boiled_strain_17+p"))]
  OTU <- as.character(OTU)
  OTU[OTU==""] <- "OTU5"
  OTU_strain <- data.frame("strain"=strain,"OTU"=OTU)
  
  # build a phyletic pattern of only chosen OGs
  best_hits_phy_pattern <- phy_pattern_my_strains_by_OTU[,as.character(sig_hits[[1]])]
  row.names(best_hits_phy_pattern) <- phy_pattern_my_strains_by_OTU$Species
  
  
  #create a phenotype df with colors by the phenotype (required from this script, made by Haim)
  my_pheno <- as.data.frame(phy_pattern_my_strains_by_OTU[,phenotype_continious])
  colnames(my_pheno) <- "phenotype"
  row.names(my_pheno) <- phy_pattern_my_strains_by_OTU$Species
  my_pheno <- my_pheno[order(my_pheno$phenotype), , drop=F]
  #OPTIONAL: if green pixels of cdl50 is less than 0 => make it 0
  #my_pheno$phenotype[my_pheno$phenotype<0] <- 0
  
  my_pheno$color <- color.factor("darkgreen",(my_pheno$phenotype-min(my_pheno$phenotype))^3, max(my_pheno$phenotype-min(my_pheno$phenotype))^3)
  
  # color tree branches by OTU

  #plotting (this will create a file)
  plot_phyletic_stretch_and_trait(PhyP.df = as.matrix(best_hits_phy_pattern), TreeFile = tree, Traits_Colors_df = my_pheno,outFile = output_file_path,dend_height = 10,dend_width = 7)

}

#! output_path="~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/heatmap_OTU2_sig_all_strains_temp.pdf"
heatmap_sig_OGs(OTU_of_interest = "all",sig_hits = subset_siginifact_OGs,phenotype_continious = "median_delta_pxls", 
                phy_pattern_my_strains = phy_pattern_my_strains, output_file_path = output_path)


#! output_path="~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/heatmap_OTU2_sig_OTU2_strains.pdf"
heatmap_sig_OGs(OTU_of_interest = "OTU2",sig_hits = subset_siginifact_OGs,phenotype_continious = "median_delta_pxls", phy_pattern_my_strains = phy_pattern_my_strains, output_file_path = output_path)



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

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/OTU2_sig_OG_growth_scatter_plot.pdf")
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


#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure3/OTU2_sig_OG_growth_scatter_plot_all_strains.pdf")
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
