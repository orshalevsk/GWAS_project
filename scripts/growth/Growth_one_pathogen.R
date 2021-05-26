library(ggplot2)
library(dabestr)

#loading April exp, and subsetting to treatments Control, 1, 2 and 3.
default_dir <- "~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS/results/Green_pixels/"
setwd(default_dir)
# the "temp" file is the right one!!
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

#scatter plot, no boiled_pathogens

#!pdf("/Volumes/small_projects/oshalev/GWAS_project/Figures/full_Exp_2019/all_strains/green_pixels_scatter_plot_over_time.pdf")
#scatter_weight_df_subset <- scatter_weight_df[scatter_weight_df$dpi%in%c("1","4","7","6") & scatter_weight_df$treatment%in%c("pathogen_only","boiled_strain_17+p","control","30+p"),]
scatter_weight_df_subset <- scatter_weight_df[scatter_weight_df$dpi%in%c("0","1","4","5","7","6") ,]


levels(scatter_weight_df_subset$treatment)

colors_scatter_by_controls <- rep("#E0E0E0",length(scatter_weight_df_subset$treatment))
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df_subset$treatment))=="control",6)] <- "#000000"
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df_subset$treatment))=="boiled_strain_17+p",6)] <- "#7570b3"
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df_subset$treatment))=="pathogen_only",6)] <- "#d95f02"
colors_scatter_by_controls[rep(levels(factor(scatter_weight_df_subset$treatment))=="30+p",6)] <- "#1b9e77"

#! pdf("/Users/oshalev/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure2/Growth_scatter_plot_protective_strains_coninfections_P4.C9.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_subset) +
  geom_line() +
  geom_line(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), 
            data = scatter_weight_df_subset[scatter_weight_df_subset$treatment%in%c("pathogen_only","control","30+p","boiled_strain_17+p"),]) +
  #geom_pointrange(aes(ymin=pxls_average-sem, ymax=pxls_average+sem)) + #for sem
  geom_point()+
  geom_point(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), 
             data = scatter_weight_df_subset[scatter_weight_df_subset$treatment%in%c("pathogen_only","control","30+p","boiled_strain_17+p"),])+
  scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  theme(legend.position="none")
#! dev.off()


# creating 7-0 pxls dataset, to bin by OTU
dataset_GWAS$pxls_7dpi_minus_0dpi[dataset_GWAS$dpi==7] <- dataset_GWAS$pxls[dataset_GWAS$dpi==7]-dataset_GWAS$pxls[dataset_GWAS$dpi==0]
dataset_GWAS_delta_pixels <- dataset_GWAS[dataset_GWAS$dpi==7,]
dataset_GWAS_delta_pixels <- dataset_GWAS_delta_pixels[!dataset_GWAS_delta_pixels$treatment%in% c("boiled_strain_17+p","pathogen_only"),] # remove boiled protector and pathogen


# adding OTU
dataset_GWAS_delta_pixels$OTU <- metadata$O_T_U[match(dataset_GWAS_delta_pixels$strain,metadata$strain)]
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)



dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="control"] <- "Bacteria free" #adding "OTU" level to control for subsequent plotting
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$OTU==""] <- "OTU5" # adding OTU5
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="boiled_strain_17+p"] <- "Pathogen+boiled_protector"
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="pathogen_only"] <- "Pathogen"
dataset_GWAS_delta_pixels$OTU <- gsub(pattern = "OTU",replacement = "ATUE",x = dataset_GWAS_delta_pixels$OTU)
dataset_GWAS_delta_pixels$OTU <- factor(dataset_GWAS_delta_pixels$OTU, levels = c("ATUE2", "ATUE3", "ATUE4", "ATUE5", "ATUE6","ATUE7","ATUE8","ATUE9","ATUE10","ATUE11","ATUE_NEW","Bacteria free","Pathogen+boiled_protector","Pathogen"))


#! pdf("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS_few_pathogens_by_protective/figures/1st_rep_P4.C9_7dpi-0dpi.pdf")
ggplot(data = dataset_GWAS_delta_pixels[dataset_GWAS_delta_pixels$OTU %in% c("ATUE2","ATUE3","ATUE4","Bacteria free","Pathogen"),], aes(x = OTU, y=pxls_7dpi_minus_0dpi))+
  geom_boxplot()+
  #geom_point()+
  theme_bw()
#! dev.off()


two.group.unpaired <- dabest(dataset_GWAS_delta_pixels , OTU, pxls_7dpi_minus_0dpi, 
         # The idx below passes "Control" as the control group, 
         # and "Group1" as the test group. The mean difference
         # will be computed as mean(Group1) - mean(Control1).
         idx = list(c("Bacteria free","ATUE2", "ATUE3", "ATUE4", "ATUE5", "ATUE6","ATUE7","ATUE8","ATUE9","ATUE10","ATUE11","ATUE_NEW")),
         paired = FALSE) %>%
  mean_diff()

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure2/Protection_by_OTU.pdf")
plot(two.group.unpaired)
#! dev.off()


two.group.unpaired <- dabest(dataset_GWAS_delta_pixels , OTU, pxls_7dpi_minus_0dpi, 
                             # The idx below passes "Control" as the control group, 
                             # and "Group1" as the test group. The mean difference
                             # will be computed as mean(Group1) - mean(Control1).
                             idx = list(c("Bacteria free","ATUE2", "ATUE3", "ATUE4")),
                             paired = FALSE) %>%
  mean_diff()

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure2/Protection_by_OTU_main_OTUs.pdf")
plot(two.group.unpaired)
#! dev.off()


