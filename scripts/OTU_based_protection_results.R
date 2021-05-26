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

# creating 7-0 pxls dataset, to bin by OTU
dataset_GWAS$pxls_7dpi_minus_0dpi[dataset_GWAS$dpi==7] <- dataset_GWAS$pxls[dataset_GWAS$dpi==7]-dataset_GWAS$pxls[dataset_GWAS$dpi==0]
dataset_GWAS_delta_pixels <- dataset_GWAS[dataset_GWAS$dpi==7,]
#dataset_GWAS_delta_pixels <- dataset_GWAS_delta_pixels[!dataset_GWAS_delta_pixels$treatment%in% c("boiled_strain_17+p","pathogen_only"),] # remove boiled protector and pathogen


# adding OTU
dataset_GWAS_delta_pixels$OTU <- metadata$O_T_U[match(dataset_GWAS_delta_pixels$strain,metadata$strain)]
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)



dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="control"] <- "Bacteria free" #adding "OTU" level to control for subsequent plotting
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$OTU==""] <- "OTU5" # adding OTU5
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="boiled_strain_17+p"] <- "Pathogen+boiled_protector"
dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$treatment=="pathogen_only"] <- "Pathogen"
dataset_GWAS_delta_pixels$OTU <- gsub(pattern = "OTU",replacement = "ATUE",x = dataset_GWAS_delta_pixels$OTU)
dataset_GWAS_delta_pixels$OTU <- factor(dataset_GWAS_delta_pixels$OTU, levels = c("Bacteria free","ATUE2", "ATUE3", "ATUE4", "ATUE5", "ATUE6","ATUE7","ATUE8","ATUE9","ATUE10","ATUE11","ATUE_NEW","Pathogen+boiled_protector","Pathogen"))


ggplot(data = dataset_GWAS_delta_pixels, aes(x = OTU, y=pxls_7dpi_minus_0dpi))+
  geom_violin()+
  #geom_point()+
  theme_bw()



fit <- stan_glm(formula = pxls_7dpi_minus_0dpi ~ OTU, data = dataset_GWAS_delta_pixels)
fit$stan_summary


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

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure2/Protection_by_OTU.pdf", useDingbats = F, width = 15, height = 10)
plot(two.group.unpaired)
#! dev.off()


two.group.unpaired <- dabest(dataset_GWAS_delta_pixels , OTU, pxls_7dpi_minus_0dpi, 
                             # The idx below passes "Control" as the control group, 
                             # and "Group1" as the test group. The mean difference
                             # will be computed as mean(Group1) - mean(Control1).
                             idx = list(c("Bacteria free","ATUE2", "ATUE3", "ATUE4")),
                             paired = FALSE) %>%
  mean_diff()

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure2/Protection_by_OTU_main_OTUs.pdf", useDingbats = F)
plot(two.group.unpaired)
#! dev.off()