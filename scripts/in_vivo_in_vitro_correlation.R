library(ggplot2)
library(dabestr)
library(rstanarm)
library(ggplot2)
library(dabestr)
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


# creating 7-0 pxls dataset, to bin by OTU
dataset_GWAS$pxls_7dpi_minus_0dpi[dataset_GWAS$dpi==7] <- dataset_GWAS$pxls[dataset_GWAS$dpi==7]-dataset_GWAS$pxls[dataset_GWAS$dpi==0]
dataset_GWAS_delta_pixels <- dataset_GWAS[dataset_GWAS$dpi==7,]
dataset_GWAS_delta_pixels <- dataset_GWAS_delta_pixels[!dataset_GWAS_delta_pixels$treatment%in% c("boiled_strain_17+p","pathogen_only"),] # remove boiled protector and pathogen


# adding OTU
dataset_GWAS_delta_pixels$OTU <- metadata$O_T_U[match(dataset_GWAS_delta_pixels$strain,metadata$strain)]
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)

#removing control treatments
dataset_GWAS_delta_pixels <- dataset_GWAS_delta_pixels[!dataset_GWAS_delta_pixels$treatment %in% c("control","boiled_strain_17+p","pathogen_only"),]


dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$OTU==""] <- "OTU5" # adding OTU5
dataset_GWAS_delta_pixels$OTU <- gsub(pattern = "OTU",replacement = "ATUE",x = dataset_GWAS_delta_pixels$OTU)
dataset_GWAS_delta_pixels$OTU <- factor(dataset_GWAS_delta_pixels$OTU, levels = c("ATUE2", "ATUE3", "ATUE4", "ATUE5", "ATUE6","ATUE7","ATUE8","ATUE9","ATUE10","ATUE11","ATUE_NEW","Bacteria free","Pathogen+boiled_protector","Pathogen"))


for (strain in unique(dataset_GWAS_delta_pixels$strain)){
  dataset_GWAS_delta_pixels$median_pxls[dataset_GWAS_delta_pixels$strain==strain] <- median(dataset_GWAS_delta_pixels$pxls_7dpi_minus_0dpi[dataset_GWAS_delta_pixels$strain==strain])
}

dataset_GWAS_delta_pixels_median <- unique(dataset_GWAS_delta_pixels[,c("strain","OTU","median_pxls")])
in_vitro <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/in_vitro.csv")

dataset_GWAS_delta_pixels_median$in_vitro <- as.numeric(in_vitro$in_vitro[match(dataset_GWAS_delta_pixels_median$strain, in_vitro$strain)])
dataset_GWAS_delta_pixels_median$in_vitro[is.na(dataset_GWAS_delta_pixels_median$in_vitro)] <- 0


plot(dataset_GWAS_delta_pixels_median$median_pxls, dataset_GWAS_delta_pixels_median$in_vitro)

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/in_vitro_in_vivo_cor.pdf", useDingbats = F)
ggplot(data = dataset_GWAS_delta_pixels_median, aes(x = median_pxls, y = in_vitro)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()
#! dev.off()

summary(lm(dataset_GWAS_delta_pixels_median$median_pxls ~ dataset_GWAS_delta_pixels_median$in_vitro))

cor.test(dataset_GWAS_delta_pixels_median$median_pxls,dataset_GWAS_delta_pixels_median$in_vitro,method = "pearson")
cor.test(dataset_GWAS_delta_pixels_median$median_pxls,dataset_GWAS_delta_pixels_median$in_vitro,method = "spearman")
cor.test(dataset_GWAS_delta_pixels_median$median_pxls,dataset_GWAS_delta_pixels_median$in_vitro,method = "kendall")

dataset_GWAS_delta_pixels_median$in_vitro_binary[dataset_GWAS_delta_pixels_median$in_vitro > 0] <- 1
dataset_GWAS_delta_pixels_median$in_vitro_binary[dataset_GWAS_delta_pixels_median$in_vitro == 0] <- 0

dataset_GWAS_delta_pixels_median$in_vitro_binary <- as.factor(dataset_GWAS_delta_pixels_median$in_vitro_binary)


ggplot(data = dataset_GWAS_delta_pixels_median[dataset_GWAS_delta_pixels_median$OTU=="ATUE2",], aes(x = in_vitro_binary, y = median_pxls)) + 
  geom_point() +
  geom_boxplot() +
  theme_classic()

dataset_GWAS_delta_pixels_median_ATUE2 <- dataset_GWAS_delta_pixels_median[dataset_GWAS_delta_pixels_median$OTU=="ATUE2",]

wilcox.test(dataset_GWAS_delta_pixels_median_ATUE2$median_pxls[dataset_GWAS_delta_pixels_median_ATUE2$in_vitro_binary==1],dataset_GWAS_delta_pixels_median_ATUE2$median_pxls[dataset_GWAS_delta_pixels_median_ATUE2$in_vitro_binary==0])


ggplot(data = dataset_GWAS_delta_pixels_median, aes(x = in_vitro_binary, y = median_pxls)) + 
  geom_point() +
  geom_boxplot() +
  theme_classic()


wilcox.test(dataset_GWAS_delta_pixels_median$median_pxls[dataset_GWAS_delta_pixels_median$in_vitro_binary==1],dataset_GWAS_delta_pixels_median$median_pxls[dataset_GWAS_delta_pixels_median$in_vitro_binary==0])

fit <- stan_glm(formula = median_pxls ~ in_vitro_binary, data = dataset_GWAS_delta_pixels_median)

grid = dataset_GWAS_delta_pixels_median %>%
  data_grid(in_vitro_binary)

fits = grid %>%
  add_fitted_draws(fit)

preds = grid %>%
  add_predicted_draws(fit)



#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/In_vitro_in_vivo_stan_model.pdf", useDingbats = F)
dataset_GWAS_delta_pixels_median %>%
  ggplot(aes(y = in_vitro_binary, x = median_pxls)) +
  stat_interval(aes(x = .prediction), data = preds) +
  stat_pointinterval(aes(x = .value), data = fits, .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()

'''
dataset_GWAS_delta_pixels_median$in_vitro <- as.factor(dataset_GWAS_delta_pixels_median$in_vitro)
fit1 <- stan_glm(formula = median_pxls ~ in_vitro, data = dataset_GWAS_delta_pixels_median)


grid = dataset_GWAS_delta_pixels_median %>%
  data_grid(in_vitro)

fits = grid %>%
  add_fitted_draws(fit1)

preds = grid %>%
  add_predicted_draws(fit1)



#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/In_vitro_in_vivo_stan_model.pdf", useDingbats = F)
dataset_GWAS_delta_pixels_median %>%
  ggplot(aes(y = in_vitro, x = median_pxls)) +
  stat_interval(aes(x = .prediction), data = preds) +
  stat_pointinterval(aes(x = .value), data = fits, .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()
'''