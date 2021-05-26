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


dataset_GWAS$pxls_7dpi_minus_0dpi[dataset_GWAS$dpi==7] <- dataset_GWAS$pxls[dataset_GWAS$dpi==7]-dataset_GWAS$pxls[dataset_GWAS$dpi==0]
dataset_GWAS_delta_pixels <- dataset_GWAS[dataset_GWAS$dpi==7,]


table(dataset_GWAS_delta_pixels$treatment[dataset_GWAS_delta_pixels$pxls_7dpi_minus_0dpi<(-4000)])


#dataset_GWAS_protective <- dataset_GWAS_delta_pixels[dataset_GWAS_delta_pixels$treatment %in% c("control", "pathogen_only", "97+p", "boiled_strain_17+p", "8+p", "30+p", "32+p"),]

dataset_GWAS_protective <- dataset_GWAS_delta_pixels[dataset_GWAS_delta_pixels$treatment %in% c("control", "97+p", "8+p", "30+p"),]

dataset_GWAS_protective$treatment <- factor(dataset_GWAS_protective$treatment, levels=c("control", "97+p", "8+p", "30+p"))

levels(dataset_GWAS_protective$treatment) <- c("control","pathogen", "pathogen +\nheat-killed protective", "protective")

summary(lm(formula = pxls_7dpi_minus_0dpi ~ treatment ,data = dataset_GWAS_protective))

ggplot(data = dataset_GWAS_protective, aes(x = treatment, y= pxls_7dpi_minus_0dpi))+
  #geom_boxplot()+
  geom_point() +
  theme_classic()


fit_sig_hits <- stan_glm(formula = pxls_7dpi_minus_0dpi ~ treatment, data = dataset_GWAS_protective)

grid = dataset_GWAS_protective %>%
  data_grid(treatment)

fits = grid %>%
  add_fitted_draws(fit_sig_hits)

preds = grid %>%
  add_predicted_draws(fit_sig_hits)



#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/controls_stats.pdf", useDingbats = F)
dataset_GWAS_protective %>%
  ggplot(aes(y = treatment, x = pxls_7dpi_minus_0dpi)) +
  stat_interval(aes(x = .prediction), data = preds) +
  stat_pointinterval(aes(x = .value), data = fits, .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()


min_control <- min(dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="control"])

dataset_GWAS_protective$dead_alive[dataset_GWAS_protective$pxls_7dpi_minus_0dpi >= min_control] <- "normal"
dataset_GWAS_protective$dead_alive[dataset_GWAS_protective$pxls_7dpi_minus_0dpi < min_control & dataset_GWAS_protective$pxls_7dpi_minus_0dpi > 0 ] <- "sick"
dataset_GWAS_protective$dead_alive[dataset_GWAS_protective$pxls_7dpi_minus_0dpi < 0] <- "dead"
dataset_GWAS_protective$dead_alive <- factor(dataset_GWAS_protective$dead_alive, levels = c("dead","sick","normal"))

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/distribution_dead_sick_normal_control_treatments_0_dead_min_normal_sick.pdf", useDingbats = F)
ggplot(data = dataset_GWAS_protective, aes(x = treatment, fill=dead_alive))+
  geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("darkred","darkgoldenrod1","grey70")) +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  theme_classic()
#! dev.off()






two.group.unpaired <- dabest(dataset_GWAS_protective , treatment, pxls_7dpi_minus_0dpi, 
                             # The idx below passes "Control" as the control group, 
                             # and "Group1" as the test group. The mean difference
                             # will be computed as mean(Group1) - mean(Control1).
                             idx = list(c("control","pathogen_only","boiled_strain_17+p","30+p")),
                             paired = FALSE) %>%
  mean_diff()

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure2/Protection_by_OTU.pdf")
plot(two.group.unpaired)
#! dev.off()

t.test(dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="30+p"],dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="pathogen_only"])
t.test(dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="control"],dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="pathogen_only"])




mean(dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="control"])
mean(dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="pathogen_only"])



mean(dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="pathogen_only"])
sd(dataset_GWAS_protective$pxls_7dpi_minus_0dpi[dataset_GWAS_protective$treatment=="pathogen_only"])

