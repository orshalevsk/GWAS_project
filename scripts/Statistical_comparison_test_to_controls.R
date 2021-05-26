library(ggplot2)
library(dabestr)
library(rstanarm)
library(bayesplot)

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


# adding OTU
dataset_GWAS_delta_pixels$OTU <- metadata$O_T_U[match(dataset_GWAS_delta_pixels$strain,metadata$strain)]
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)

#renaming control treatments
#dataset_GWAS_delta_pixels$strain[dataset_GWAS_delta_pixels$treatment %in% c("boiled_strain_17+p")] <- "boiled+pathogen"
#dataset_GWAS_delta_pixels$strain[dataset_GWAS_delta_pixels$treatment %in% c("pathogen_only")] <- "pathogen"

#removing pathogen/boiled_pathogen treatment

dataset_GWAS_delta_pixels <- dataset_GWAS_delta_pixels[!dataset_GWAS_delta_pixels$treatment %in% c("boiled_strain_17+p", "pathogen_only"),]


dataset_GWAS_delta_pixels$strain[is.na(dataset_GWAS_delta_pixels$strain)] <- "a_control"


dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$OTU==""] <- "OTU5" # adding OTU5
dataset_GWAS_delta_pixels$OTU <- gsub(pattern = "OTU",replacement = "ATUE",x = dataset_GWAS_delta_pixels$OTU)
dataset_GWAS_delta_pixels$OTU <- factor(dataset_GWAS_delta_pixels$OTU, levels = c("ATUE2", "ATUE3", "ATUE4", "ATUE5", "ATUE6","ATUE7","ATUE8","ATUE9","ATUE10","ATUE11","ATUE_NEW","Bacteria free","Pathogen+boiled_protector","Pathogen"))
dataset_GWAS_delta_pixels$strain <- factor(dataset_GWAS_delta_pixels$strain)

#fit1 <- summary(lm(formula = pxls_7dpi_minus_0dpi ~ strain + low_concentration, data = dataset_GWAS_delta_pixels))


fit <- stan_glm(formula = pxls_7dpi_minus_0dpi ~ strain + plate,data =  dataset_GWAS_delta_pixels,iter=10000)
#posterior <- as.matrix(fit)
#plot(fit, regex_pars = c("strain"))


fit_df <- as.data.frame(fit$stan_summary[,c("mean","2.5%","97.5%")])
fit_df$strain <- rownames(fit_df)
fit_df <- fit_df[!fit_df$strain %in% c("plate","sigma","mean_PPD","log-posterior","(Intercept)"),]
level_order <- fit_df$strain[order(fit_df$mean)]
fit_df$strain <- factor(fit_df$strain,levels = level_order)
fit_df$strain_function[fit_df$`2.5%`>0] <- "growth_promotion"
fit_df$strain_function[fit_df$`2.5%`<0 & fit_df$`97.5%`>0] <- "protection"
fit_df$strain_function[fit_df$`97.5%`<0] <- "no protection"



#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Mean_difference_to_control_only_test_strains.pdf", width = 10, height = 10, useDingbats = F)
ggplot(data = fit_df, aes(x=strain,y=mean, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(aes(group=strain, color=strain_function),position=position_dodge(width=1)) +
  theme_classic() +
  scale_color_manual(values = c("#1b9e77","#d95f02","black")) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Mean difference to control (green pixels)") + 
  xlab("Strain")
#! dev.off()

#adding control to the OTU
dataset_GWAS_delta_pixels$OTU <- as.character(dataset_GWAS_delta_pixels$OTU)
dataset_GWAS_delta_pixels$OTU[is.na(dataset_GWAS_delta_pixels$OTU)] <- "a_control"
dataset_GWAS_delta_pixels$OTU <- factor(dataset_GWAS_delta_pixels$OTU)

fit1 <- stan_glm(formula = pxls_7dpi_minus_0dpi ~ OTU,data =  dataset_GWAS_delta_pixels,iter=10000)
#posterior <- as.matrix(fit11)
#plot(fit11, regex_pars = c("strain"))


fit1_df <- as.data.frame(fit1$stan_summary[,c("mean","2.5%","97.5%")])
fit1_df$strain <- rownames(fit1_df)
fit1_df <- fit1_df[!fit1_df$strain %in% c("plate","sigma","mean_PPD","log-posterior","(Intercept)"),]
level_order <- fit1_df$strain[order(fit1_df$mean)]
fit1_df$strain <- factor(fit1_df$strain,levels = level_order)
fit1_df$strain_function[fit1_df$`2.5%`>0] <- "growth_promotion"
fit1_df$strain_function[fit1_df$`2.5%`<0 & fit1_df$`97.5%`>0] <- "protection"
fit1_df$strain_function[fit1_df$`97.5%`<0] <- "no protection"



#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Mean_difference_to_control_only_test_strains.pdf", width = 10, height = 10, useDingbats = F)
ggplot(data = fit1_df, aes(x=strain,y=mean, ymin = `2.5%`, ymax = `97.5%`))+
  geom_pointrange(aes(group=strain, color=strain_function),position=position_dodge(width=1)) +
  theme_classic() +
  scale_color_manual(values = c("#d95f02","black")) +
  geom_hline(yintercept=0, linetype="dashed") +
  ylab("Mean difference to control (green pixels)") + 
  xlab("Strain")
#! dev.off()


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

dataset_GWAS_delta_pixels$OTU <- factor(dataset_GWAS_delta_pixels$OTU, levels=c("a_control","ATUE2", "ATUE3", "ATUE4", "ATUE5", "ATUE6","ATUE7","ATUE8","ATUE9","ATUE10","ATUE11","ATUE_NEW"))

grid = dataset_GWAS_delta_pixels %>%
  data_grid(OTU)

fits = grid %>%
  add_fitted_draws(fit1)

preds = grid %>%
  add_predicted_draws(fit1)


minimum_control <- quantile(fits$.value[fits$OTU=="a_control"],probs = (1/6))
maximum_control <- quantile(fits$.value[fits$OTU=="a_control"],probs = (5/6))


#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/OTU_wise_protection_plot.pdf", useDingbats = F)
dataset_GWAS_delta_pixels[!dataset_GWAS_delta_pixels$OTU=="a_control",] %>%
  ggplot(aes(y = OTU, x = pxls_7dpi_minus_0dpi)) +
  stat_interval(aes(x = .prediction), data = preds[preds$OTU!="a_control",]) +
  stat_pointinterval(aes(x = .value), data = fits[fits$OTU!="a_control",], .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  geom_vline(xintercept = minimum_control, linetype="dotted", 
             color = "blue")+
  geom_vline(xintercept = maximum_control, linetype="dotted", 
             color = "blue")+
  theme_classic()
#! dev.off()

#print stats to include in the text
quantile(fits$.value[fits$OTU=="a_control"],probs = (1/6))
quantile(fits$.value[fits$OTU=="a_control"],probs = (1/2))
quantile(fits$.value[fits$OTU=="a_control"],probs = (5/6))

for (OTU in c("ATUE2","ATUE3","ATUE4")){
  print(OTU)
  print(quantile(fits$.value[fits$OTU==OTU],probs = 0.025))
  print(quantile(fits$.value[fits$OTU==OTU],probs = 0.5))
  print(quantile(fits$.value[fits$OTU==OTU],probs = 0.975))
  
}




#summary(lm(formula = pxls_7dpi_minus_0dpi ~ OTU, data = dataset_GWAS_delta_pixels))





