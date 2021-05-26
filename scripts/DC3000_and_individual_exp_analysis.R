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
results <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/results_20210318.csv")
results$Position <- paste(results$Plate, results$Pot, sep = "_")

metadata <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/plan/first_Exp_G2/metadata_treamtnets_G2_first_exp_18032021.csv")
metadata$Position <- paste("p", metadata$Position, sep = "")

final_data <- merge(x = results, y = metadata, by = "Position")
final_data$Green_pixels
final_data$Treatment <- as.character(final_data$Treatment)
final_data$Treatment[final_data$Treatment==" Control"] <- " control"
final_data$Treatment <- gsub(pattern = " ", replacement = "", x = final_data$Treatment)
final_data$Treatment <- as.factor(final_data$Treatment)

final_data$delta_Green_pixels[final_data$day==7] <- final_data$Green_pixels[final_data$day==7]-final_data$Green_pixels[final_data$day==0]

final_data_delta <- final_data[final_data$day==7,]

final_data$delta_Green_pixels

ggplot(data = final_data[final_data$day==7,], aes(x = Treatment, y = delta_Green_pixels))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_boxplot()



final_data_patho <- final_data[(grepl(x = final_data$Treatment, pattern = "patho")) | final_data$Treatment == "control" | final_data$Treatment == "Pathogen_(p4.C9)",]

ggplot(data = final_data_patho[final_data_patho$day==7,], aes(x = Treatment, y = delta_Green_pixels))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_boxplot()



##scatter plot the data
final_data$treatment <- final_data$Treatment
final_data$dpi <- final_data$day
final_data$pxls <- final_data$Green_pixels

source("~/ownCloud/My papers/bacterial_GWAS_paper/scripts/data_to_scatter_plot.R")
scatter_weight_df <- scatter_plot_data(final_data)
scatter_weight_df <- scatter_weight_df[!is.na(scatter_weight_df$treatment),]


#plotting general results, colors by control / treatment

library(RColorBrewer)


# pathogen and individual infections
scatter_weight_df_patho <- scatter_weight_df[scatter_weight_df$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),]


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/individual_subset_scatterplot.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_line(linetype="dashed", 
            data = scatter_weight_df_patho[scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# with control not dotted.
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/individual_subset_scatterplot_control_visual1.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# DC3000 and co-infections
scatter_weight_df_patho <- scatter_weight_df[(grepl(x = scatter_weight_df$treatment, pattern = "DC3000")) | scatter_weight_df$treatment=="control",]

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/DC3000_subset_scatterplot.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_line(linetype="dashed", 
            data = scatter_weight_df_patho[scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# with control not dotted.
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/DC3000_subset_scatterplot_control_visual1.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()





# pathogen and co-infections
scatter_weight_df_patho <- scatter_weight_df[(grepl(x = scatter_weight_df$treatment, pattern = "patho")) | scatter_weight_df$treatment == "control" | scatter_weight_df$treatment == "Pathogen_(p4.C9)",]

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/coinfection_subset_scatterplot.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_line(linetype="dashed", 
            data = scatter_weight_df_patho[scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# with control not dotted.
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/coinfection_subset_scatterplot_control_visual1.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()

#### bayesian stats
final_data$delta_Green_pixels[final_data$day==7] <- final_data$Green_pixels[final_data$day==7]-final_data$Green_pixels[final_data$day==0]

final_data_delta <- final_data[final_data$day == 7,]

for (Treatment in unique(final_data_delta$Treatment)){
  #print(Treatment)
  treat_sd <- sd(x = final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment], na.rm = T)
  treat_mean <- mean(x = final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment], na.rm = T)
  max_pxls <- treat_mean + 2.5*treat_sd
  min_pxls <- treat_mean - 2.5*treat_sd
  if (any(final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] > max_pxls)){
    print(paste("too high", Treatment))
  }
  if (any(final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] < min_pxls)){
    print(paste("too low", Treatment))
  }
  final_data_delta <- final_data_delta[!final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] > max_pxls,]
  final_data_delta <- final_data_delta[!final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] < min_pxls,]
  
}

fit_sig_hits <- stan_glm(formula = delta_Green_pixels ~ treatment, data = final_data[final_data$day == 7,])

strain_pxls_7_0_stan_mod_df <- as.data.frame(fit_sig_hits$stan_summary)
strain_pxls_7_0_stan_mod_df$strain <- rownames(strain_pxls_7_0_stan_mod_df)
strain_pxls_7_0_stan_mod_df <- strain_pxls_7_0_stan_mod_df[!strain_pxls_7_0_stan_mod_df$strain %in% c("sigma","mean_PPD","log-posterior"),]

strain_pxls_7_0_stan_mod_df$strain <- as.factor(gsub(pattern = "treatment", replacement = "", x = as.character(strain_pxls_7_0_stan_mod_df$strain)))
strain_pxls_7_0_stan_mod_df$strain <- factor(strain_pxls_7_0_stan_mod_df$strain, levels = strain_pxls_7_0_stan_mod_df$strain[order(strain_pxls_7_0_stan_mod_df$`50%`,decreasing = T)])

grid = final_data_delta %>%
  data_grid(treatment)

fits = grid %>%
  add_fitted_draws(fit_sig_hits)

preds = grid %>%
  add_predicted_draws(fit_sig_hits)


treatment_order <- as.character(levels(strain_pxls_7_0_stan_mod_df$strain))
treatment_order[treatment_order=="(Intercept)"] <- "control"
final_data_delta$treatment <- factor(final_data_delta$treatment, levels = treatment_order)
preds$treatment <- factor(preds$treatment, levels = treatment_order)
fits$treatment <- factor(fits$treatment, levels = treatment_order)


# Strating with DC3000
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/DC3000_subset_bayesian.pdf", useDingbats = F)
ggplot(aes(x = strain, y=`50%`, ymin = `2.5%`, ymax = `97.5%`), 
       data = strain_pxls_7_0_stan_mod_df[(grepl(x = strain_pxls_7_0_stan_mod_df$strain, pattern = "DC3000")),])+
  geom_pointrange() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0, linetype="dashed") 
#! dev.off()


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/DC3000_subset_bayesian1.pdf", useDingbats = F)
final_data_delta[grepl(x = final_data_delta$treatment, pattern = "DC3000") | final_data_delta$treatment == "control",] %>%
  ggplot(aes(y = treatment, x = delta_Green_pixels)) +
  stat_interval(aes(x = .prediction), data = preds[grepl(x = preds$treatment, pattern = "DC3000") | preds$treatment == "control",]) +
  stat_pointinterval(aes(x = .value), data = fits[grepl(x = fits$treatment, pattern = "DC3000") | fits$treatment == "control",], .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()


# Now individual strains
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/individual_subset_bayesian.pdf", useDingbats = F)
ggplot(aes(x = strain, y=`50%`, ymin = `2.5%`, ymax = `97.5%`), 
       data = strain_pxls_7_0_stan_mod_df[strain_pxls_7_0_stan_mod_df$strain %in% c("GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),])+
  geom_pointrange() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0, linetype="dashed") 
#! dev.off()


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/individual_subset_bayesian1.pdf", useDingbats = F)
final_data_delta[final_data_delta$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),] %>%
  ggplot(aes(y = treatment, x = delta_Green_pixels)) +
  stat_interval(aes(x = .prediction), data = preds[preds$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),]) +
  stat_pointinterval(aes(x = .value), data = fits[fits$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),], .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()

# Now co-infections
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/coinfection_subset_bayesian.pdf", useDingbats = F)
ggplot(aes(x = strain, y=`50%`, ymin = `2.5%`, ymax = `97.5%`), 
       data = strain_pxls_7_0_stan_mod_df[(grepl(x = strain_pxls_7_0_stan_mod_df$strain, pattern = "patho")) | strain_pxls_7_0_stan_mod_df$strain == "control" | strain_pxls_7_0_stan_mod_df$strain == "Pathogen_(p4.C9)",])+
  geom_pointrange() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0, linetype="dashed") 
#! dev.off()

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_Ham_PCR_20210318/coinfection_subset_bayesian1.pdf", useDingbats = F)
final_data_delta[(grepl(x = final_data_delta$treatment, pattern = "patho")) | final_data_delta$treatment == "control" | final_data_delta$treatment == "Pathogen_(p4.C9)",] %>%
  ggplot(aes(y = treatment, x = delta_Green_pixels)) +
  stat_interval(aes(x = .prediction), data = preds[(grepl(x = preds$treatment, pattern = "patho")) | preds$treatment == "control" | preds$treatment == "Pathogen_(p4.C9)",]) +
  stat_pointinterval(aes(x = .value), data = fits[(grepl(x = fits$treatment, pattern = "patho")) | fits$treatment == "control" | fits$treatment == "Pathogen_(p4.C9)",], .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()























#### second exp

results <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/results_20210325.csv")
results$Position <- paste(results$Plate, results$Pot, sep = "_")

metadata <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/plan/second_exp_G2_and_Haim/G2_w_Haim_eight_ranodmized_final.csv")
metadata$Position <- paste("p", metadata$Position, sep = "")



final_data <- merge(x = results, y = metadata, by = "Position")
final_data$Green_pixels
final_data$Treatment <- as.character(final_data$Treatment)
final_data$Treatment[final_data$Treatment==" Control"] <- " control"
final_data$Treatment <- gsub(pattern = " ", replacement = "", x = final_data$Treatment)
final_data$Treatment <- as.factor(final_data$Treatment)

final_data_Haim <- final_data[final_data$Treatment %in% c("control", "p1.E1", "p12.A8" ,"p11.A2", "p25.C10", "p27.F3", "p7.A6", "p7.A9", "p8.H5", "Pathogen_alternative_(DC3000)"),]
#! write.table(x = final_data_Haim, file = "~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/Haim_results.csv", row.names = F)

final_data$delta_Green_pixels[final_data$day==7] <- final_data$Green_pixels[final_data$day==7]-final_data$Green_pixels[final_data$day==0]


ggplot(data = final_data[final_data$day==7,], aes(x = Treatment, y = delta_Green_pixels))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_boxplot()




##scatter plot the data
final_data$treatment <- final_data$Treatment
final_data$dpi <- final_data$day
final_data$pxls <- final_data$Green_pixels

source("~/ownCloud/My papers/bacterial_GWAS_paper/scripts/data_to_scatter_plot.R")
scatter_weight_df <- scatter_plot_data(final_data)
scatter_weight_df <- scatter_weight_df[!is.na(scatter_weight_df$treatment),]


#plotting general results, colors by control / treatment

library(RColorBrewer)

# Haim subset
scatter_weight_df_patho <- scatter_weight_df[scatter_weight_df$treatment %in% c("control", "p1.E1", "p12.A8" ,"p11.A2", "p25.C10", "p27.F3", "p7.A6", "p7.A9", "p8.H5"),]

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/Haim_subset_scatterplot.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_line(linetype="dashed", 
            data = scatter_weight_df_patho[scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# with control not dotted.
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/Haim_subset_scatterplot_control_visual1.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()



# pathogen and individual infections
scatter_weight_df_patho <- scatter_weight_df[scatter_weight_df$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),]

scatter_weight_df_patho$only_sig_ko <- scatter_weight_df_patho$treatment
scatter_weight_df_patho$treatment <- factor(scatter_weight_df_patho$treatment, levels = c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"))



#colors_plot <- brewer.pal(12, "Set3")
#colors_plot <- c(colors_plot,"grey")
#colors_plot <- c("#000000","#fdc086","#a6611a", "#386cb0", "#fdc086", "grey50")
colors_plot <- c("#000000", "#a6611a","#386cb0","#7cae00" ,"#af0dd6","#edf903")


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_scatterplot.pdf", useDingbats = F)
ggplot(data = scatter_weight_df_patho, aes(x = dpi, y = pxls_average, group = treatment, color = treatment)) +
  #geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),]) +
  #geom_point(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),]) +
  #geom_line(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),]) +
  #geom_point(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),])+
  geom_line() +
  geom_point() +
  scale_color_manual(values = colors_plot) +
  theme_classic() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") 
#! dev.off()

colors_plot <- c("#000000", "#a6611a","#386cb0","#7cae00" ,"#af0dd6","#edf903")

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_scatterplot_visual1.pdf", useDingbats = F)
ggplot(data = scatter_weight_df_patho, aes(x = dpi, y = pxls_average, group = treatment, color = treatment)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_line(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000"),]) +
  geom_point(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000"),])+
  
  scale_color_manual(values = colors_plot) +
  theme_classic() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") 
#! dev.off()

scatter_weight_df_patho$only_sig_ko <- factor(scatter_weight_df_patho$only_sig_ko)
levels(scatter_weight_df_patho$only_sig_ko) <- c("", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "", "")

colors_plot <- c("#000000", "grey","grey" ,"grey","#386cb0")

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/individual_subset_scatterplot_visual2.pdf", useDingbats = F, width = 10, height = 10)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line() +
  geom_point() +
  #geom_line(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00003884_12+DC3000"),]) +
  #geom_point(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00003884_12+DC3000"),])+
  scale_color_manual(values = colors_plot) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  #scale_color_brewer(palette = "Set1") +
  directlabels::geom_dl(aes(label = only_sig_ko), method = list(directlabels::dl.combine("last.points")), cex = 0.8) 
#! dev.off()


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/individual_subset_scatterplot.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_line(linetype="dashed", 
            data = scatter_weight_df_patho[scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# with control not dotted.
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/individual_subset_scatterplot_control_visual1.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# DC3000 and co-infections
scatter_weight_df_patho <- scatter_weight_df[(grepl(x = scatter_weight_df$treatment, pattern = "DC3000")) | scatter_weight_df$treatment=="control",]

scatter_weight_df_patho$only_sig_ko <- scatter_weight_df_patho$treatment
scatter_weight_df_patho$treatment <- factor(scatter_weight_df_patho$treatment, levels = c("control","Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00000450_7+DC3000", "GC00000032_87-GC50_54+DC3000", "GC00003884_12+DC3000"))



#colors_plot <- brewer.pal(12, "Set3")
#colors_plot <- c(colors_plot,"grey")
#colors_plot <- c("#000000","#fdc086","#a6611a", "#386cb0", "#fdc086", "grey50")
colors_plot <- c("#000000", "#a6611a","#386cb0","#7cae00" ,"#af0dd6","#edf903")


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_scatterplot.pdf", useDingbats = F)
ggplot(data = scatter_weight_df_patho, aes(x = dpi, y = pxls_average, group = treatment, color = treatment)) +
  #geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),]) +
  #geom_point(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),]) +
  #geom_line(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),]) +
  #geom_point(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "p5.F2"),])+
  geom_line() +
  geom_point() +
  scale_color_manual(values = colors_plot) +
  theme_classic() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") 
#! dev.off()

colors_plot <- c("#000000", "#a6611a","#386cb0","#7cae00" ,"#af0dd6","#edf903")

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_scatterplot_visual1.pdf", useDingbats = F)
ggplot(data = scatter_weight_df_patho, aes(x = dpi, y = pxls_average, group = treatment, color = treatment)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_line(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000"),]) +
  geom_point(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000"),])+

  scale_color_manual(values = colors_plot) +
  theme_classic() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") 
#! dev.off()

scatter_weight_df_patho$only_sig_ko <- factor(scatter_weight_df_patho$only_sig_ko)
levels(scatter_weight_df_patho$only_sig_ko) <- c("", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "", "")
scatter_weight_df_patho$treatment <- factor(scatter_weight_df_patho$treatment, levels = c("control","Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00000450_7+DC3000", "GC00000032_87-GC50_54+DC3000", "GC00003884_12+DC3000"))


colors_plot <- c("#000000", "#f03b20","#63e5cc","grey" ,"grey","grey")

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_scatterplot_visual3.pdf", useDingbats = F, width = 10, height = 10)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(linetype = "dashed", size =0.7) +
  geom_point(size = 2) +
  geom_line(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00003884_12+DC3000"),], size = 0.7) +
  geom_point(data = scatter_weight_df_patho[scatter_weight_df_patho$treatment %in% c("control", "Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00003884_12+DC3000"),], size = 2)+
  scale_color_manual(values = colors_plot) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  #scale_color_brewer(palette = "Set1") +
  directlabels::geom_dl(aes(label = only_sig_ko), method = list(directlabels::dl.combine("last.points")), cex = 0.8) 
#! dev.off()


# with control not dotted.
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_scatterplot_control_visual1.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# pathogen and co-infections
scatter_weight_df_patho <- scatter_weight_df[(grepl(x = scatter_weight_df$treatment, pattern = "patho")) | scatter_weight_df$treatment == "control" | scatter_weight_df$treatment == "Pathogen_(p4.C9)",]

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/coinfection_subset_scatterplot.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho[!scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_line(linetype="dashed", 
            data = scatter_weight_df_patho[scatter_weight_df_patho$treatment%in%c("control"),]) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()


# with control not dotted.
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/coinfection_subset_scatterplot_control_visual1.pdf", useDingbats = F)
ggplot(aes(x=dpi, y=pxls_average, group=treatment, color=treatment), data = scatter_weight_df_patho) +
  geom_line(data = scatter_weight_df_patho) +
  geom_point()+
  #scale_color_manual(values = colors_scatter_by_controls) +
  theme_bw() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  scale_color_brewer(palette = "Set1")
#! dev.off()




#### bayesian stats
final_data$delta_Green_pixels[final_data$day==7] <- final_data$Green_pixels[final_data$day==7]-final_data$Green_pixels[final_data$day==0]
final_data$treatment <- final_data$Treatment

final_data_delta <- final_data[final_data$day == 7,]

#for (Treatment in unique(final_data_delta$Treatment)){
#  #print(Treatment)
#  treat_sd <- sd(x = final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment], na.rm = T)
#  treat_mean <- mean(x = final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment], na.rm = T)
#  max_pxls <- treat_mean + 2.5*treat_sd
#  min_pxls <- treat_mean - 2.5*treat_sd
#  if (any(final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] > max_pxls)){
#    print(paste("too high", Treatment))
#  }
#  if (any(final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] < min_pxls)){
#    print(paste("too low", Treatment))
#  }
#  final_data_delta <- final_data_delta[!final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] > max_pxls,]
#  final_data_delta <- final_data_delta[!final_data_delta$delta_Green_pixels[final_data_delta$Treatment == Treatment] < min_pxls,]
#  
#}



fit_sig_hits <- stan_glm(formula = delta_Green_pixels ~ treatment, data = final_data_delta)

strain_pxls_7_0_stan_mod_df <- as.data.frame(fit_sig_hits$stan_summary)
strain_pxls_7_0_stan_mod_df$strain <- rownames(strain_pxls_7_0_stan_mod_df)
strain_pxls_7_0_stan_mod_df <- strain_pxls_7_0_stan_mod_df[!strain_pxls_7_0_stan_mod_df$strain %in% c("sigma","mean_PPD","log-posterior"),]

strain_pxls_7_0_stan_mod_df$strain <- as.factor(gsub(pattern = "treatment", replacement = "", x = as.character(strain_pxls_7_0_stan_mod_df$strain)))
strain_pxls_7_0_stan_mod_df$strain <- factor(strain_pxls_7_0_stan_mod_df$strain, levels = strain_pxls_7_0_stan_mod_df$strain[order(strain_pxls_7_0_stan_mod_df$`50%`,decreasing = T)])

grid = final_data_delta %>%
  data_grid(treatment)

fits = grid %>%
  add_fitted_draws(fit_sig_hits)

preds = grid %>%
  add_predicted_draws(fit_sig_hits)


treatment_order <- as.character(levels(strain_pxls_7_0_stan_mod_df$strain))
treatment_order[treatment_order=="(Intercept)"] <- "control"
final_data_delta$treatment <- factor(final_data_delta$treatment, levels = treatment_order)
preds$treatment <- factor(preds$treatment, levels = treatment_order)
fits$treatment <- factor(fits$treatment, levels = treatment_order)


# Strating with DC3000
strain_pxls_7_0_stan_mod_df$sig_diff <- strain_pxls_7_0_stan_mod_df$`97.5%`<0

strain_pxls_7_0_stan_mod_df_DC3000 <- strain_pxls_7_0_stan_mod_df[(grepl(x = strain_pxls_7_0_stan_mod_df$strain, pattern = "DC3000")),]
strain_pxls_7_0_stan_mod_df_DC3000$strain <- factor(x = strain_pxls_7_0_stan_mod_df_DC3000$strain, 
       levels = c("Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00000450_7+DC3000", "GC00000032_87-GC50_54+DC3000", "GC00003884_12+DC3000"))

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_bayesian.pdf", width = 6,height = 10, useDingbats = F)
#! pdf("/Users/oshalev/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/FigureS5/DC3000_subset_bayesian_scaled.pdf", width = 6, height = 10, useDingbats = F)
ggplot(data = strain_pxls_7_0_stan_mod_df_DC3000, aes(x = strain, y = `50%`, ymin = `2.5%`, ymax = `97.5%`, color = sig_diff))+
  geom_pointrange(aes(group = strain), position = position_dodge(width=1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("black", "#d95f02")) +
  coord_cartesian(ylim = c(-80000,5000))
#! dev.off()


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_bayesian.pdf", useDingbats = F)
ggplot(aes(x = strain, y=`50%`, ymin = `2.5%`, ymax = `97.5%`), 
       data = strain_pxls_7_0_stan_mod_df[(grepl(x = strain_pxls_7_0_stan_mod_df$strain, pattern = "DC3000")),])+
  geom_pointrange() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0, linetype="dashed") 
#! dev.off()


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/DC3000_subset_bayesian1.pdf", useDingbats = F)
final_data_delta[grepl(x = final_data_delta$treatment, pattern = "DC3000") | final_data_delta$treatment == "control",] %>%
  ggplot(aes(y = treatment, x = delta_Green_pixels)) +
  stat_interval(aes(x = .prediction), data = preds[grepl(x = preds$treatment, pattern = "DC3000") | preds$treatment == "control",]) +
  stat_pointinterval(aes(x = .value), data = fits[grepl(x = fits$treatment, pattern = "DC3000") | fits$treatment == "control",], .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()

colors_plot <- c("#000000", "#a6611a","#386cb0","#fdc086" ,"grey80")


# now comparison to pathogen, rather than contorl

final_data_DC3000_lm <- final_data_delta[final_data_delta$treatment %in% c("Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00000450_7+DC3000", "GC00000032_87-GC50_54+DC3000", "GC00003884_12+DC3000"),]
final_data_DC3000_lm$treatment <- factor(x = final_data_DC3000_lm$treatment, levels = c("Pathogen_alternative_(DC3000)", "Strain_17+DC3000", "GC00000450_7+DC3000", "GC00000032_87-GC50_54+DC3000", "GC00003884_12+DC3000"))

lm_dc30000 <- stan_glm(formula = delta_Green_pixels ~ treatment, data = final_data_DC3000_lm)

lm_dc30000$stan_summary[,c(4,7,10)]


strain_pxls_7_0_stan_mod_df <- as.data.frame(fit_sig_hits$stan_summary)
strain_pxls_7_0_stan_mod_df$strain <- rownames(strain_pxls_7_0_stan_mod_df)
strain_pxls_7_0_stan_mod_df <- strain_pxls_7_0_stan_mod_df[!strain_pxls_7_0_stan_mod_df$strain %in% c("sigma","mean_PPD","log-posterior"),]

strain_pxls_7_0_stan_mod_df$strain <- as.factor(gsub(pattern = "treatment", replacement = "", x = as.character(strain_pxls_7_0_stan_mod_df$strain)))
strain_pxls_7_0_stan_mod_df$strain <- factor(strain_pxls_7_0_stan_mod_df$strain, levels = strain_pxls_7_0_stan_mod_df$strain[order(strain_pxls_7_0_stan_mod_df$`50%`,decreasing = T)])




# Now individual strains
strain_pxls_7_0_stan_mod_df_individual <- strain_pxls_7_0_stan_mod_df[strain_pxls_7_0_stan_mod_df$strain %in% c("GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),]
strain_pxls_7_0_stan_mod_df_individual$strain <- factor(x = strain_pxls_7_0_stan_mod_df_individual$strain, 
                                                    levels = c("Strain_17", "GC00000450_7", "GC00000032_87-GC50_54", "GC00003884_12"))

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/individual_subset_bayesian.pdf", width = 6,height = 10, useDingbats = F)
ggplot(data = strain_pxls_7_0_stan_mod_df_individual, aes(x = strain, y = `50%`, ymin = `2.5%`, ymax = `97.5%`, color = sig_diff))+
  geom_pointrange(aes(group = strain), position = position_dodge(width=1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("black", "#d95f02"))
#! dev.off()

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/individual_subset_bayesian.pdf", useDingbats = F)
ggplot(aes(x = strain, y=`50%`, ymin = `2.5%`, ymax = `97.5%`), 
       data = strain_pxls_7_0_stan_mod_df[strain_pxls_7_0_stan_mod_df$strain %in% c("GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),])+
  geom_pointrange() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0, linetype="dashed") 
#! dev.off()


#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/individual_subset_bayesian1.pdf", useDingbats = F)
final_data_delta[final_data_delta$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),] %>%
  ggplot(aes(y = treatment, x = delta_Green_pixels)) +
  stat_interval(aes(x = .prediction), data = preds[preds$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),]) +
  stat_pointinterval(aes(x = .value), data = fits[fits$treatment %in% c("control", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Strain_17"),], .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()

# Now co-infections
#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/coinfection_subset_bayesian.pdf", useDingbats = F)
ggplot(aes(x = strain, y=`50%`, ymin = `2.5%`, ymax = `97.5%`), 
       data = strain_pxls_7_0_stan_mod_df[(grepl(x = strain_pxls_7_0_stan_mod_df$strain, pattern = "patho")) | strain_pxls_7_0_stan_mod_df$strain == "control" | strain_pxls_7_0_stan_mod_df$strain == "Pathogen_(p4.C9)",])+
  geom_pointrange() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0, linetype="dashed") 
#! dev.off()

#! pdf("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/LAST_EXPs_20201/results_G2_and_haim_20210325/coinfection_subset_bayesian1.pdf", useDingbats = F)
final_data_delta[(grepl(x = final_data_delta$treatment, pattern = "patho")) | final_data_delta$treatment == "control" | final_data_delta$treatment == "Pathogen_(p4.C9)",] %>%
  ggplot(aes(y = treatment, x = delta_Green_pixels)) +
  stat_interval(aes(x = .prediction), data = preds[(grepl(x = preds$treatment, pattern = "patho")) | preds$treatment == "control" | preds$treatment == "Pathogen_(p4.C9)",]) +
  stat_pointinterval(aes(x = .value), data = fits[(grepl(x = fits$treatment, pattern = "patho")) | fits$treatment == "control" | fits$treatment == "Pathogen_(p4.C9)",], .width = c(.66, .95), position = position_nudge(y = -0.3)) +
  geom_point(aes(alpha=0.2)) +
  scale_color_brewer() +
  theme_classic()
#! dev.off()
