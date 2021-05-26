library(ggplot2)
library(RColorBrewer)
library(rstanarm)
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
library(directlabels)



results_summary <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._experiment/results/summary/processed_area_in_pixels_summary.csv")
meatadata <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._experiment/randomized_13-plates_bacterial_K.O._October_2020_AFTER_INFECTIONS.csv")
treatments_meatadata <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._experiment/Treatments_metadata.csv")

full_results_summary <- merge(x = results_summary, y = meatadata, by = "position")
full_results_summary <- merge(x = full_results_summary, y = treatments_meatadata, by = "Treatment")

full_results_summary$Treatment <- as.factor(full_results_summary$Treatment)

scatter_plot_df <- data.frame("dpi"=NULL, "med_pxls"=NULL,"Treatment"=NULL, "Treatment_details"=NULL)
for (Treatment in unique(full_results_summary$Treatment)){
  for (dpi in unique(full_results_summary$day)){
    med_pxls <- median(full_results_summary$Green_pixels[full_results_summary$Treatment==Treatment & full_results_summary$day==dpi])
    Treatment_details <- unique(full_results_summary$Treatment_details[full_results_summary$Treatment == Treatment])
    temp_df_scatter <- data.frame(dpi,med_pxls,Treatment, Treatment_details)
    scatter_plot_df <- rbind(scatter_plot_df, temp_df_scatter)
  }
}




scatter_plot_df <- scatter_plot_df[scatter_plot_df$Treatment_details!="WT (p5.F2)",]

Treatment_details_levels <- c("Control", "Pathogen only", "WT only (p5.F2)", "strain_17",
                                               "GC00000450_7", "GC00000032_87-GC50_54", "GC00003884_12",
                                               "GC00000009_145", "GC00003215_1", "GC00004535_r1_1", "GC00005419_3",
                                               "GC00005987_r1_1")

scatter_plot_df$Treatment_details <- factor(scatter_plot_df$Treatment_details, levels = Treatment_details_levels)

scatter_plot_df$treat_by_group <- scatter_plot_df$Treatment_details
levels(scatter_plot_df$treat_by_group) <- c("Control", "Pathogen", "p5.F2",
                                            "Pathogen + p5.F2", rep("Pathogen + p5.F2 K.O.", 8))

scatter_plot_df$only_sig_ko <- scatter_plot_df$Treatment_details
levels(scatter_plot_df$only_sig_ko) <- c("", "", "",
                                            "", "GC00000450_7", "GC00000032_87 / \n GC00000050_54", "GC00003884_12",
                                         "", "", "", "", "")


#colors_plot <- brewer.pal(12, "Set3")
#colors_plot <- c(colors_plot,"grey")
#colors_plot <- c("#000000","#fdc086","#a6611a", "#386cb0", "#fdc086", "grey50")
colors_plot <- c("#000000", "#a6611a","#386cb0","#fdc086" ,"grey80")

scatter_plot_df$dpi <- as.factor(scatter_plot_df$dpi)
#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/Scatter_plot_KO_updated.pdf", useDingbats = F)
ggplot(data = scatter_plot_df, aes(x = dpi, y = med_pxls, group = Treatment, color = treat_by_group)) +
  geom_line() +
  geom_point()+
  geom_line(data = scatter_plot_df[scatter_plot_df$Treatment_details %in% c("Control", "Pathogen", "p5.F2",
                                                     "Pathogen + p5.F2"),]) +
  geom_point(data = scatter_plot_df[scatter_plot_df$Treatment_details %in% c("Control", "Pathogen", "p5.F2",
                                                                      "Pathogen + p5.F2"),])+
  scale_color_manual(values = colors_plot) +
  theme_classic() +
  ylab(label = "median pixels") +
  xlab(label = "dpi") +
  geom_dl(aes(label = only_sig_ko), method = list(dl.combine("last.points")), cex = 0.8) 
#! dev.off()









full_results_summary$Green_pixels_difference[full_results_summary$day == 7] <- full_results_summary$Green_pixels[full_results_summary$day == 7] - full_results_summary$Green_pixels[full_results_summary$day == 1]
full_results_summary_difference <- full_results_summary[full_results_summary$day == 7,]


unique(full_results_summary_difference$Treatment_details)
full_results_summary_difference <- full_results_summary_difference[!(full_results_summary_difference$Treatment_details %in% c("WT only (p5.F2)", "WT (p5.F2)")),]
full_results_summary_difference$Treatment_details <- factor(full_results_summary_difference$Treatment_details, levels = c("Control", "Pathogen only", "strain_17",
                                                                     "GC00003215_1", "GC00004535_r1_1", "GC00000009_145",
                                                                     "GC00005987_r1_1", "GC00005419_3", "GC00000032_87-GC50_54",
                                                                     "GC00000450_7", "GC00003884_12"))

fit_sig_hits <- stan_glm(formula = Green_pixels_difference ~ Treatment_details , data = full_results_summary_difference, seed = 12345)
fit_sig_hits_df <- as.data.frame(fit_sig_hits$stan_summary)
fit_sig_hits_df$treatment <- rownames(fit_sig_hits_df)
fit_sig_hits_df <- fit_sig_hits_df[!fit_sig_hits_df$treatment %in% c("sigma", "log-posterior", "mean_PPD", "(Intercept)"), c('treatment','2.5%', '50%', '97.5%')]
fit_sig_hits_df$treatment <- factor(fit_sig_hits_df$treatment)

levels(fit_sig_hits_df$treatment) <- gsub('Treatment_details', '', levels(fit_sig_hits_df$treatment))
fit_sig_hits_df$treatment <- factor(fit_sig_hits_df$treatment, levels = c("Pathogen only", "strain_17", "GC00000450_7", 
                                             "GC00000032_87-GC50_54", "GC00003884_12", "GC00003215_1",
                                             "GC00004535_r1_1", "GC00000009_145",
                                             "GC00005987_r1_1", "GC00005419_3"))

levels(fit_sig_hits_df$treatment) <- c("Pathogen", "Pathogen + p5.F2", "Pathogen + GC00000450_7", 
                                       "Pathogen + GC00000032_87 // GC00000050_54", "Pathogen + GC00003884_12", "Pathogen + GC00003215_1",
                                       "Pathogen + GC00004535_r1_1", "Pathogen + GC00000009_145",
                                       "Pathogen + GC00005987_r1_1", "Pathogen + GC00005419_3")



fit_sig_hits_df <- fit_sig_hits_df[fit_sig_hits_df$treatment %in% c("Pathogen", "Pathogen + p5.F2", "Pathogen + GC00000450_7", 
                                 "Pathogen + GC00000032_87 // GC00000050_54", "Pathogen + GC00003884_12"),]

# Plotting the effect of each treatment, all in comparison to Control.
# At the end i am producing 95% CI for all effects.

fit_sig_hits_df$sig_diff <- fit_sig_hits_df$`97.5%`<0

#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/KO_glm_to_control_only_KO_of_interest.pdf", width = 6,height = 10, useDingbats = F)
ggplot(data = fit_sig_hits_df, aes(x = treatment, y = `50%`, ymin = `2.5%`, ymax = `97.5%`, color = sig_diff))+
  geom_pointrange(aes(group = treatment), position = position_dodge(width=1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("black", "#d95f02"))
#! dev.off()




#! pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/Scatter_plot_KO.pdf", useDingbats = F)
ggplot(data = full_results_summary_difference, aes(x = Treatment_details, y = Green_pixels_difference, group = Treatment_details, color = Treatment_details)) +
  geom_boxplot() +
  geom_point()+
  theme_classic() +
  ylab(label = "median pixels") +
  xlab(label = "dpi")
#! dev.off()



# now in comparison to the pathogen only, rather than control

pathogen_summary_difference <- full_results_summary_difference[full_results_summary_difference$Treatment_details!="Control",]

fit_sig_hits <- stan_glm(formula = Green_pixels_difference ~ Treatment_details , data = pathogen_summary_difference, seed = 12345)
fit_sig_hits_df <- as.data.frame(fit_sig_hits$stan_summary)
fit_sig_hits_df$treatment <- rownames(fit_sig_hits_df)
fit_sig_hits_df <- fit_sig_hits_df[!fit_sig_hits_df$treatment %in% c("sigma", "log-posterior", "mean_PPD", "(Intercept)"), c('treatment','2.5%', '50%', '97.5%')]
fit_sig_hits_df$treatment <- factor(fit_sig_hits_df$treatment)
fit_sig_hits_df
