library(dplyr)
library(ggplot2)
library(growthcurver)
library(rstanarm)

iron_exp_results_p1 <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._iron_exp/results_plates_LB/plate1/iron_exp_results_p1_20201125.csv")
iron_exp_results_p2 <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._iron_exp/results_plates_LB/plate2/iron_exp_results_p2_20201126.csv")
iron_exp_results_p3 <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._iron_exp/results_plates_LB/plate3/iron_exp_plate3_20201127.csv")

iron_exp_results_p3 <- iron_exp_results_p3[,!(colnames(iron_exp_results_p3) %in% "H5")]

# transform the data
iron_exp_results_p1_t <- (tidyr::gather(iron_exp_results_p1))
iron_exp_results_p1_t <- iron_exp_results_p1_t[iron_exp_results_p1_t$value!="",]
iron_exp_results_p1_t$time <- as.numeric(iron_exp_results_p1_t$value[iron_exp_results_p1_t$key=="Time..s."])
iron_exp_results_p1_t <- iron_exp_results_p1_t[iron_exp_results_p1_t$key!="Time..s.",]
iron_exp_results_p1_t$plate <- "p1"

iron_exp_results_p2_t <- (tidyr::gather(iron_exp_results_p2))
iron_exp_results_p2_t <- iron_exp_results_p2_t[iron_exp_results_p2_t$value!="",]
iron_exp_results_p2_t$time <- as.numeric(iron_exp_results_p2_t$value[iron_exp_results_p2_t$key=="Time..s."])
iron_exp_results_p2_t <- iron_exp_results_p2_t[iron_exp_results_p2_t$key!="Time..s.",]
iron_exp_results_p2_t$plate <- "p2"

iron_exp_results_p3_t <- (tidyr::gather(iron_exp_results_p3))
iron_exp_results_p3_t <- iron_exp_results_p3_t[iron_exp_results_p3_t$value!="",]
iron_exp_results_p3_t$time <- as.numeric(iron_exp_results_p3_t$value[iron_exp_results_p3_t$key=="Time..s."])
iron_exp_results_p3_t <- iron_exp_results_p3_t[iron_exp_results_p3_t$key!="Time..s.",]
iron_exp_results_p3_t$plate <- "p3"


iron_exp_t <- rbind(iron_exp_results_p1_t,iron_exp_results_p2_t)

iron_exp_metadata <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._iron_exp/plate1_ranodmized_edited.csv")
iron_exp_metadata_p3 <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/K.O._iron_exp/plate3_ranodmized_high_concentrations_edited.csv")

iron_exp_results_df <- merge(x = iron_exp_metadata, y = iron_exp_t, by.x = "well", by.y = "key")
iron_exp_results_df_p3 <- merge(x = iron_exp_metadata_p3, y = iron_exp_results_p3_t, by.x = "well", by.y = "key")

iron_exp_results_df <- rbind(iron_exp_results_df, iron_exp_results_df_p3)

iron_exp_results_df$value <- as.numeric(iron_exp_results_df$value)


iron_exp_results_df$condition <- factor(iron_exp_results_df$condition, levels = c("0_uM", "50_uM", "100_uM", "150_uM", "200_uM", 
                                                                                  "250_uM", "300_uM", "400_uM", "500_uM", "600_uM", "700_uM"))

iron_exp_results_df$treatment <- factor(x = iron_exp_results_df$treatment, 
       levels = c("Strain_17", "WT", "GC00000450_7", "GC00000032_87-GC50_54", "GC00003884_12", "Pathogen_(p4.C9)", ""))


# subset to first batch only (up to 300_uM)
iron_exp_results_df_subset <- iron_exp_results_df[iron_exp_results_df$condition %in% c("0_uM", "50_uM", "100_uM","150_uM", "200_uM", "250_uM", "300_uM"),]

#!pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/iron_exp_growth_curve_no_GC00003884_12_subset.pdf"
#!   , useDingbats = F, width = 15, height = 10)
ggplot(data = iron_exp_results_df_subset[!(iron_exp_results_df_subset$treatment %in% c("", "WT", "GC00003884_12")) &
                                    iron_exp_results_df_subset$time<35000,],
       aes(x = time, y = log(value), color = treatment, group = treatment, fill = treatment))+
  geom_point(alpha=0.3) +
  geom_smooth() +
  stat_smooth(geom='line', alpha=0.5, se=FALSE) +
  scale_color_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3", "#d95f02")) +
  scale_fill_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3", "#d95f02")) +
  facet_wrap(. ~ condition) +
  theme_classic()
#!dev.off()

#!pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/iron_exp_growth_curve_no_GC00003884_12_GC00000032_87-GC50_54_subset.pdf"
#!    , useDingbats = F, width = 15, height = 10)
ggplot(data = iron_exp_results_df_subset[!(iron_exp_results_df_subset$treatment %in% c("", "WT", "GC00003884_12", "GC00000032_87-GC50_54")) &
                                    iron_exp_results_df_subset$time<35000,],
       aes(x = time, y = log(value), color = treatment, group = treatment, fill = treatment))+
  geom_point(alpha=0.3) +
  geom_smooth() +
  stat_smooth(geom='line', alpha=0.5, se=FALSE) +
  scale_color_manual(values = c("#b2df8a", "#1f78b4","#d95f02")) +
  scale_fill_manual(values = c("#b2df8a", "#1f78b4","#d95f02")) +
  facet_wrap(. ~ condition) +
  theme_classic()
#!dev.off()

#!pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/iron_exp_growth_curve_no_GC00003884_12_GC00000450_7_subset.pdf"
#!    , useDingbats = F, width = 15, height = 10)
ggplot(data = iron_exp_results_df_subset[!(iron_exp_results_df_subset$treatment %in% c("", "WT", "GC00003884_12", "GC00000450_7")) &
                                    iron_exp_results_df_subset$time<35000,],
       aes(x = time, y = log(value), color = treatment, group = treatment, fill = treatment))+
  geom_point(alpha=0.3) +
  geom_smooth() +
  stat_smooth(geom='line', alpha=0.5, se=FALSE) +
  scale_color_manual(values = c("#b2df8a", "#a6cee3","#d95f02")) +
  scale_fill_manual(values = c("#b2df8a", "#a6cee3","#d95f02")) +
  facet_wrap(. ~ condition) +
  theme_classic()
#!dev.off()



#!pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/iron_exp_growth_curve_all_subset.pdf"
#!   , useDingbats = F, width = 15, height = 10)
ggplot(data = iron_exp_results_df_subset[!(iron_exp_results_df_subset$treatment %in% c("", "WT")) &
                                    iron_exp_results_df_subset$time<35000,],
       aes(x = time, y = log(value), color = treatment, group = treatment, fill = treatment))+
  geom_point(alpha=0.3) +
  geom_smooth() +
  stat_smooth(geom='line', alpha=0.5, se=FALSE) +
  scale_color_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3", "#c6dbef" ,"#d95f02")) +
  scale_fill_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3","#c6dbef" ,"#d95f02")) +
  facet_wrap(. ~ condition) +
  theme_classic()
#!dev.off()


#!pdf("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/iron_exp_growth_curve_WT_s17_pathogen_subset.pdf"
#!      , useDingbats = F, width = 15, height = 10)
ggplot(data = iron_exp_results_df_subset[!(iron_exp_results_df_subset$treatment %in% c("", "GC00003884_12", "GC00000450_7", "GC00000032_87-GC50_54")) &
                                    iron_exp_results_df_subset$time<35000,],
       aes(x = time, y = log(value), color = treatment, group = treatment, fill = treatment))+
  geom_point(alpha=0.3) +
  geom_smooth() +
  stat_smooth(geom='line', alpha=0.5, se=FALSE) +
  scale_color_manual(values = c("#b2df8a", "#005a32", "#d95f02")) +
  scale_fill_manual(values = c("#b2df8a", "#005a32", "#d95f02")) +
  facet_wrap(. ~ condition) +
  theme_classic()
#!dev.off()





########growth curver part

growthcurve_per_plate <- function(data=NULL, trim=35000){
  data <- data[data$Time..s.!="",]
  colnames(data)[colnames(data)=="Time..s."] <- "time"
  
  for (i in 1:length(data)){
    data[,i] <- as.numeric(as.character(data[,i]))
  }
  
  data_growth <- SummarizeGrowthByPlate(data, t_trim = trim, bg_correct = "min")
  return(data_growth)
}


conclude_growth_by_conditions <- function(trim_time = 35000, conditions_subset = c(0,100,150,200,250,300,50,400,500,600,700), treatments_to_exclude = c("","GC00003884_12", "WT")){
  iron_exp_results_p1_growth <- growthcurve_per_plate(data=iron_exp_results_p1, trim = trim_time)
  iron_exp_results_p2_growth <- growthcurve_per_plate(data=iron_exp_results_p2, trim = trim_time)
  iron_exp_results_p3_growth <- growthcurve_per_plate(data=iron_exp_results_p3, trim = trim_time)
  
  iron_exp_results_p1_growth$plate <- "p1"
  iron_exp_results_p2_growth$plate <- "p2"
  iron_exp_results_p3_growth$plate <- "p3"
  
  iron_exp_results_growth <- rbind(iron_exp_results_p1_growth, iron_exp_results_p2_growth)
  iron_exp_results_growth_df <- merge(x = iron_exp_metadata, y = iron_exp_results_growth , by.x = "well", by.y = "sample")
  iron_exp_results_growth_df_p3 <- merge(x = iron_exp_metadata_p3, y = iron_exp_results_p3_growth , by.x = "well", by.y = "sample")
  
  iron_exp_results_growth_df <- rbind(iron_exp_results_growth_df,iron_exp_results_growth_df_p3)
  levels(iron_exp_results_growth_df$condition) <- c(0,100,150,200,250,300,50,400,500,600,700)
  iron_exp_results_growth_df$condition <- as.numeric(as.character(iron_exp_results_growth_df$condition))
  
  file_dir <- "~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Figure4/growth_properties/"
  conditions_filename <- paste(c(min(conditions_subset),max(conditions_subset)), collapse = "_")
  iron_exp_results_growth_df <- iron_exp_results_growth_df[!(iron_exp_results_growth_df$treatment %in% treatments_to_exclude) & 
                               iron_exp_results_growth_df$condition %in% conditions_subset,]
  # plot the regression plot

  
  reg_plot <- ggplot(data = iron_exp_results_growth_df, 
         aes(x = condition, y = auc_l, color = treatment)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "glm", aes(fill = treatment), alpha = 0.2) +
    theme_classic()

  filename <- paste(paste(file_dir, "regression_plot", sep = ""), conditions_filename, sep = "_")
  pdf(filename , useDingbats = F, width = 10, height = 10)
  print(reg_plot)
  dev.off()
  
  iron_exp_results_growth_df$treatment <- factor(x = iron_exp_results_growth_df$treatment, levels = c("Strain_17", "GC00000032_87-GC50_54", "GC00000450_7", "GC00003884_12", "Pathogen_(p4.C9)", "WT", ""))
  fit <- stan_glm(data = iron_exp_results_growth_df, 
                  formula = auc_l ~ treatment*condition, iter = 15000)
  
  # first plot only the interactions of treatment by condition
  coeffcients_interaction_factors <- rownames(fit$stan_summary)[grepl(pattern = ":condition", x = rownames(fit$stan_summary))]
  fit_interactions <- as.data.frame(fit$stan_summary[coeffcients_interaction_factors,c("2.5%","50%","97.5%")])
  fit_interactions$condition <- rownames(fit_interactions)
  fit_interactions$condition <- gsub(pattern = "treatment", replacement = "", x = fit_interactions$condition)
  fit_interactions$condition <- gsub(pattern = ":condition", replacement = "", x = fit_interactions$condition)
  fit_interactions$sig_dif <- fit_interactions[,"97.5%"] < 0
  
  interactions_plot <- ggplot(data = fit_interactions, 
                              aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`, color = sig_dif))+
    geom_pointrange(position=position_dodge(width=1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values = c("black", "#d95f02"))
  
  
  filename <- paste(paste(file_dir, "interaction_plot", sep = ""), conditions_filename, sep = "_")
  pdf(filename , useDingbats = F)
  print(interactions_plot)
  dev.off()
  
  # now plot only the treatment factor in the model
  coeffcients_treatment <- rownames(fit$stan_summary)[grepl(pattern = "treatment", x = rownames(fit$stan_summary))]
  coeffcients_treatment <- coeffcients_treatment[!grepl(":condition",coeffcients_treatment)]
  fit_treatments <- as.data.frame(fit$stan_summary[coeffcients_treatment,c("2.5%","50%","97.5%")])
  fit_treatments$condition <- rownames(fit_treatments)
  fit_treatments$condition <- gsub(pattern = "treatment", replacement = "", x = fit_treatments$condition)
  fit_treatments$sig_dif <- fit_treatments[,"97.5%"] < 0
  
  treatment_plot <- ggplot(data = fit_treatments, 
                           aes(x=condition,y=`50%`, ymin = `2.5%`, ymax = `97.5%`, color = sig_dif))+
    geom_pointrange(position=position_dodge(width=1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values = c("black", "#d95f02"))

  
  filename <- paste(paste(file_dir, "treatment_coefficient_plot", sep = ""), conditions_filename, sep = "_")
  pdf(filename , useDingbats = F)
  print(treatment_plot)
  dev.off()
}

#run with all conditions
condition_chosen <- c(0,100,150,200,250,300,50,400,500,600,700)
conclude_growth_by_conditions()


#run up to 300 (first batch)
condition_chosen <- c(0,100,150,200,250,300,50)
conclude_growth_by_conditions(conditions_subset = condition_chosen)










