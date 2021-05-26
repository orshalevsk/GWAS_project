'''
The correlation between weight and green pixels
'''
library (ggplot2)


weight <- read.csv2("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Plant_GWAS/experiments/11_genotypes_trial_Nov_2018/results/weight/GWAS_weigt_20181207.csv")
pixels <- read.csv2("/Users/oshalev/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Plant_GWAS/experiments/11_genotypes_trial_Nov_2018/results/Green_pixels_area/processed_area_pixels_final.csv")

final_df <- merge(x = weight, y = pixels, by = "position")

final_df$weight <- as.numeric(as.character(final_df$weight))

cor.test(final_df$pixels_Talia_7dpi, final_df$weight)

final_df <- final_df[!final_df$weight > (mean(final_df$weight) + sd(final_df$weight)*2.5),]
final_df <- final_df[!final_df$pixels_Talia_7dpi > (mean(final_df$pixels_Talia_7dpi) + sd(final_df$pixels_Talia_7dpi)*2.5),]


#! pdf ("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/Cor_weight_pixels.pdf", useDingbats = F)
ggplot(data = final_df, aes(x = weight, y = pixels_Talia_7dpi))+
  geom_point()+
  geom_smooth(method = "lm") + 
  theme_classic()
#! dev.off()

