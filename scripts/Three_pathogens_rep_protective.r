library(ggplot2)
dataframe <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS_few_pathogens_by_protective_rep2_OCT_2019/plan/full_dataframe_GWAS_protective_oct_2019.csv")
results <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/GWAS_project/Bacterial_rescue_GWAS_few_pathogens_by_protective_rep2_OCT_2019/results/results_dataframe.csv")
full_pathogens_test <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/first_protective_exp.csv")

results$dpi7dpi0 <- results$Green_pixels[results$dpi==7]-results$Green_pixels[results$dpi==1]
results_median <- results[results$dpi==7,]
final_dataframe <- merge(x=dataframe,y=results_median,by="position")


final_dataframe <- final_dataframe[,colnames(final_dataframe) %in% c("strain","pathogen","dpi7dpi0", "O_T_U")]
colnames(final_dataframe) <- c("pathogen","test_strain","OTU","median_delta_pixels")
final_dataframe$exp <- "1"
final_dataframe$median_delta_pixels <- final_dataframe$median_delta_pixels*0.12

colnames(full_pathogens_test) <- c("X", "test_strain", "median_delta_pixels","OTU","first_exp","pathogen")
full_pathogens_test <- full_pathogens_test[,colnames(full_pathogens_test) %in% c("test_strain", "median_delta_pixels","OTU","pathogen")]
full_pathogens_test$exp <- "2"

full_pathogens_test <- full_pathogens_test[,match(colnames(full_pathogens_test),colnames(final_dataframe))]
full_pathogens_test$OTU <- as.character(full_pathogens_test$OTU)



final_dataframe$median_delta_pixels[final_dataframe$pathogen=="p25.A12" & final_dataframe$test_strain=="control"]

full_pathogens_test$test_strain <- as.character(full_pathogens_test$test_strain)
full_pathogens_test$test_strain[full_pathogens_test$pathogen=="p25.A12" & full_pathogens_test$test_strain=="p25.A12"] <- "control"
full_pathogens_test$test_strain[full_pathogens_test$pathogen=="p25.C2" & full_pathogens_test$test_strain=="p25.C2"] <- "control"
full_pathogens_test$test_strain[full_pathogens_test$pathogen=="p4.C9" & full_pathogens_test$test_strain=="p4.C9"] <- "control"
full_pathogens_test$test_strain[full_pathogens_test$pathogen=="p4.C9_first" & full_pathogens_test$test_strain=="p4.C9"] <- "control"



both_exp <- rbind(full_pathogens_test,final_dataframe)
levels(both_exp$pathogen) <- c("p25.A12","p25.C2","p4.C9","p4.C9","control")





View(both_exp)



ggplot(aes(x=pathogen,y=median_delta_pixels), data = both_exp)+
  geom_boxplot(aes(fill=test_strain)) + 
  facet_grid(OTU ~ .)

results$Green_pixels[results$dpi==1]

