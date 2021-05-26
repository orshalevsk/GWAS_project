###This one is to return the GWAS database, as median per dpi, so that a scatterplot can be nicley made

scatter_plot_data <- function(GWAS_data){
  # now lets scatter plot this, using mean or median
  MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
  }
  
  #now lets create the dataframe that can show mean/median + error bars
  scatter_weight_df <- data.frame(t(rep(NA,5)))
  full_dataset <- GWAS_data
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
  scatter_weight_df$OTU <- full_dataset$O_T_U[match(scatter_weight_df$treatment,full_dataset$treatment)]
  scatter_weight_df$strain <- full_dataset$strain[match(scatter_weight_df$treatment,full_dataset$treatment)]
  return(scatter_weight_df)
}

