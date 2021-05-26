#This one is just to load the data from the GWAS , and re-order it



#loading April exp, and subsetting to treatments Control, 1, 2 and 3.
database_return <- function(){
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
  
  dataset_GWAS <- dataset_GWAS[dataset_GWAS$dpi!="7_nolid",]
  
  dataset_GWAS$strain[dataset_GWAS$treatment=="control"] <- "control"
  dataset_GWAS$strain[dataset_GWAS$treatment=="pathogen_only"] <- "p4.C9_pathogen"
  dataset_GWAS$strain[dataset_GWAS$treatment=="boiled_strain_17+p"] <- "p5.F2_boiled"
  
  dataset_GWAS <- merge(x= dataset_GWAS, y= metadata[,c("year","location","lant","leaf","O_T_U","strain")], by = "strain", all.x=TRUE)
  dataset_GWAS$strain <- as.factor(dataset_GWAS$strain)
  
  return(dataset_GWAS)
}

