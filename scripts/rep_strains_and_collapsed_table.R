dataset_GWAS_delta_pixels$strain[dataset_GWAS_delta_pixels$OTU=="OTU5"]
my_strains <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Pseudomonas_Phyletic_Pattern_all_geneClusters.Clusterd_by_Jaccard_cd_hit_like.0.99_EDIT_MY_STRAINS.csv")
my_strains$cluster_representative

strains_df <- data.frame("rep_strain"=NULL, "collapsed_strains"=NULL,"cluster_size"=NULL, "OTU"=NULL)
for (strain in unique(dataset_GWAS_delta_pixels$strain)){
  if (strain %in% unique(my_strains$cluster_representative)){
    strains <- my_strains[my_strains$cluster_representative==strain,c("X",paste("X.",c(1:141),sep = ""))]
    strains <- as.character(unlist(strains))
    strains <- strains[strains!=""]
    strains <- paste(strains, collapse = "|")
    cluster_size <- my_strains$cluster_size[my_strains$cluster_representative==strain]
    OTU <- unique(dataset_GWAS_delta_pixels$OTU[dataset_GWAS_delta_pixels$strain==strain])
    temp_df <- data.frame("rep_strain"=strain, "collapsed_strains"=strains,"cluster_size"=cluster_size, "OTU"=OTU)
    strains_df <- rbind(strains_df,temp_df)
  }
  else {
    next
  }
}

strains_df_commensals <- strains_df[strains_df$OTU!="OTU5",]
#! write.csv2(x = strains_df_commensals, file = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/rep_and_collapsed_strains.csv")
