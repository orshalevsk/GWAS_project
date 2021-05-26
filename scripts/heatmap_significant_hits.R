#########plotting heatmap, using Haims script. a function so i can do it repeatedly


heatmap_sig_OGs <- function(OTU_of_interest ="all", sig_hits, phenotype_continious, output_file_path, phy_pattern_my_strains){
  ScriptPath="/Volumes/small_projects/hashkenazy/Pseudomonas/WithOr/Scripts/"
  source (paste (ScriptPath,"plot_phyletic_and_trait.R",sep=""))
  
  #store all possible trees path
  tree_phyl_all_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_my_strains_aligned_no_p4.C9.fasta.treefile"
  tree_phyl_non_OTU2_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_non_OTU2_strains_aligned.fasta.treefile"
  tree_phyl_OTU2_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_OTU2_strains_aligned.fasta.treefile"
  tree_phyl_OTU4_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_OTU4_strains_aligned.fasta.treefile"
  tree_phyl_OTU3_strains <- "/Volumes/small_projects/oshalev/GWAS_project/core_genome_OTU3_strains_aligned.fasta.treefile"
  
  
  #choose OTU to work with
  if (OTU_of_interest=="all"){
    tree=  tree_phyl_all_strains
    OTU <- paste("OTU",c(11,2,3,4,6,7,8,9,10,5,"_NEW"),sep = "")
  } else if (OTU_of_interest=="OTU2"){
    tree=  tree_phyl_OTU2_strains
    OTU <- "OTU2"
    
  } else if (OTU_of_interest=="OTU3"){
    tree=  tree_phyl_OTU3_strains
    OTU <- "OTU3"
    
  } else if (OTU_of_interest=="OTU4"){
    tree=  tree_phyl_OTU4_strains
    OTU <- "OTU4"
    
  }
  
  strains_by_OTU <- unique(dataset_GWAS$strain[dataset_GWAS$O_T_U %in% OTU])
  strains_by_OTU <- strains_by_OTU[!is.na(strains_by_OTU)]
  phy_pattern_my_strains_by_OTU <- phy_pattern_my_strains[phy_pattern_my_strains$Species %in% strains_by_OTU,]

  # build a phyletic pattern of only chosen OGs
  best_hits_phy_pattern <- phy_pattern_my_strains_by_OTU[,as.character(sig_hits)]
  row.names(best_hits_phy_pattern) <- phy_pattern_my_strains_by_OTU$Species
  
  
  #create a phenotype df with colors by the phenotype (required from this script, made by Haim)
  my_pheno <- as.data.frame(phy_pattern_my_strains_by_OTU[,phenotype_continious])
  colnames(my_pheno) <- "phenotype"
  row.names(my_pheno) <- phy_pattern_my_strains_by_OTU$Species
  my_pheno <- my_pheno[order(my_pheno$phenotype), , drop=F]
  #OPTIONAL: if green pixels of cdl50 is less than 0 => make it 0
  #my_pheno$phenotype[my_pheno$phenotype<0] <- 0
  
  my_pheno$color <- color.factor("green",(my_pheno$phenotype-min(my_pheno$phenotype))^3, max(my_pheno$phenotype-min(my_pheno$phenotype))^3)

  #plotting (this will create a file)
  plot_phyletic_stretch_and_trait(PhyP.df = as.matrix(best_hits_phy_pattern), TreeFile = tree, Traits_Colors_df = my_pheno,outFile = output_file_path)
  
  
}
