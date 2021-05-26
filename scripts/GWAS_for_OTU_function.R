###treeWAS function for my 103 strains GWAS
treeWAS_by_OTU <- function(tree_path, phyletic_pattern, phenotype,output_file, correction_pval=F){
  #library needed ape
  library(ape)
  library(treeWAS)
  library(adegenet)
  # first - read the tree
  tree_phyl_all_strains <- tree_path
  tree <- read.tree(file = tree_phyl_all_strains)
  tree <- keep.tip(tree,as.character(phyletic_pattern$Species))
  
  
  #limiting only to OGs with at least 2 occurences, and with at least one strain who miss it
  number_of_occurences_per_gene<- colSums(phyletic_pattern[,!names(phyletic_pattern) %in% c("cdl50","auc_l","dpi7_0median","Species")])
  at_least=names(number_of_occurences_per_gene[number_of_occurences_per_gene>1 & number_of_occurences_per_gene<length(phyletic_pattern$Species)])
  phyletic_pattern <- phyletic_pattern[,c("Species",at_least,"cdl50","auc_l","dpi7_0median")]
  
  
  #create a matrix of OGs only
  phy_pattern_my_strains_OGs <- phyletic_pattern[,grepl(x = names(phyletic_pattern) ,pattern = "GC")]
  row.names(phy_pattern_my_strains_OGs) <- phyletic_pattern$Species
  #making vector for the phenotype
  Traits_vector <- phyletic_pattern[,phenotype]
  names(Traits_vector) <- phyletic_pattern$Species
  if (correction_pval==F){
    treewas_results <- treeWAS(snps = as.matrix(phy_pattern_my_strains_OGs), phen = Traits_vector,tree=tree,
                                      test = c("terminal", "simultaneous", "subsequent"),filename.plot = output_file)
    
  }else {
    treewas_results <- treeWAS(snps = as.matrix(phy_pattern_my_strains_OGs), phen = Traits_vector,tree=tree,
                                      test = c("terminal", "simultaneous", "subsequent"),filename.plot = output_file,p.value.correct = F, p.value =correction_pval/length(phy_pattern_my_strains_OGs))
    
  }
  return(treewas_results)
}
