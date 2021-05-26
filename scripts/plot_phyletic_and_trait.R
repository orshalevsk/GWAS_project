library(gplots)
library(dendextend)

# function for ploting stretches heatmap and traits as SideBar optionally orderd by phylogenetic tree (edited from Haim)
plot_phyletic_stretch_and_trait <- function(PhyPFile=NULL,Traits_Colors_df,outFile=NULL,categories=NULL,pallete=NULL,TreeFile=NULL,sep=",",subset_genes=NULL,PhyP.df=NULL, dend_height=10, dend_width=12, strains_to_label="") {
  
  # PhyP.df =  rep_strains_phyletic_pattern
  # Traits_Colors_df = OTU_colors_df
  # outFile = PhyleticPatten_with_OTU_order_byTree
  # categories = myOTU_rep
  # pallete = mypalette_rep
  # TreeFile =RepStrains_core_Tree
  
  
  # PhyPFile="/Volumes/small_projects/hashkenazy/Pseudomonas/WithOr/Jacc_le_0.9_ANI_grteq_99_acces1/Jacc_between_rep_strains_0.99_Jacc_less_0.9_ANI_grt_99_access1.OG_MSA_and_phylletic_modules.ONLY_REP_STRAINS.FinalStretches.grt1_and_missing_Phyletic.percent.csv"
  # sep=","
  # Traits_Colors_df=OTU_cat_df
  # PhyPFile=full_in
  
  # PhyP.df =  rep_strains_phyletic_pattern
  # Traits_Colors_df = OTU_colors_df
  # outFile = PhyleticPatten_with_OTU_order_byTree
  # categories = myOTU_rep
  # pallete = mypalette_rep
  # TreeFile =RepStrains_core_Tree
  Phyletic_Pattern=data.frame()
  if (is.null(PhyPFile)&is.null(PhyP.df))
  {
    print ("ERROR! at a dataframe or a file with the phyletic pattern is mandatory for ploting...")      
  }
  if (!is.null(PhyP.df))
  {
    Phyletic_Pattern=PhyP.df
  }
  else
  {
    Phyletic_Pattern = read.table(PhyPFile,sep=sep,as.is = T,header = T) # comma delimited
    row.names(Phyletic_Pattern)=Phyletic_Pattern[,1] # assume the species names are the first column
    Phyletic_Pattern=Phyletic_Pattern[,-1]
  }
  
  if (!is.null(subset_genes))
  {
    Phyletic_Pattern=Phyletic_Pattern[,subset_genes]
  }
  
  names(Traits_Colors_df)=c("trait","color")
  
  # preparing the general things for the heatmap ploting
  if (!is.null(outFile)) {
    pdf(outFile,height = dend_height,width = dend_width)
  }
  
  
  Breaks=seq(0,1,by=0.1)
  Colors=c("#ffffff",brewer.pal(9,"PuRd"))
  
  if(!is.null(TreeFile)){
    message("TreeFile was provided, so order the data according the tree :-)")
    library(ape)
    Tree=read.tree(TreeFile)
    # resolve multichotomic (multi2di) and make ultrametric
    Tree_brlen <- compute.brlen(multi2di(Tree), method="Grafen") #so that branches have all same length
    #turn the phylo tree to a dendrogram object
    hc <- as.hclust(Tree_brlen) #Compulsory step as as.dendrogram doesn't have a  method for phylo objects.
    dend <- as.dendrogram(hc)
    
    
    #OTU_strain$OTU2 <- as.factor(OTU_strain$OTU=="OTU2")
    #color_easy = c("black","red")[OTU_strain$OTU2]
    
    #cols_branches <- color_easy
    #dend <- color_branches(dend, k = length(color_easy)-1, col = cols_branches)
    #dend <- color_branches(dend, k = 4, col = cols_branches)
    #col_labels <- get_leaves_branches_col(dend)
    #col_labels <- col_labels[order(order.dendrogram(dend))]
    
    #initiate cols with all black
    #cols <- rep('black', length(hc$labels))
    #turn red the specified rows in tf
    #cols[hc$labels %in% strains_to_label] <- 'red'
 
    
    # order the data by the tree
    Trait_Ordered_By_Tree=merge(x=data.frame(species=Tree$tip),y=Traits_Colors_df,by.x="species",by.y="row.names",sort = F)
    RowSideCols=as.character(Trait_Ordered_By_Tree$color)
    Phyletic_Pattern_Ordered_by_Tree=merge(x=data.frame(species=Tree$tip),y=Phyletic_Pattern,by.x="species",by.y="row.names",sort = F)
    row.names(Phyletic_Pattern_Ordered_by_Tree)=Phyletic_Pattern_Ordered_by_Tree$species
    Phyletic_Pattern_Ordered_by_Tree=Phyletic_Pattern_Ordered_by_Tree[,-1]
    p=heatmap.2(x = as.matrix(Phyletic_Pattern_Ordered_by_Tree),
                Rowv=dend,dendrogram='row',
                col=Colors,breaks = Breaks,
                trace = 'none', density.info="none",
                RowSideColors=RowSideCols,
                margin=c(25,5),
                srtCol=45) 
  }
  else {
    # order traits by PhyleticPattern
    Trait_Ordered_By_PhyleticP=merge(x=data.frame(species=row.names(Phyletic_Pattern)),y=Traits_Colors_df,by.x="species",by.y="row.names",sort = F)
    RowSideCols=as.character(Trait_Ordered_By_PhyleticP$color)
    p=heatmap.2(x = as.matrix(Phyletic_Pattern),
                col=Colors,breaks = Breaks,
                trace = 'none', density.info="none",
                RowSideColors=RowSideCols,
                margin=c(25,5),
                srtCol=45)  
  }
  if (!is.null(categories)&!is.null(pallete))
  {
    p+legend("topright",      
             legend = categories,
             col = pallete, 
             lty= 1,             
             lwd = 5,           
             cex=.7)
  }
  if (!is.null(outFile))
  {
    dev.off()
  }
  
  # return(Phyletic_Pattern_Ordered_by_Tree,Trait_Ordered_By_Tree)
}  
#dev.off()