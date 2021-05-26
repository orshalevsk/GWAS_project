'''
Venn diagram for the positive correlations
'''
library(VennDiagram)
library(RColorBrewer)

positive_sig_hits <- read.csv2("~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/TableS2/other_related_tables/all_treeWAS_results_2nd_replicate_annotated_edited_proximity_BY_OTU_POSITIVE_COR.csv")

all <- as.character(positive_sig_hits$Var1[positive_sig_hits$OTU=="all"])
all <- all[all!="GC00002948_r1_1"]
all <- unique(all)
OTU3 <-  as.character(positive_sig_hits$Var1[positive_sig_hits$OTU=="OTU3"])
OTU4 <- as.character(positive_sig_hits$Var1[positive_sig_hits$OTU=="OTU4"])


myCol <- brewer.pal(3, "Pastel2")


venn.diagram(x = list("all"=all, "OTU3"=OTU3, "OTU4"=OTU4), scaled = TRUE, 
             filename = "~/ownCloud/My papers/bacterial_GWAS_paper/figures/Final_figures/FigureS3/temp.pdf",
             fill = myCol,
             cat.default.pos = "outer",
             sub.cex = 20,
             main.cex = 20,
             cat.pos = c(-27, 27, 135))
