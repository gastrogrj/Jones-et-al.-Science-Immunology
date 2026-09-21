# Figure 3A and the mouse scRNA-seq panels S3A/B.
# 3A is reconstructed from the current legend: DSS cluster 4 versus cluster 2
# (the other Ccr2+/Ly6c2+ monocyte cluster). Confirm the comparator against the
# original DE export; no original Figure 3A DE code was present in the repository.
# S3A requires the original BLOOD-ONLY embedding, not the integrated Fig 4 UMAP.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(readr); library(patchwork)
})
source("mouse_analysis_helpers.R")
set.seed(1)
OUTPUT_DIR <- "output/Figure_3_mouse"
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
COLON_OBJECT <- "output/Figure_2_S2/combined_colon_Csf1r_processed.rds"
if(!file.exists(COLON_OBJECT)) COLON_OBJECT <- "data/objects/colon_final.rds"
if(!file.exists(COLON_OBJECT)) stop("Run Figure_2_and_S2_revised.R or supply colon_final.rds.")
obj <- label_colon(join_rna(readRDS(COLON_OBJECT),normalise=TRUE))

FIG3A_COMPARATOR <- "2"
iam <- cell_groups(obj,"4","DSS")
other <- cell_groups(obj,FIG3A_COMPARATOR,"DSS")
de <- de_table(obj,iam,other)
de$minus_log10_adjusted_P <- -log10(pmax(de$p_val_adj,1e-300))
de$significant <- de$p_val_adj<.05
write_csv(de,file.path(OUTPUT_DIR,"Fig3A_cluster4_vs_cluster2_DSS_DE.csv"))
write_csv(tibble(cell=c(iam,other),group=c(rep("IAM_cluster4",length(iam)),
  rep("Comparator_cluster2",length(other)))),file.path(OUTPUT_DIR,"Fig3A_cell_groups.csv"))
label_genes <- c("Cxcl9","Cxcl10","Acod1","Socs1","Il1a","Tnip1","H2-Eb1",
  "Pstpip2","Slpi","Gbp6","Ifit2","Ifit1","Mrc1","C1qc","C1qa","Tnf",
  "Arg2","Chil3","Il23a","Il12a")
p <- ggplot(de,aes(log2FC,minus_log10_adjusted_P))+geom_point(aes(colour=significant),size=.5)+
  geom_text_repel(data=de[de$gene%in%label_genes,],aes(label=gene),size=7/ggplot2::.pt,
    fontface="italic",max.overlaps=Inf,seed=1)+
  scale_colour_manual(values=c("FALSE"="grey70","TRUE"="#38A7D3"),guide="none")+
  labs(x="Average log2 fold change (cluster 4 / cluster 2)",y="-log10 adjusted P")+
  theme_classic(base_size=7)
save_figure(file.path(OUTPUT_DIR,"Fig3A_IAM_vs_other_monocytes"),p,85,70)

# S3B: current panel displays the colonic IAM cluster split by condition.
# The source selection was absent; make this explicit and export the barcodes.
s3b_cells <- cell_groups(obj,"4")
s3b <- subset(obj,cells=s3b_cells)
s3b$timepoint <- factor(s3b$timepoint,levels=c("Naive","DSS"))
p <- VlnPlot(s3b,features="Acod1",assay="RNA",group.by="timepoint",pt.size=.1)+
  labs(x=NULL,y="Acod1 RNA log1p(counts per 10,000)")+NoLegend()+theme_classic(base_size=7)
save_figure(file.path(OUTPUT_DIR,"FigS3B_Acod1_colon_IAM"),p,55,55)
write_csv(tibble(cell=s3b_cells),file.path(OUTPUT_DIR,"FigS3B_cell_selection.csv"))

BLOOD_ONLY_OBJECT <- "data/objects/blood_only_S3A.rds"
if(file.exists(BLOOD_ONLY_OBJECT)) {
  blood <- join_rna(readRDS(BLOOD_ONLY_OBJECT),normalise=TRUE)
  if(!"umap"%in%names(blood@reductions)) stop("S3A blood-only UMAP is missing.")
  if("tissue"%in%names(blood[[]]) && any(blood$tissue!="Blood"))
    stop("S3A object contains non-blood cells.")
  p <- feature_grid(blood,c("Ly6c2","Acod1"),ncol=2)
  save_figure(file.path(OUTPUT_DIR,"FigS3A_blood_Ly6c2_Acod1"),p,95,50)
  export_cell_manifest(blood,file.path(OUTPUT_DIR,"FigS3A_blood_cells.csv"))
} else message("S3A needs the original blood-only Seurat object at ",BLOOD_ONLY_OBJECT,
               ". Its upstream workflow was not supplied; the Figure 4 embedding is not substituted.")
export_session(OUTPUT_DIR)
