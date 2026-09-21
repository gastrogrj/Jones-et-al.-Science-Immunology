# Export the final Figure 2 cell set to Scanpy without converting assays or
# accidentally treating SCT residuals as RNA expression. Run after Figure 2.
suppressPackageStartupMessages({library(Seurat);library(readr);library(tibble);library(dplyr)})
source("mouse_analysis_helpers.R")
out <- "output/Figure_S2E/input"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
obj <- label_colon(join_rna(readRDS("output/Figure_2_S2/combined_colon_Csf1r_processed.rds"),normalise=TRUE))
Matrix::writeMM(GetAssayData(obj,assay="RNA",layer="data"),file.path(out,"logcounts.mtx"))
Matrix::writeMM(GetAssayData(obj,assay="RNA",layer="counts"),file.path(out,"counts.mtx"))
write_csv(tibble(gene=rownames(obj[["RNA"]])),file.path(out,"genes.csv"))
export_cell_manifest(obj,file.path(out,"cells.csv"))
# Source code's Seurat comparator: raw cluster 3 vs rest, all conditions.
a<-de_table(obj,cell_groups(obj,"4"),cell_groups(obj,c("1","2","3","5")),min_pct=.25,only_positive=TRUE)
a<-a |> filter(p_val_adj<.05,log2FC>0) |> arrange(desc(log2FC),gene)
write_csv(head(a,200),file.path(out,"Seurat_IAM_top200.csv"))
export_session(dirname(out))
