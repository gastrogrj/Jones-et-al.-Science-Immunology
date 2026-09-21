# The original notebook did not contain the displayed method-comparison panel.
# Compute it from the two actual top-200 marker exports, never fixed values.
suppressPackageStartupMessages({library(readr);library(dplyr);library(ggplot2);library(ggrepel)})
source("mouse_analysis_helpers.R")
out <- "output/Figure_S2E"
s <- read_csv(file.path(out,"input/Seurat_IAM_top200.csv"),show_col_types=FALSE)
p <- read_csv(file.path(out,"Scanpy_IAM_top200.csv"),show_col_types=FALSE)
stopifnot(all(c("gene","avg_log2FC")%in%names(s)),all(c("gene","avg_log2FC")%in%names(p)))
m <- inner_join(select(s,gene,Seurat_log2FC=avg_log2FC),select(p,gene,Scanpy_log2FC=avg_log2FC),by="gene")
if(nrow(m)<3) stop("Fewer than three shared markers; inspect the selected clusters.")
a <- cor.test(m$Seurat_log2FC,m$Scanpy_log2FC,method="pearson")
b <- cor.test(m$Seurat_log2FC,m$Scanpy_log2FC,method="spearman",exact=FALSE)
write_csv(m,file.path(out,"FigS2E_shared_markers.csv"))
write_csv(data.frame(shared_genes=nrow(m),Seurat_markers=nrow(s),Scanpy_markers=nrow(p),
  Pearson_r=unname(a$estimate),Pearson_P=a$p.value,Spearman_rho=unname(b$estimate),Spearman_P=b$p.value),
  file.path(out,"FigS2E_concordance_statistics.csv"))
g <- ggplot(m,aes(Seurat_log2FC,Scanpy_log2FC))+geom_point(size=.7)+geom_smooth(method="lm",se=FALSE,linewidth=.3)+
  geom_text_repel(data=m[order(m$Scanpy_log2FC,decreasing=TRUE)[seq_len(min(10,nrow(m)))],],
    aes(label=gene),size=7/ggplot2::.pt,seed=1)+
  labs(x="Seurat average log2FC",y="Scanpy log2FC",subtitle=sprintf("Shared genes = %d; Pearson r = %.2f; Spearman rho = %.2f",
    nrow(m),unname(a$estimate),unname(b$estimate)))+theme_classic(base_size=7)
save_figure(file.path(out,"FigS2E_marker_concordance"),g,85,70)
export_session(out)
