# Mouse panels 5B/C/D/H. Run after Figure_2_and_S2_revised.R.
# 5C/D preserve the aggregated-expression implementation in the original Figure 2
# These are descriptive analyses of pooled libraries, not replicate-level tests.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(fgsea); library(decoupleR); library(progeny); library(patchwork)
})
source("mouse_analysis_helpers.R")
check_packages(c("fmsb","xml2"))
set.seed(1)
OUTPUT_DIR <- "output/Figure_5_mouse"
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
OBJECT <- "output/Figure_2_S2/combined_colon_Csf1r_processed.rds"
if(!file.exists(OBJECT)) OBJECT <- "data/objects/colon_final.rds"
if(!file.exists(OBJECT)) stop("Run Figure 2 or supply data/objects/colon_final.rds.")
obj <- label_colon(join_rna(readRDS(OBJECT),normalise=TRUE))

# 5B: IAM versus all other Csf1r+ clusters in DSS, as specified by the legend.
# Do not confuse this with the monocyte-only comparator in Figure 3A.
iam <- cell_groups(obj,"4","DSS")
others <- cell_groups(obj,c("1","2","3","5"),"DSS")
de <- de_table(obj,iam,others,min_pct=.25)
write_csv(de,file.path(OUTPUT_DIR,"Fig5B_IAM_vs_other_DSS_clusters_DE.csv"))
write_csv(tibble(cell=c(iam,others),group=c(rep("IAM",length(iam)),rep("Other_Csf1r",length(others)))),
          file.path(OUTPUT_DIR,"Fig5B_cell_groups.csv"))
# The original enrichment recipe and GMT version were not supplied. These
# declared values are reconstruction settings, not authenticated historical ones.
REACTOME_GMT <- "data/gmt/m2.cp.reactome.v2025.1.Mm.symbols.gmt"
if(file.exists(REACTOME_GMT)) {
  pathways <- fgsea::gmtPathways(REACTOME_GMT)
  ranked <- de |> filter(is.finite(log2FC)) |> arrange(desc(log2FC),gene)
  ranks <- setNames(ranked$log2FC,ranked$gene)
  write_csv(tibble(gene=names(ranks),log2FC=as.numeric(ranks)),file.path(OUTPUT_DIR,"Fig5B_ranked_genes.csv"))
  enrichment <- fgsea::fgseaMultilevel(pathways=pathways,stats=ranks,minSize=10,maxSize=2000,
                                      scoreType="std",eps=0,nPermSimple=10000)
  saveRDS(enrichment,file.path(OUTPUT_DIR,"Fig5B_Reactome_fgsea.rds"))
  ee <- as.data.frame(enrichment)
  ee$leadingEdge <- vapply(ee$leadingEdge,paste,collapse=";",FUN.VALUE=character(1))
  write_csv(ee,file.path(OUTPUT_DIR,"Fig5B_Reactome_all.csv"))
  display_names <- c("REACTOME_INTERFERON_ALPHA_BETA_SIGNALING","REACTOME_INTERLEUKIN_1_SIGNALING",
    "REACTOME_INTERFERON_GAMMA_SIGNALING","REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION","REACTOME_TOLL_LIKE_RECEPTOR_CASCADES")
  displayed <- ee |> filter(pathway%in%display_names,padj<.05) |> arrange(NES)
  missing <- setdiff(display_names,displayed$pathway)
  if(length(missing)) warning("Displayed Fig5B pathways absent or not FDR<.05 in this rerun: ",paste(missing,collapse=", "))
  if(nrow(displayed)) {
    displayed$label <- gsub("_"," ",sub("^REACTOME_","",displayed$pathway))
    displayed$label <- factor(displayed$label,levels=displayed$label)
    p<-ggplot(displayed,aes(NES,label))+geom_col(fill="#E9797E")+
      labs(x="Normalised enrichment score",y=NULL)+theme_classic(base_size=7)
    save_figure(file.path(OUTPUT_DIR,"Fig5B_Reactome_selected_pathways"),p,100,65)
  }
  interferon <- "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"
  if(interferon%in%names(pathways)) {
    p<-fgsea::plotEnrichment(pathways[[interferon]],ranks)+labs(title="IFN alpha/beta signalling")+
      theme_classic(base_size=7)
    save_figure(file.path(OUTPUT_DIR,"Fig5B_IFN_enrichment_curve"),p,80,55)
  }
  # Exact gene labels visible in the submitted leading-edge panel; verify that
  # they remain leading-edge members in the computed enrichment result.
  genes <- c("Socs1","Ifit1","Rsad2","Ifit3","Oasl1","Ifit2","Isg15","Usp18","Irf1","Irf7")
  le <- enrichment$leadingEdge[match(interferon,enrichment$pathway)]
  if(length(le) && !is.na(match(interferon,enrichment$pathway))) {
    h<-de |> filter(gene%in%genes,gene%in%le[[1]])
    if(nrow(h)!=length(genes)) warning("Some displayed Figure5B genes are absent from this rerun's leading edge.")
    h$gene<-factor(h$gene,levels=rev(genes))
    p<-ggplot(h,aes("IAM vs others",gene,fill=log2FC))+geom_tile()+
      scale_fill_gradientn(colours=c("#440154","#21918c","#fde725"),name="Average log2FC")+
      labs(x=NULL,y=NULL)+theme_minimal(base_size=7)+theme(panel.grid=element_blank(),axis.text.y=element_text(face="italic"))
    save_figure(file.path(OUTPUT_DIR,"Fig5B_IFN_leading_edge"),p,45,70)
    write_csv(h,file.path(OUTPUT_DIR,"Fig5B_displayed_leading_edge.csv"))
  }
  write_csv(tibble(file=REACTOME_GMT,md5=unname(tools::md5sum(REACTOME_GMT))),file.path(OUTPUT_DIR,"GMT_checksum.csv"))
} else message("Fig5B enrichment requires ",REACTOME_GMT,"; 5C/D/H continue.")

# Explicit averaging in linear normalised space, matching AverageExpression's
# behaviour for layer='data'. This is a cluster/library mean, not count-sum
# pseudobulk and not a mouse biological replicate.
obj$activity_group <- paste0("C",obj$figure_cluster,"-",obj$timepoint)
means <- as.matrix(AverageExpression(obj,assays="RNA",group.by="activity_group",layer="data",verbose=FALSE)$RNA)
write_csv(as.data.frame(means) |> rownames_to_column("gene"),file.path(OUTPUT_DIR,"cluster_condition_mean_expression.csv"))
dss_columns <- paste0("C",as.character(1:5),"-DSS")
stopifnot(all(dss_columns%in%colnames(means)))

# 5C: save the exact network used, so a future OmniPath update is not invisible.
NETWORK_FILE <- "data/networks/collectri_mouse.csv"
if(file.exists(NETWORK_FILE)) {
  network <- read_csv(NETWORK_FILE,show_col_types=FALSE)
} else {
  network <- decoupleR::get_collectri(organism="mouse",split_complexes=FALSE)
  message("Downloaded current CollecTRI network. Archive the exported CSV for the final release.")
}
stopifnot(all(c("source","target","mor")%in%names(network)))
write_csv(network,file.path(OUTPUT_DIR,"collectri_mouse_used.csv"))
ulm <- decoupleR::run_ulm(mat=means,network=network,.source="source",.target="target",.mor="mor",minsize=5)
write_csv(ulm,file.path(OUTPUT_DIR,"Fig5C_ULM_all_scores.csv"))
tf_order <- c("Stat3","Irf5","Irf4","Irf2","Irf9","Irf8")
cluster_order <- c("4","1","2","3","5")
dd<-ulm |> filter(condition%in%dss_columns,source%in%tf_order) |>
  mutate(cluster=sub("-DSS$","",sub("^C","",condition)),
         cluster=factor(cluster,levels=cluster_order),tf=factor(source,levels=rev(tf_order)))
if(nrow(dd)!=30) stop("Expected all six Figure 5C TFs across five DSS clusters.")
p<-ggplot(dd,aes(cluster,tf,fill=score))+geom_tile()+
  scale_fill_gradient2(low="#3b4cc0",mid="white",high="#b40426",midpoint=0,name="TF activity\nULM t statistic")+
  labs(x="Figure 2 cluster",y=NULL)+theme_minimal(base_size=7)+theme(panel.grid=element_blank())
save_figure(file.path(OUTPUT_DIR,"Fig5C_mouse_TF_activity"),p,60,70)

# 5D: scaling is fitted across ALL naive and DSS cluster profiles, then DSS is
# selected, preserving the original code. Do not scale DSS profiles alone.
scores<-progeny::progeny(expr=means,organism="Mouse",top=500,scale=TRUE,perm=1)
saveRDS(scores,file.path(OUTPUT_DIR,"Fig5D_PROGENy_all_scores.rds"))
write_csv(as.data.frame(scores) |> rownames_to_column("profile"),file.path(OUTPUT_DIR,"Fig5D_PROGENy_all_scores.csv"))
ss<-scores[dss_columns,,drop=FALSE]
rownames(ss)<-as.character(1:5)
pathway_order<-c("Hypoxia","Estrogen","NFkB","TNFa","JAK-STAT","EGFR","MAPK",
                 "TGFb","PI3K","VEGF","p53","Trail","Androgen","WNT")
export_progeny_radar_7pt(ss,file.path(OUTPUT_DIR,"Fig5D_mouse_7pt"),
  colours=FIGURE_CLUSTER_COLOURS,pathway_order=pathway_order,width_mm=58,height_mm=62)

# 5H: original source selects IAM cluster across the colon object, with no DSS
# restriction. Retain this selection; use Spearman as in the current legend.
hh<-FetchData(obj,vars=c("Il1b","Acod1"),cells=cell_groups(obj,"4"),layer="data")
result<-safe_spearman(hh$Il1b,hh$Acod1)
write_csv(result,file.path(OUTPUT_DIR,"Fig5H_Spearman.csv"))
write_csv(rownames_to_column(hh,"cell"),file.path(OUTPUT_DIR,"Fig5H_cell_expression.csv"))
p<-ggplot(hh,aes(Il1b,Acod1))+geom_point(colour="#38A7D3",size=.3,alpha=.6)+
  geom_density_2d(colour="grey30",linewidth=.2)+
  annotate("text",x=-Inf,y=Inf,hjust=-.05,vjust=1.1,size=7/ggplot2::.pt,
    label=paste0("Spearman rho = ",sprintf("%.2f",result$rho),"\nP ",format.pval(result$p_value,digits=3)))+
  labs(x="Il1b RNA log1p(counts per 10,000)",y="Acod1 RNA log1p(counts per 10,000)")+
  theme_classic(base_size=7)
save_figure(file.path(OUTPUT_DIR,"Fig5H_Il1b_Acod1_Spearman"),p,70,70)
writeLines(c("5B: reconstructed from legend; DSS IAM vs other DSS Csf1r+ clusters, RNA Wilcoxon, signed log2FC ranks, fgseaMultilevel, minSize=10, maxSize=2000, BH FDR<0.05.",
 "5C: CollecTRI ULM on cluster-condition linear-normalised means; six prespecified TFs, DSS shown.",
 "5D: Mouse PROGENy, top=500, scale=TRUE, perm=1; scale all condition-cluster profiles, plot DSS only; radar rings are within-pathway range percentages.",
 "5H: all colon IAM cluster cells; RNA data layer; two-sided Spearman, no multiple-testing correction for this single pair."),
 file.path(OUTPUT_DIR,"analysis_settings.txt"))
export_session(OUTPUT_DIR)
