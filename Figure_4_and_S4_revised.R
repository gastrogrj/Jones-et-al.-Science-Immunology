# Jones et al. - Figure 4 and Supplementary Figure 4
#
# This script contains the blood/colon murine scRNA-seq analyses used in
# Fig. 4A-B and 4D-F, and Fig. S4B, S4E and S4F. Flow-cytometry panels
# (Fig. 4C/G and Fig. S4A/C/D/G) were generated in FlowJo/Prism and are not
# recreated here.
#
# Run from the repository root. All paths are relative and generated files are
# written beneath output/Figure_4_S4.

suppressPackageStartupMessages({
  library(DoubletFinder)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(igraph)
  library(monocle3)
  library(patchwork)
  library(readr)
  library(scales)
  library(Seurat)
  library(SeuratWrappers)
  library(SoupX)
  library(splines)
  library(tibble)
  library(tidyr)
  library(UCell)
})

source("mouse_analysis_helpers.R")
set.seed(1)
options(future.globals.maxSize = 30 * 1024^3)

# -----------------------------------------------------------------------------
# Inputs, analysis constants and outputs
# -----------------------------------------------------------------------------

COLON_NAIVE_10X_DIR <- file.path("data", "mouse_scRNAseq", "G1")
COLON_DSS_10X_DIR <- file.path("data", "mouse_scRNAseq", "G2")
BLOOD_NAIVE_10X_DIR <- file.path("data", "mouse_scRNAseq", "B1")
BLOOD_DSS_10X_DIR <- file.path("data", "mouse_scRNAseq", "B2")

OUTPUT_DIR <- file.path("output", "Figure_4_S4")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PARENT_CSF1R_CLUSTERS <- c("0", "2", "3", "5")
FIRST_EXCLUDED_CLUSTER <- "7"
SECOND_EXCLUDED_CLUSTER <- "10"
IAM_SEURAT_CLUSTER <- "8"

# Corrected for RNA-normalised input. This can alter singlet calls; use the
# original processed objects for exact figure reproduction. The historical
# merged-library DoubletFinder workflow is retained, with its limitation in README.
DOUBLET_FINDER_USES_SCT <- FALSE

# The historical Figure_4.R applied Harmony by tissue at both integration
# stages. This constant makes that consequential choice explicit.
HARMONY_GROUP <- "tissue"

save_panel <- function(filename,plot,width,height) {
  save_figure(file.path(OUTPUT_DIR,sub("\\.pdf$","",filename)),plot,width*25.4,height*25.4)
}

check_10x_directory <- function(path) {
  required <- c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")
  missing <- required[!dir.exists(file.path(path, required))]
  if (length(missing) > 0) {
    stop(
      "Missing 10x directories under ", path, ": ",
      paste(missing, collapse = ", ")
    )
  }
}

# -----------------------------------------------------------------------------
# SoupX and Seurat preprocessing used for Figure 4/S4
# -----------------------------------------------------------------------------

run_soupx <- function(sample_dir) {
  filtered_counts <- read_gene_expression(
    file.path(sample_dir, "filtered_feature_bc_matrix")
  )
  raw_counts <- read_gene_expression(
    file.path(sample_dir, "raw_feature_bc_matrix")
  )

  soup_channel <- SoupX::SoupChannel(raw_counts, filtered_counts)
  preliminary <- CreateSeuratObject(filtered_counts) |>
    SCTransform(verbose = FALSE, return.only.var.genes = FALSE) |>
    RunPCA(verbose = FALSE) |>
    RunUMAP(dims = 1:30, verbose = FALSE) |>
    FindNeighbors(dims = 1:30, verbose = FALSE) |>
    FindClusters(verbose = FALSE)

  soup_channel <- SoupX::setClusters(
    soup_channel,
    setNames(preliminary$seurat_clusters, colnames(preliminary))
  )
  soup_channel <- SoupX::setDR(
    soup_channel,
    Embeddings(preliminary, "umap")
  )
  soup_channel <- SoupX::autoEstCont(soup_channel)
  SoupX::adjustCounts(soup_channel, roundToInt = TRUE)
}

make_sample_object <- function(counts, sample_id, timepoint, tissue) {
  object <- CreateSeuratObject(counts, project = sample_id)
  object$Sample_ID <- sample_id
  object$timepoint <- timepoint
  object$tissue <- tissue
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^Rp[Sl]")
  subset(object, subset = nFeature_RNA > 200 & percent.mt < 15)
}

FINAL_OBJECT <- "data/objects/blood_colon_final.rds"
SIGNATURE_OBJECT <- "data/objects/blood_colon_before_cluster10.rds"
if(file.exists(FINAL_OBJECT)) {
  analysis_object <- join_rna(readRDS(FINAL_OBJECT),normalise=TRUE)
  combined_after_first_filter <- if(file.exists(SIGNATURE_OBJECT)) readRDS(SIGNATURE_OBJECT) else NULL
} else {
  warning("Raw rerun: QC >200 genes / <15% mitochondrial.")
input_directories <- c(
  COLON_NAIVE_10X_DIR,
  COLON_DSS_10X_DIR,
  BLOOD_NAIVE_10X_DIR,
  BLOOD_DSS_10X_DIR
)
invisible(lapply(input_directories, check_10x_directory))

colon_naive <- make_sample_object(
  run_soupx(COLON_NAIVE_10X_DIR),
  "Colon_Naive",
  "Naive",
  "Colon"
)
colon_dss <- make_sample_object(
  run_soupx(COLON_DSS_10X_DIR),
  "Colon_DSS",
  "DSS",
  "Colon"
)
blood_naive <- make_sample_object(
  run_soupx(BLOOD_NAIVE_10X_DIR),
  "Blood_Naive",
  "Naive",
  "Blood"
)
blood_dss <- make_sample_object(
  run_soupx(BLOOD_DSS_10X_DIR),
  "Blood_DSS",
  "DSS",
  "Blood"
)

combined_all <- merge(
  blood_naive,
  y = c(blood_dss, colon_naive, colon_dss),
  add.cell.ids = c("Blood_Naive", "Blood_DSS", "Colon_Naive", "Colon_DSS"),
  project = "COMBINED"
)

sample_objects <- SplitObject(combined_all, split.by = "Sample_ID")
sample_objects <- lapply(sample_objects, function(x) {
  x |>
    NormalizeData(verbose = FALSE) |>
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = 4000,
      verbose = FALSE
    )
})

integration_features <- SelectIntegrationFeatures(sample_objects)
# Record the command log required by DoubletFinder's RNA branch without
# replacing the original shared feature set used for PCA.
combined_all <- join_rna(combined_all, normalise = TRUE)
combined_all <- FindVariableFeatures(combined_all, selection.method = "vst",
                                    nfeatures = 4000, verbose = FALSE)
VariableFeatures(combined_all) <- integration_features

combined_all <- combined_all |>
  ScaleData(features = integration_features, verbose = FALSE) |>
  RunPCA(features = integration_features, verbose = FALSE) |>
  RunHarmony(group.by.vars = HARMONY_GROUP, verbose = FALSE) |>
  RunUMAP(reduction = "harmony", dims = 1:40, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:40, verbose = FALSE) |>
  FindClusters(resolution = 0.1, verbose = FALSE)

combined_all[["RNA"]] <- JoinLayers(combined_all[["RNA"]])

expected_doublets <- round(ncol(combined_all) * 0.07)
combined_all <- DoubletFinder::doubletFinder(
  combined_all,
  PCs = 1:10,
  pN = 0.25,
  pK = 0.09,
  nExp = expected_doublets,
  sct = DOUBLET_FINDER_USES_SCT
)

doublet_column <- grep(
  "^DF.classifications",
  colnames(combined_all[[]]),
  value = TRUE
)
if (length(doublet_column) != 1) {
  stop("Expected exactly one DoubletFinder classification column.")
}
combined_all <- subset(
  combined_all,
  cells = colnames(combined_all)[combined_all[[]][[doublet_column]] == "Singlet"]
)

# Direct in-memory subsetting replaces the historical CSV write/read round-trip.
Idents(combined_all) <- "seurat_clusters"
combined_csf1r <- subset(combined_all, idents = PARENT_CSF1R_CLUSTERS)

combined_csf1r <- combined_csf1r |>
  SCTransform(
    vars.to.regress = c("percent.mt", "percent.ribo"),
    return.only.var.genes = FALSE,
    verbose = FALSE
  ) |>
  RunPCA(verbose = FALSE) |>
  RunHarmony(
    assay.use = "SCT",
    group.by.vars = HARMONY_GROUP,
    verbose = FALSE
  ) |>
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindClusters(resolution = 0.6, verbose = FALSE)

# Remove cluster 7, repeat dimensional reduction/clustering, then remove
# cluster 10 exactly as in the supplied Figure_4.R.
Idents(combined_csf1r) <- "seurat_clusters"
keep_after_first_filter <- setdiff(levels(Idents(combined_csf1r)), FIRST_EXCLUDED_CLUSTER)
combined_after_first_filter <- subset(combined_csf1r, idents = keep_after_first_filter) |>
  RunPCA(verbose = FALSE) |>
  RunHarmony(
    assay.use = "SCT",
    group.by.vars = HARMONY_GROUP,
    verbose = FALSE
  ) |>
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindClusters(resolution = 0.8, verbose = FALSE)

Idents(combined_after_first_filter) <- "seurat_clusters"
keep_after_second_filter <- setdiff(
  levels(Idents(combined_after_first_filter)),
  SECOND_EXCLUDED_CLUSTER
)
analysis_object <- subset(
  combined_after_first_filter,
  idents = keep_after_second_filter
)

if (!IAM_SEURAT_CLUSTER %in% levels(Idents(analysis_object))) {
  stop(
    "Expected IAM Seurat cluster ", IAM_SEURAT_CLUSTER,
    " was not present after filtering."
  )
}

DefaultAssay(analysis_object) <- "RNA"
analysis_object[["RNA"]] <- JoinLayers(analysis_object[["RNA"]])

saveRDS(
  analysis_object,
  file.path(OUTPUT_DIR, "combined_blood_colon_Csf1r_processed.rds"),
  compress = "xz"
)


saveRDS(combined_after_first_filter,file.path(OUTPUT_DIR,"blood_colon_before_cluster10.rds"),compress="xz")

}
if(!all(c("tissue","timepoint","seurat_clusters")%in%colnames(analysis_object[[]])))
  stop("Blood/colon object requires tissue, timepoint and seurat_clusters.")
saveRDS(analysis_object,file.path(OUTPUT_DIR,"combined_blood_colon_Csf1r_processed.rds"),compress="xz")
export_cell_manifest(analysis_object,file.path(OUTPUT_DIR,"blood_colon_cell_manifest.csv"))

# -----------------------------------------------------------------------------
# Figure 4A: blood/colon UMAP, clusters and defining features
# -----------------------------------------------------------------------------

p_fig4a_tissue <- DimPlot(
  analysis_object,
  group.by = "tissue",
  raster = TRUE,
  pt.size = 0.35
) +
  theme_void()

p_fig4a_clusters <- DimPlot(
  analysis_object,
  group.by = "seurat_clusters",
  raster = TRUE,
  pt.size = 0.35
) +
  NoLegend() +
  theme_void()

fig4a_features <- c("Ly6c2", "Ccr2", "Cd163", "H2-Aa", "Cxcl9", "Itgax")
stopifnot(all(fig4a_features %in% rownames(analysis_object)))

p_fig4a_features <- feature_grid(analysis_object,fig4a_features,ncol=3)

save_panel(
  "Fig4A_blood_colon_UMAP.pdf",
  (p_fig4a_tissue + p_fig4a_clusters) / p_fig4a_features,
  width = 10,
  height = 7.2
)


# Figure 4B: fixed signatures from the ORIGINAL Figure_4.R (16 mouse genes,
# 14 mapped human genes), not signatures recalculated from the integrated data.
mouse_iam_signature <- c("Cxcl10","Cxcl9","Acod1","Slpi","Gbp2","Gbp5","Inhba","Sod2",
  "Il1rn","Saa3","Chil3","Clec4e","Plac8","Cfb","Cd274","Nampt")
human_iam_signature <- c("Acod1","Irf1","Stat1","Aqp9","Vcan","S100a8","Acsl1",
  "Trem1","Ccr1","Fpr1","S100a9","Cd300e","S100a10","S100a4")
check_features(analysis_object,c(mouse_iam_signature,human_iam_signature))
write_csv(tibble(gene=mouse_iam_signature),file.path(OUTPUT_DIR,"Fig4B_mouse_signature.csv"))
write_csv(tibble(gene=human_iam_signature),file.path(OUTPUT_DIR,"Fig4B_human_signature_mouse_orthologues.csv"))
DefaultAssay(analysis_object)<-"RNA"
analysis_object<-AddModuleScore(analysis_object,features=list(mouse_iam_signature),assay="RNA",
  name="Mouse_IAM_",ctrl=5,seed=1)
analysis_object<-AddModuleScore(analysis_object,features=list(human_iam_signature),assay="RNA",
  name="Human_IAM_",ctrl=5,seed=1)
p_mouse<-FeaturePlot(analysis_object,features="Mouse_IAM_1",min.cutoff=.05,max.cutoff=2,
  order=TRUE,raster=TRUE)+labs(title="Mouse Cxcl9/10 signature",colour="Module score")+theme_void(base_size=7)
p_human<-FeaturePlot(analysis_object,features="Human_IAM_1",min.cutoff=.2,max.cutoff=1,
  order=TRUE,raster=TRUE)+labs(title="Human CXCL9/10 signature",colour="Module score")+theme_void(base_size=7)
save_figure(file.path(OUTPUT_DIR,"Fig4B_mouse_human_IAM_module_scores"),p_mouse+p_human,85,50)
write_csv(tibble(cell=colnames(analysis_object),mouse=analysis_object$Mouse_IAM_1,
                 human=analysis_object$Human_IAM_1),file.path(OUTPUT_DIR,"Fig4B_module_scores.csv"))

# S4B used the object BEFORE cluster 10 removal in the original script.
if(!is.null(combined_after_first_filter)) {
  signature_object<-combined_after_first_filter
  Idents(signature_object)<-"seurat_clusters"
  DefaultAssay(signature_object)<-"SCT"
  markers<-FindMarkers(signature_object,ident.1="8",assay="SCT",
    logfc.threshold=.25,min.pct=.25,test.use="wilcox") |> rownames_to_column("gene")
  signature_sizes<-c(20,50,100,200)
  # Original sig_list uses adjusted P < .05 and ranks by fold change;
  # the separate >.5 signature variable was not used to build these four lists.
  selected<-markers |> filter(p_val_adj<.05,avg_log2FC>0) |> arrange(desc(avg_log2FC))
  if(nrow(selected)<200) stop("Fewer than 200 positive significant S4B markers.")
  signatures<-setNames(lapply(signature_sizes,function(n) head(selected$gene,n)),paste0("C8_sig",signature_sizes))
  write_csv(markers,file.path(OUTPUT_DIR,"FigS4B_cluster8_markers.csv"))
  write_csv(bind_rows(lapply(names(signatures),function(n) tibble(signature=n,gene=signatures[[n]]))),
    file.path(OUTPUT_DIR,"FigS4B_signature_genes.csv"))
  signature_object<-AddModuleScore(signature_object,features=signatures,assay="SCT",
    name="Audit_sig",nbin=24,ctrl=100,seed=1)
  add_cols<-paste0("C8_sig",signature_sizes)
  for(i in seq_along(add_cols)) signature_object[[add_cols[i]]]<-signature_object[[paste0("Audit_sig",i)]][,1]
  signature_object<-join_rna(signature_object)
  uu<-as.data.frame(UCell::ScoreSignatures_UCell(
    mat=GetAssayData(signature_object,assay="RNA",layer="data"),features=signatures))
  uc_cols<-paste0(names(signatures),"_UCell")
  if(!setequal(rownames(uu),colnames(signature_object))) stop("UCell cell barcodes do not match.")
  # Explicit barcode alignment prevents correlations between mismatched cells.
  uu<-uu[colnames(signature_object),uc_cols,drop=FALSE]
  signature_object<-AddMetaData(signature_object,metadata=uu)
  corr<-expand.grid(add_size=signature_sizes,ucell_size=signature_sizes)
  corr$r<-mapply(function(a,u) cor(signature_object[[]][[paste0("C8_sig",a)]],
    signature_object[[]][[paste0("C8_sig",u,"_UCell")]],method="pearson"),corr$add_size,corr$ucell_size)
  write_csv(corr,file.path(OUTPUT_DIR,"FigS4B_AddModuleScore_UCell_correlations.csv"))
  export_cell_manifest(signature_object,file.path(OUTPUT_DIR,"FigS4B_cell_scores.csv"))
  pp<-FeaturePlot(signature_object,features=add_cols,ncol=2,min.cutoff="q05",order=TRUE,raster=TRUE)&theme_void(base_size=7)
  pc<-ggplot(corr,aes(factor(add_size),factor(ucell_size),fill=r))+geom_tile(colour="white")+
    geom_text(aes(label=sprintf("%.2f",r)),size=7/ggplot2::.pt)+
    scale_fill_gradient(low="white",high="red",limits=c(0,1),name="Pearson r")+coord_fixed()+
    labs(x="AddModuleScore size",y="UCell size")+theme_minimal(base_size=7)
  save_figure(file.path(OUTPUT_DIR,"FigS4B_signature_sizes_and_concordance"),pp+pc,150,90)
} else message("S4B needs data/objects/blood_colon_before_cluster10.rds from the original analysis.")

# -----------------------------------------------------------------------------
# Reproducible Monocle3 trajectory helpers
# -----------------------------------------------------------------------------


ROOTS_FILE <- "data/annotations/trajectory_roots.csv"
TRAJECTORY_OBJECTS <- "data/objects/Monocle3_trajectory_objects.rds"
archived_trajectories<-if(file.exists(TRAJECTORY_OBJECTS)) readRDS(TRAJECTORY_OBJECTS) else list()
root_table<-if(file.exists(ROOTS_FILE)) read_csv(ROOTS_FILE,show_col_types=FALSE,
  col_types=cols(.default=col_character())) else data.frame()

run_trajectory <- function(object,graph_resolution,name) {
  if(name %in% names(archived_trajectories)) {
    cds<-archived_trajectories[[name]]
    if(!setequal(colnames(cds),colnames(object)))
      stop("Archived trajectory cell set differs from the selected Seurat object: ",name)
    if(!"gene_short_name"%in%colnames(SummarizedExperiment::rowData(cds)))
      SummarizedExperiment::rowData(cds)$gene_short_name<-rownames(cds)
    return(cds)
  }
  object<-join_rna(object)
  if(!all(c("trajectory","root_cell")%in%names(root_table)))
    stop("For exact trajectory reproduction supply archived Monocle objects or trajectory_roots.csv; see README.")
  roots<-root_table$root_cell[root_table$trajectory==name]
  if(!length(roots) || any(!roots%in%colnames(object)))
    stop("Provide valid original root-cell barcodes for ",name," in ",ROOTS_FILE)
  if(grepl("blood_colon",name) && any(object[[]][roots,"tissue"]!="Blood"))
    stop("Blood-colon roots must be annotated blood monocytes.")
  cds<-SeuratWrappers::as.cell_data_set(object)
  SummarizedExperiment::rowData(cds)$gene_short_name<-rownames(cds)
  cds<-monocle3::cluster_cells(cds,reduction_method="UMAP",resolution=graph_resolution,random_seed=1)
  cds<-monocle3::learn_graph(cds,use_partition=TRUE)
  cds<-monocle3::order_cells(cds,reduction_method="UMAP",root_cells=roots)
  cds
}
trajectory_plot <- function(cds, colour_by) {
  plot_cells(
    cds,
    color_cells_by = colour_by,
    label_cell_groups = FALSE,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    label_roots = FALSE,
    show_trajectory_graph = TRUE,
    trajectory_graph_color = "grey35",
    cell_size = 0.5
  )
}

# -----------------------------------------------------------------------------
# Figure 4D: trajectories from blood to colon, split by condition
# -----------------------------------------------------------------------------

trajectory_objects <- SplitObject(analysis_object, split.by = "timepoint")
cds_naive <- run_trajectory(trajectory_objects[["Naive"]], graph_resolution = 0.00078,name="naive_blood_colon")
cds_dss <- run_trajectory(trajectory_objects[["DSS"]], graph_resolution = 0.003,name="dss_blood_colon")

p_fig4d_naive_time <- trajectory_plot(cds_naive, "pseudotime") + ggtitle("Naive")
p_fig4d_naive_cluster <- trajectory_plot(cds_naive, "seurat_clusters") + ggtitle("Naive")
p_fig4d_dss_time <- trajectory_plot(cds_dss, "pseudotime") + ggtitle("Acute colitis")
p_fig4d_dss_cluster <- trajectory_plot(cds_dss, "seurat_clusters") + ggtitle("Acute colitis")

save_panel(
  "Fig4D_blood_colon_trajectories.pdf",
  (p_fig4d_naive_time + p_fig4d_naive_cluster) |
    (p_fig4d_dss_time + p_fig4d_dss_cluster),
  width = 12,
  height = 5.8
)

# -----------------------------------------------------------------------------
# Figure 4E-F: gene expression on the trajectory
# -----------------------------------------------------------------------------

fig4e_genes <- c("Cxcl9", "Acod1", "Il1b", "Nod2")
# The displayed Figure 4E is a Monocle gene overlay with a % of maximum
# scale, not a Seurat RNA feature plot. Use the DSS trajectory as in the original.
p_fig4e<-monocle3::plot_cells(cds_dss,genes=fig4e_genes,show_trajectory_graph=FALSE,
  label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,
  norm_method="log",scale_to_range=TRUE)
save_figure(file.path(OUTPUT_DIR,"Fig4E_trajectory_gene_features"),p_fig4e,145,45)

fig4f_genes <- c("Il1b", "Acod1")
missing_fig4f_genes <- setdiff(
  fig4f_genes,
  SummarizedExperiment::rowData(cds_dss)$gene_short_name
)
if (length(missing_fig4f_genes) > 0) {
  stop("Trajectory genes not found: ", paste(missing_fig4f_genes, collapse = ", "))
}
trend_cells_file <- "data/annotations/Fig4F_S4F_cells.csv"
if(!file.exists(trend_cells_file)) stop("Provide the original displayed-branch cell barcodes in ",trend_cells_file,
  "; the original Figure 4 script also contains a cluster-11-only fit, so all-cell trends cannot be assumed equivalent.")
trend_cells<-read_csv(trend_cells_file,show_col_types=FALSE)$cell
if(anyDuplicated(trend_cells) || !length(trend_cells) || any(!trend_cells%in%colnames(cds_dss)))
  stop("Invalid Figure 4F / S4F cell selection.")
finite_trend<-is.finite(monocle3::pseudotime(cds_dss)[trend_cells])
if(!all(finite_trend)) stop("Selected trend cells include unreachable pseudotime values.")
cds_fig4f <- cds_dss[SummarizedExperiment::rowData(cds_dss)$gene_short_name %in% fig4f_genes, trend_cells]
p_fig4f <- plot_genes_in_pseudotime(
  cds_fig4f,
  cell_size = 0.4,
  min_expr = 0,
  color_cells_by = "pseudotime",
  trend_formula = "~ splines::ns(pseudotime, df = 3)",
  label_by_short_name = TRUE
)
save_panel("Fig4F_Il1b_Acod1_over_pseudotime.pdf", p_fig4f, width = 5.5, height = 4.5)

# -----------------------------------------------------------------------------
# Figure S4E: colon-only trajectories, split by condition
# -----------------------------------------------------------------------------

colon_only <- subset(analysis_object, subset = tissue == "Colon")
colon_trajectory_objects <- SplitObject(colon_only, split.by = "timepoint")
cds_colon_naive <- run_trajectory(
  colon_trajectory_objects[["Naive"]],
  graph_resolution = 0.00078,name="naive_colon_only"
)
cds_colon_dss <- run_trajectory(
  colon_trajectory_objects[["DSS"]],
  graph_resolution = 0.003,name="dss_colon_only"
)

p_s4e_naive_time <- trajectory_plot(cds_colon_naive, "pseudotime") + ggtitle("Naive")
p_s4e_naive_cluster <- trajectory_plot(cds_colon_naive, "seurat_clusters") + ggtitle("Naive")
p_s4e_dss_time <- trajectory_plot(cds_colon_dss, "pseudotime") + ggtitle("DSS")
p_s4e_dss_cluster <- trajectory_plot(cds_colon_dss, "seurat_clusters") + ggtitle("DSS")

save_panel(
  "FigS4E_colon_only_trajectories.pdf",
  (p_s4e_naive_time + p_s4e_naive_cluster) |
    (p_s4e_dss_time + p_s4e_dss_cluster),
  width = 12,
  height = 5.8
)


for(nm in c("Naive","DSS")) {
  cc<-if(nm=="Naive") cds_colon_naive else cds_colon_dss
  keep<-is.finite(monocle3::pseudotime(cc))
  pp<-monocle3::plot_genes_in_pseudotime(cc[c("Cxcl9","Acod1"),keep],
    min_expr=0,cell_size=.3,color_cells_by="pseudotime",trend_formula="~ splines::ns(pseudotime, df = 3)")
  save_figure(file.path(OUTPUT_DIR,paste0("FigS4E_",nm,"_Cxcl9_Acod1_curves")),pp,75,55)
}

# -----------------------------------------------------------------------------
# Figure S4F: Slamf7 over pseudotime in the combined DSS blood/colon trajectory
# -----------------------------------------------------------------------------

if (!"Slamf7" %in% SummarizedExperiment::rowData(cds_dss)$gene_short_name) {
  stop("Slamf7 was not found in the DSS trajectory object.")
}
cds_slamf7 <- cds_dss[SummarizedExperiment::rowData(cds_dss)$gene_short_name == "Slamf7", trend_cells]
p_s4f <- plot_genes_in_pseudotime(
  cds_slamf7,
  cell_size = 0.5,
  min_expr = 0,
  color_cells_by = "pseudotime",
  trend_formula = "~ splines::ns(pseudotime, df = 3)",
  label_by_short_name = TRUE
)
save_panel("FigS4F_Slamf7_over_pseudotime.pdf", p_s4f, width = 5.2, height = 3.8)

saveRDS(
  list(
    naive_blood_colon = cds_naive,
    dss_blood_colon = cds_dss,
    naive_colon_only = cds_colon_naive,
    dss_colon_only = cds_colon_dss
  ),
  file.path(OUTPUT_DIR, "Monocle3_trajectory_objects.rds"),
  compress = "xz"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)

for(nm in c("naive_blood_colon","dss_blood_colon","naive_colon_only","dss_colon_only")) {
  cc<-get(c(naive_blood_colon="cds_naive",dss_blood_colon="cds_dss",
            naive_colon_only="cds_colon_naive",dss_colon_only="cds_colon_dss")[[nm]])
  write_csv(tibble(cell=colnames(cc),pseudotime=as.numeric(monocle3::pseudotime(cc)),
                   partition=as.character(monocle3::partitions(cc))),file.path(OUTPUT_DIR,paste0(nm,"_pseudotime.csv")))
}
