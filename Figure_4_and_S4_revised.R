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

# Preserves the setting in the supplied analysis. Because the merged object is
# RNA-normalised before this call and SCTransform is applied afterwards, this
# should be confirmed against the environment used for the final cell calls.
DOUBLET_FINDER_USES_SCT <- TRUE

# The historical Figure_4.R applied Harmony by tissue at both integration
# stages. This constant makes that consequential choice explicit.
HARMONY_GROUP <- "tissue"

save_panel <- function(filename, plot, width, height) {
  ggsave(
    file.path(OUTPUT_DIR, filename),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
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
  filtered_counts <- Seurat::Read10X(
    file.path(sample_dir, "filtered_feature_bc_matrix")
  )
  raw_counts <- Seurat::Read10X(
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
VariableFeatures(combined_all) <- integration_features

combined_all <- combined_all |>
  NormalizeData(verbose = FALSE) |>
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

p_fig4a_features <- FeaturePlot(
  analysis_object,
  features = fig4a_features,
  ncol = 3,
  order = TRUE,
  raster = TRUE
) &
  theme_void()

save_panel(
  "Fig4A_blood_colon_UMAP.pdf",
  (p_fig4a_tissue + p_fig4a_clusters) / p_fig4a_features,
  width = 10,
  height = 7.2
)

# -----------------------------------------------------------------------------
# IAM marker table and signatures used in Figure 4B/S4B
# -----------------------------------------------------------------------------

Idents(analysis_object) <- "seurat_clusters"
iam_markers <- FindMarkers(
  analysis_object,
  ident.1 = IAM_SEURAT_CLUSTER,
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.25,
  only.pos = TRUE
) |>
  rownames_to_column("gene")

marker_logfc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(iam_markers)
)[1]
if (is.na(marker_logfc_column)) stop("No IAM marker log-fold-change column found.")

iam_markers <- iam_markers |>
  mutate(log2FC = .data[[marker_logfc_column]]) |>
  filter(p_val_adj < 0.05, log2FC > 0.5) |>
  arrange(desc(log2FC))

write_csv(iam_markers, file.path(OUTPUT_DIR, "IAM_cluster8_markers.csv"))

signature_sizes <- c(20, 50, 100, 200)
if (nrow(iam_markers) < max(signature_sizes)) {
  stop(
    "Only ", nrow(iam_markers),
    " IAM markers passed the prespecified filters; 200 are required for Fig. S4B."
  )
}

mouse_signature_list <- lapply(signature_sizes, function(n) {
  head(iam_markers$gene, n)
})
names(mouse_signature_list) <- paste0("MouseIAM", signature_sizes)

# Human top-specificity genes represented by mouse orthologues in the object.
# This is the exact list retained in the supplied Figure_4.R.
human_iam_signature <- c(
  "Acod1", "Irf1", "Stat1", "Aqp9", "Vcan", "S100a8", "Acsl1",
  "Trem1", "Ccr1", "Fpr1", "S100a9", "Cd300e", "S100a10", "S100a4"
)
human_iam_signature <- intersect(human_iam_signature, rownames(analysis_object))
write_csv(
  tibble(gene = human_iam_signature),
  file.path(OUTPUT_DIR, "Fig4B_human_IAM_signature_mouse_orthologues.csv")
)

# -----------------------------------------------------------------------------
# Figure 4B: mouse and human IAM module scores
# -----------------------------------------------------------------------------

analysis_object <- AddModuleScore(
  analysis_object,
  features = list(mouse_signature_list[["MouseIAM20"]]),
  name = "Mouse_IAM_",
  ctrl = 5,
  seed = 1
)
analysis_object <- AddModuleScore(
  analysis_object,
  features = list(human_iam_signature),
  name = "Human_IAM_",
  ctrl = 5,
  seed = 1
)

p_fig4b_mouse <- FeaturePlot(
  analysis_object,
  features = "Mouse_IAM_1",
  order = TRUE,
  raster = TRUE
) +
  ggtitle("Mouse Cxcl9/10 signature") +
  theme_void()

p_fig4b_human <- FeaturePlot(
  analysis_object,
  features = "Human_IAM_1",
  order = TRUE,
  raster = TRUE
) +
  ggtitle("Human CXCL9/10 signature") +
  theme_void()

save_panel(
  "Fig4B_mouse_human_IAM_module_scores.pdf",
  p_fig4b_mouse + p_fig4b_human,
  width = 7.4,
  height = 3.5
)

# -----------------------------------------------------------------------------
# Figure S4B: top-20/50/100/200 signatures and UCell concordance
# -----------------------------------------------------------------------------

DefaultAssay(analysis_object) <- "SCT"
analysis_object <- AddModuleScore(
  analysis_object,
  features = mouse_signature_list,
  name = "C8_sig",
  nbin = 24,
  ctrl = 100,
  seed = 1
)

temporary_addmodule_columns <- paste0("C8_sig", seq_along(signature_sizes))
final_addmodule_columns <- paste0("C8_sig", signature_sizes)
column_indices <- match(
  temporary_addmodule_columns,
  colnames(analysis_object[[]])
)
if (anyNA(column_indices)) stop("Could not identify all AddModuleScore columns.")
colnames(analysis_object@meta.data)[column_indices] <- final_addmodule_columns

DefaultAssay(analysis_object) <- "RNA"
ucell_scores <- UCell::ScoreSignatures_UCell(
  mat = GetAssayData(analysis_object, assay = "RNA", layer = "data"),
  features = mouse_signature_list
)

expected_ucell_columns <- paste0(names(mouse_signature_list), "_UCell")
missing_ucell_columns <- setdiff(expected_ucell_columns, colnames(ucell_scores))
if (length(missing_ucell_columns) > 0) {
  stop(
    "UCell did not return expected columns: ",
    paste(missing_ucell_columns, collapse = ", ")
  )
}
analysis_object@meta.data[, expected_ucell_columns] <- ucell_scores[, expected_ucell_columns]

p_s4b_signatures <- FeaturePlot(
  analysis_object,
  features = final_addmodule_columns,
  ncol = 2,
  min.cutoff = "q05",
  order = TRUE,
  raster = TRUE
) &
  theme_void()

correlation_matrix <- outer(
  final_addmodule_columns,
  expected_ucell_columns,
  Vectorize(function(addmodule_column, ucell_column) {
    cor(
      analysis_object[[]][[addmodule_column]],
      analysis_object[[]][[ucell_column]],
      method = "pearson",
      use = "complete.obs"
    )
  })
)
rownames(correlation_matrix) <- signature_sizes
colnames(correlation_matrix) <- signature_sizes

s4b_correlation_df <- as.data.frame(as.table(correlation_matrix)) |>
  setNames(c("AddModuleScore_size", "UCell_size", "Pearson_r"))
write_csv(
  s4b_correlation_df,
  file.path(OUTPUT_DIR, "FigS4B_AddModuleScore_UCell_correlations.csv")
)

p_s4b_correlations <- ggplot(
  s4b_correlation_df,
  aes(x = AddModuleScore_size, y = UCell_size, fill = Pearson_r)
) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", Pearson_r)), size = 3.5) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), oob = squish) +
  coord_fixed() +
  labs(x = "AddModuleScore size", y = "UCell size", fill = "Pearson r") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank())

save_panel(
  "FigS4B_signature_sizes_and_concordance.pdf",
  p_s4b_signatures + p_s4b_correlations,
  width = 10,
  height = 7
)

# -----------------------------------------------------------------------------
# Reproducible Monocle3 trajectory helpers
# -----------------------------------------------------------------------------

select_root_cells <- function(object) {
  metadata <- object[[]]
  candidates <- rownames(metadata)
  if ("tissue" %in% colnames(metadata) && any(metadata$tissue == "Blood")) {
    candidates <- rownames(metadata)[metadata$tissue == "Blood"]
  }
  ly6c2 <- FetchData(object, vars = "Ly6c2", cells = candidates, layer = "data")[[1]]
  n_root <- max(20L, ceiling(length(candidates) * 0.05))
  names(sort(ly6c2, decreasing = TRUE))[seq_len(min(n_root, length(ly6c2)))]
}

principal_node_for_cells <- function(cds, root_cells) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  root_cells <- intersect(root_cells, rownames(closest_vertex))
  if (length(root_cells) == 0) stop("No trajectory root cells were retained in the CDS.")

  root_vertex_index <- as.numeric(
    names(which.max(table(closest_vertex[root_cells, 1])))
  )
  igraph::V(principal_graph(cds)[["UMAP"]])$name[root_vertex_index]
}

run_trajectory <- function(object, graph_resolution) {
  DefaultAssay(object) <- "RNA"
  object[["RNA"]] <- JoinLayers(object[["RNA"]])
  root_cells <- select_root_cells(object)

  cds <- SeuratWrappers::as.cell_data_set(object)
  if (!"gene_short_name" %in% colnames(rowData(cds))) {
    rowData(cds)$gene_short_name <- rownames(rowData(cds))
  }
  cds <- cluster_cells(
    cds,
    reduction_method = "UMAP",
    resolution = graph_resolution
  )
  cds <- learn_graph(cds, use_partition = TRUE)
  cds <- order_cells(
    cds,
    reduction_method = "UMAP",
    root_pr_nodes = principal_node_for_cells(cds, root_cells)
  )
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
cds_naive <- run_trajectory(trajectory_objects[["Naive"]], graph_resolution = 0.00078)
cds_dss <- run_trajectory(trajectory_objects[["DSS"]], graph_resolution = 0.003)

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
p_fig4e <- FeaturePlot(
  analysis_object,
  features = fig4e_genes,
  ncol = 4,
  order = TRUE,
  raster = TRUE
) &
  theme_void()
save_panel("Fig4E_trajectory_gene_features.pdf", p_fig4e, width = 12, height = 3.2)

fig4f_genes <- c("Il1b", "Acod1")
missing_fig4f_genes <- setdiff(
  fig4f_genes,
  rowData(cds_dss)$gene_short_name
)
if (length(missing_fig4f_genes) > 0) {
  stop("Trajectory genes not found: ", paste(missing_fig4f_genes, collapse = ", "))
}
cds_fig4f <- cds_dss[rowData(cds_dss)$gene_short_name %in% fig4f_genes, ]
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
  graph_resolution = 0.00078
)
cds_colon_dss <- run_trajectory(
  colon_trajectory_objects[["DSS"]],
  graph_resolution = 0.003
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

# -----------------------------------------------------------------------------
# Figure S4F: Slamf7 over pseudotime in the combined DSS blood/colon trajectory
# -----------------------------------------------------------------------------

if (!"Slamf7" %in% rowData(cds_dss)$gene_short_name) {
  stop("Slamf7 was not found in the DSS trajectory object.")
}
cds_slamf7 <- cds_dss[rowData(cds_dss)$gene_short_name == "Slamf7", ]
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
