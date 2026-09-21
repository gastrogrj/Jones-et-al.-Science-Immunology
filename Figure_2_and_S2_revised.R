# Jones et al. - Figure 2 and Supplementary Figure 2
#
# This script contains only the murine scRNA-seq analyses used in Fig. 2 and
# Fig. S2. Flow-cytometry panels (Fig. 2G-I; Fig. S2A-B and S2G) were generated
# in FlowJo/Prism and are not recreated here. Fig. S2E was generated with the
# supplied Figure_S2E_scanpy.py, followed by Figure_S2E_concordance.R.
#
# Run from the repository root. All paths are relative and all generated files
# are written beneath output/Figure_2_S2.

suppressPackageStartupMessages({
  library(DoubletFinder)
  library(dplyr)
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
  library(harmony)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
  library(Seurat)
  library(SoupX)
  library(stringr)
  library(tibble)
  library(tidyr)
})

source("mouse_analysis_helpers.R")
set.seed(1)
options(future.globals.maxSize = 30 * 1024^3)

# -----------------------------------------------------------------------------
# Inputs, analysis constants and outputs
# -----------------------------------------------------------------------------

NAIVE_10X_DIR <- file.path("data", "mouse_scRNAseq", "G1")
DSS_10X_DIR <- file.path("data", "mouse_scRNAseq", "G2")

# To reproduce existing figures, use the final original Seurat objects.
# Raw reprocessing is used only when these files are absent; see README.
FINAL_COLON_OBJECT <- file.path("data","objects","colon_final.rds")
FULL_COLON_OBJECT <- file.path("data","objects","colon_full_postQC.rds")

MOUSE_REACTOME_GMT <- file.path(
  "data", "gmt", "m2.cp.reactome.v2025.1.Mm.symbols.gmt"
)

OUTPUT_DIR <- file.path("output", "Figure_2_S2")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

FINAL_RESOLUTION <- 0.32
PARENT_CSF1R_CLUSTERS <- c("0", "1")
IAM_FIGURE_CLUSTER <- "4"       # Seurat cluster 3, displayed as cluster 4
COMPARATOR_FIGURE_CLUSTER <- "2" # Seurat cluster 1, displayed as cluster 2

# The input to DoubletFinder is RNA log-normalised, so sct must be FALSE.
# This correction can change singlet calls. The historical merged-library
# workflow is retained for provenance; validate it against the archived cells.
# Per-library doublet detection is a separate sensitivity analysis, not silently
# substituted into this reproduction script.
DOUBLET_FINDER_USES_SCT <- FALSE

FIGURE_CLUSTER_LEVELS <- as.character(1:5)
# Cluster colours and raw-to-displayed labels are loaded from the shared helper.

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
# SoupX and Seurat preprocessing used for Figure 2/S2
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

make_sample_object <- function(counts, sample_name, timepoint) {
  object <- CreateSeuratObject(counts, project = sample_name)
  object$sample <- sample_name
  object$timepoint <- timepoint
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^Rp[Sl]")
  object[["percent.Xist"]] <- PercentageFeatureSet(object, pattern = "^Xist")
  subset(object, subset = nFeature_RNA > 500 & percent.mt < 10)
}

if(file.exists(FINAL_COLON_OBJECT)) {
  combined_csf1r <- label_colon(join_rna(readRDS(FINAL_COLON_OBJECT),normalise=TRUE))
  combined_all <- if(file.exists(FULL_COLON_OBJECT)) join_rna(readRDS(FULL_COLON_OBJECT)) else NULL
} else {
  warning("Raw preprocessing will rerun cell calls and clustering; compare cell identities with the manuscript before replacing figures.")
check_10x_directory(NAIVE_10X_DIR)
check_10x_directory(DSS_10X_DIR)

naive <- make_sample_object(
  run_soupx(NAIVE_10X_DIR),
  sample_name = "Naive_1",
  timepoint = "Naive"
)
dss <- make_sample_object(
  run_soupx(DSS_10X_DIR),
  sample_name = "DSS_1",
  timepoint = "DSS"
)

combined_all <- merge(
  naive,
  y = dss,
  add.cell.ids = c("Naive_1", "DSS_1"),
  project = "COMBINED"
)

sample_objects <- SplitObject(combined_all, split.by = "sample")
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
# DoubletFinder's RNA branch reads the merged object's Seurat command log.
# Record variable-feature settings here, then restore the original shared set.
combined_all <- join_rna(combined_all, normalise = TRUE)
combined_all <- FindVariableFeatures(combined_all, selection.method = "vst",
                                    nfeatures = 4000, verbose = FALSE)
VariableFeatures(combined_all) <- integration_features

combined_all <- combined_all |>
  ScaleData(features = integration_features, verbose = FALSE) |>
  RunPCA(features = integration_features, verbose = FALSE) |>
  RunHarmony(group.by.vars = "sample", verbose = FALSE) |>
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = FALSE) |>
  FindClusters(resolution = 0.1, verbose = FALSE)

combined_all[["RNA"]] <- JoinLayers(combined_all[["RNA"]])

# Retain historical merged-library parameters; correct SCT flag for RNA input.
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

# Csf1r+ mononuclear phagocytes were parent clusters 0 and 1 in the submitted
# analysis. The historical CSV write/read round-trip is unnecessary.
Idents(combined_all) <- "seurat_clusters"
combined_csf1r <- subset(combined_all, idents = PARENT_CSF1R_CLUSTERS)

combined_csf1r <- combined_csf1r |>
  SCTransform(
    vars.to.regress = c("percent.mt", "percent.ribo", "percent.Xist"),
    return.only.var.genes = FALSE,
    verbose = FALSE
  ) |>
  RunPCA(verbose = FALSE) |>
  RunHarmony(
    assay.use = "SCT",
    group.by.vars = "sample",
    verbose = FALSE
  ) |>
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindClusters(resolution = FINAL_RESOLUTION, verbose = FALSE)

combined_csf1r <- label_colon(join_rna(combined_csf1r,normalise=TRUE))

saveRDS(
  combined_csf1r,
  file.path(OUTPUT_DIR, "combined_colon_Csf1r_processed.rds"),
  compress = "xz"
)


saveRDS(combined_all,file.path(OUTPUT_DIR,"colon_full_postQC.rds"),compress="xz")

}
if(!is.null(combined_all)) {
# Figure S2C uses the post-QC full-cell object before Csf1r+ enrichment.
DefaultAssay(combined_all) <- "RNA"
p_s2c_clusters <- DimPlot(combined_all, label = FALSE, raster = TRUE) + NoLegend()
p_s2c_csf1r <- feature_grid(combined_all,"Csf1r",ncol=1)
save_panel(
  "FigS2C_full_UMAP_and_Csf1r.pdf",
  p_s2c_clusters + p_s2c_csf1r,
  width = 7.2,
  height = 3.4
)


} else message("S2C needs data/objects/colon_full_postQC.rds; other panels continue.")
saveRDS(combined_csf1r,file.path(OUTPUT_DIR,"combined_colon_Csf1r_processed.rds"),compress="xz")
export_cell_manifest(combined_csf1r,file.path(OUTPUT_DIR,"colon_cell_manifest.csv"))

# -----------------------------------------------------------------------------
# Figure 2A-B: UMAP and defining features
# -----------------------------------------------------------------------------

p_fig2a <- DimPlot(
  combined_csf1r,
  group.by = "figure_cluster",
  split.by = "timepoint",
  cols = FIGURE_CLUSTER_COLOURS,
  raster = TRUE,
  pt.size = 0.35
) +
  NoLegend() +
  theme_void()

save_panel("Fig2A_Csf1r_UMAP.pdf", p_fig2a, width = 7.4, height = 3.5)

fig2b_features <- c(
  "Ly6c2", "Itgal", "H2-Aa", "Ccr2",
  "Mrc1", "C1qa", "Itgax", "Cx3cr1"
)
stopifnot(all(fig2b_features %in% rownames(combined_csf1r)))

p_fig2b <- feature_grid(combined_csf1r,fig2b_features,ncol=2)
save_figure(file.path(OUTPUT_DIR,"Fig2B_cluster_defining_features"),p_fig2b,70,85)

# -----------------------------------------------------------------------------
# Figure 2C-D: top cluster markers and proportions
# -----------------------------------------------------------------------------

Idents(combined_csf1r) <- "figure_cluster"
all_markers <- FindAllMarkers(
  combined_csf1r,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0,
  test.use = "wilcox"
) |>
  as_tibble()

logfc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(all_markers))[1]
if (is.na(logfc_column)) stop("No log-fold-change column returned by FindAllMarkers.")

top10_markers <- all_markers |>
  filter(p_val_adj < 0.05) |>
  group_by(cluster) |>
  arrange(desc(.data[[logfc_column]]), .by_group = TRUE) |>
  slice_head(n = 10) |>
  ungroup()

write_csv(
  top10_markers,
  file.path(OUTPUT_DIR, "Fig2C_top10_markers_by_cluster.csv")
)

# If the final panel used a saved marker order, supply gene + cluster columns.
gene_order_file <- "data/annotations/Fig2C_final_genes.csv"
if(file.exists(gene_order_file)) {
  top10_markers <- read_csv(gene_order_file,show_col_types=FALSE)
  stopifnot(all(c("gene","cluster")%in%names(top10_markers)))
  check_features(combined_csf1r,top10_markers$gene)
}
write_csv(top10_markers,file.path(OUTPUT_DIR,"Fig2C_displayed_genes.csv"))
combined_csf1r <- ScaleData(
  combined_csf1r,
  assay = "RNA",
  features = unique(top10_markers$gene),
  verbose = FALSE
)

p_fig2c <- DoHeatmap(
  combined_csf1r,
  assay = "RNA",
  features = unique(top10_markers$gene),
  group.by = "figure_cluster",
  disp.min = -2.5, disp.max = 2.5,
  group.colors = unname(FIGURE_CLUSTER_COLOURS[as.character(1:5)]),
  raster = TRUE
) +
  guides(colour="none") +
  labs(fill="Scaled RNA expression (z-score)") +
  theme(axis.text.y=element_text(face="italic",size=7),legend.position="bottom",
        legend.title=element_text(size=7),legend.text=element_text(size=7))

save_figure(file.path(OUTPUT_DIR,"Fig2C_top10_marker_heatmap"),p_fig2c,90,110)
write_csv(all_markers,file.path(OUTPUT_DIR,"Fig2C_all_cluster_markers.csv"))

cluster_proportions <- combined_csf1r[[]] |>
  count(timepoint, figure_cluster, name = "n_cells") |>
  group_by(timepoint) |>
  mutate(percent_of_condition = 100 * n_cells / sum(n_cells)) |>
  ungroup()

write_csv(
  cluster_proportions,
  file.path(OUTPUT_DIR, "Fig2D_cluster_proportions.csv")
)

p_fig2d <- ggplot(
  cluster_proportions,
  aes(x = figure_cluster, y = percent_of_condition, fill = figure_cluster)
) +
  geom_col(width = 0.76) +
  facet_wrap(~timepoint, ncol = 1) +
  scale_fill_manual(values = FIGURE_CLUSTER_COLOURS, guide = "none") +
  labs(x = "Cluster", y = "Proportion of sample (%)") +
  theme_classic(base_size = 10)

save_panel("Fig2D_cluster_proportions.pdf", p_fig2d, width = 4.2, height = 5.2)

# -----------------------------------------------------------------------------
# Figure 2E: IAM marker dot plot (cluster 4 versus all other clusters)
# -----------------------------------------------------------------------------

fig2e_genes <- c("Ccrl2", "Cd274", "Slamf7")
fig2e_data <- FetchData(
  combined_csf1r,
  vars = c(fig2e_genes, "figure_cluster", "timepoint"), layer="data"
) |>
  filter(timepoint == "DSS") |>
  mutate(
    group = if_else(
      as.character(figure_cluster) == IAM_FIGURE_CLUSTER,
      "Cluster 4",
      "Other clusters"
    )
  ) |>
  pivot_longer(all_of(fig2e_genes), names_to = "gene", values_to = "expression") |>
  group_by(group, gene) |>
  summarise(
    percent_expressed = 100 * mean(expression > 0),
    average_expression = mean(expression),
    .groups = "drop"
  ) |>
  group_by(gene) |>
  mutate(average_expression_scaled = scale_or_zero(average_expression)) |>
  ungroup() |>
  mutate(group = factor(group, levels = c("Cluster 4", "Other clusters")))

p_fig2e <- ggplot(
  fig2e_data,
  aes(x = gene, y = group, size = percent_expressed, fill = average_expression_scaled)
) +
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(2, 9), name = "% expressed") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Average\nexpression\n(scaled)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1))

save_panel("Fig2E_IAM_marker_dotplot.pdf", p_fig2e, width = 4.2, height = 3.2)

# -----------------------------------------------------------------------------
# Figure 2F: cell-level Spearman correlations in cluster 4
# -----------------------------------------------------------------------------

correlation_x <- c("Cxcl9", "Cxcl10")
correlation_y <- c("Slamf7", "Cd274", "Ccrl2")

iam_cells <- WhichCells(combined_csf1r, idents = IAM_FIGURE_CLUSTER)
correlation_expression <- FetchData(
  combined_csf1r,
  vars = c(correlation_x, correlation_y),
  cells = iam_cells,
  layer = "data"
)

fig2f_correlations <- tidyr::crossing(
  gene_x = correlation_x,
  gene_y = correlation_y
) |>
  mutate(statistics = map2(
    gene_x,
    gene_y,
    function(x, y) safe_spearman(
      correlation_expression[[x]],
      correlation_expression[[y]]
    )
  )) |>
  unnest(statistics) |>
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH"),
    significance = if_else(p_value < 0.05, "*", ""),
    gene_x = factor(gene_x, levels = correlation_x),
    gene_y = factor(gene_y, levels = rev(correlation_y))
  )

write_csv(
  fig2f_correlations,
  file.path(OUTPUT_DIR, "Fig2F_cluster4_Spearman_correlations.csv")
)

p_fig2f <- ggplot(
  fig2f_correlations,
  aes(x = gene_x, y = gene_y, fill = rho)
) +
  geom_tile(colour = "white") +
  geom_text(aes(label = significance), size = 5) +
  scale_fill_gradient(limits = c(0, 0.5), low = "white", high = "red", oob = squish) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = expression(rho)) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic")
  )

save_panel("Fig2F_cluster4_Spearman_heatmap.pdf", p_fig2f, width = 3.6, height = 3.6)

# -----------------------------------------------------------------------------
# Figure S2D: stability of the Cxcl9/10+ cluster over clustering resolutions
# -----------------------------------------------------------------------------

STABILITY_GRAPH <- "SCT_snn"
if (!STABILITY_GRAPH %in% names(combined_csf1r@graphs))
  stop("S2D needs the original SCT_snn graph in colon_final.rds; do not substitute a newly constructed graph.")
reference_object <- FindClusters(combined_csf1r,graph.name=STABILITY_GRAPH,
                                resolution=0.25,verbose=FALSE)
reference_members <- split(colnames(reference_object),as.character(Idents(reference_object)))
overlap <- vapply(reference_members,function(z) length(intersect(z,iam_cells))/length(union(z,iam_cells)),numeric(1))
reference_iam_cells <- reference_members[[which.max(overlap)]]
resolution_values <- c(0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)

resolution_stability <- map_dfr(resolution_values, function(resolution) {
  object_at_resolution <- FindClusters(
    combined_csf1r,
    graph.name = STABILITY_GRAPH,
    resolution = resolution,
    verbose = FALSE
  )
  cluster_cells <- split(
    colnames(object_at_resolution),
    as.character(Idents(object_at_resolution))
  )
  jaccard <- vapply(cluster_cells, function(cells) {
    length(intersect(cells, reference_iam_cells)) /
      length(union(cells, reference_iam_cells))
  }, numeric(1))
  best_cluster <- names(which.max(jaccard))
  selected_cells <- cluster_cells[[best_cluster]]

  tibble(
    reference_resolution = 0.25,
    reference_cluster = names(reference_members)[which.max(overlap)],
    selected_cluster = best_cluster,
    purity_percent = 100*length(intersect(selected_cells,reference_iam_cells))/length(selected_cells),
    retention_percent = 100*length(intersect(selected_cells,reference_iam_cells))/length(reference_iam_cells),
    jaccard = max(jaccard),
    resolution = resolution,
    identity = c("Cxcl9/10+ cells", "All other cells"),
    n_cells = c(
      length(intersect(selected_cells, reference_iam_cells)),
      length(setdiff(selected_cells, reference_iam_cells))
    )
  ) |>
    mutate(percent = 100 * n_cells / sum(n_cells))
})

write_csv(
  resolution_stability,
  file.path(OUTPUT_DIR, "FigS2D_cluster_resolution_stability.csv")
)

p_s2d <- ggplot(
  resolution_stability,
  aes(x = factor(resolution), y = percent, fill = identity)
) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = c("Cxcl9/10+ cells" = "#2C7FB8", "All other cells" = "#E43D30")) +
  labs(x = "Resolution", y = "% of selected cluster", fill = NULL) +
  theme_classic(base_size = 9)

save_panel("FigS2D_cluster_resolution_stability.pdf", p_s2d, width = 5.8, height = 3.4)

# S2F and S2I Slamf7 tests/plots are generated together after Cd274 grouping.
# S2H is generated only by Figure_S2H_revised.R.

# -----------------------------------------------------------------------------
# Figure S2I: Cd274-high versus Cd274-low DSS monocytes
# -----------------------------------------------------------------------------

DefaultAssay(combined_csf1r) <- "RNA"
s2i_metadata <- combined_csf1r[[]]
s2i_cells <- rownames(s2i_metadata)[
  s2i_metadata$timepoint == "DSS" &
    as.character(s2i_metadata$figure_cluster) %in%
      c(IAM_FIGURE_CLUSTER, COMPARATOR_FIGURE_CLUSTER)
]

cd274_data <- FetchData(
  combined_csf1r,
  vars = "Cd274",
  cells = s2i_cells,
  layer = "data"
)
cd274_values <- setNames(cd274_data[["Cd274"]], rownames(cd274_data))

lower_quartile <- unname(quantile(cd274_values, 0.25, na.rm = TRUE))
upper_quartile <- unname(quantile(cd274_values, 0.75, na.rm = TRUE))
cd274_low_cells <- names(cd274_values)[cd274_values <= lower_quartile]
cd274_high_cells <- names(cd274_values)[cd274_values >= upper_quartile]
stopifnot(length(intersect(cd274_low_cells, cd274_high_cells)) == 0)

cd274_group <- rep(NA_character_, ncol(combined_csf1r))
names(cd274_group) <- colnames(combined_csf1r)
cd274_group[cd274_low_cells] <- "Cd274Lo"
cd274_group[cd274_high_cells] <- "Cd274Hi"
combined_csf1r$Cd274_group <- factor(
  cd274_group[colnames(combined_csf1r)],
  levels = c("Cd274Lo", "Cd274Hi")
)

s2i_de <- FindMarkers(
  combined_csf1r,
  ident.1 = cd274_high_cells,
  ident.2 = cd274_low_cells,
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
) |>
  rownames_to_column("gene")

s2i_logfc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(s2i_de))[1]
if (is.na(s2i_logfc_column)) stop("No log-fold-change column returned for Fig. S2I.")
s2i_de <- s2i_de |>
  mutate(
    log2FC = .data[[s2i_logfc_column]],
    minus_log10_padj = -log10(pmax(p_val_adj, 1e-300))
  )

write_csv(s2i_de, file.path(OUTPUT_DIR, "FigS2I_Cd274Hi_vs_Lo_DE.csv"))

# Balanced cell-level heatmap: 15 significant genes in each direction.
s2i_significant <- s2i_de |> filter(p_val_adj < 0.05, is.finite(log2FC))
s2i_low_genes <- s2i_significant |>
  filter(log2FC < 0) |>
  arrange(p_val_adj, log2FC) |>
  slice_head(n = 15) |>
  pull(gene)
s2i_high_genes <- s2i_significant |>
  filter(log2FC > 0) |>
  arrange(p_val_adj, desc(log2FC)) |>
  slice_head(n = 15) |>
  pull(gene)
s2i_heatmap_genes <- unique(c(s2i_low_genes, s2i_high_genes))

set.seed(1)
s2i_heatmap_cells <- c(
  sample(cd274_low_cells, min(200, length(cd274_low_cells))),
  sample(cd274_high_cells, min(200, length(cd274_high_cells)))
)

combined_csf1r <- ScaleData(
  combined_csf1r,
  assay = "RNA",
  features = s2i_heatmap_genes,
  verbose = FALSE
)
p_s2i_heatmap <- DoHeatmap(
  combined_csf1r,
  assay = "RNA",
  features = s2i_heatmap_genes,
  cells = s2i_heatmap_cells,
  group.by = "Cd274_group",
  disp.min = -2.5,
  disp.max = 2.5,
  raster = TRUE
) +
  guides(colour="none") +
  labs(fill="Scaled RNA expression (z-score)") +
  theme(axis.text.y = element_text(face = "italic", size = 7)) +
  labs(caption=paste0("Analysed cells: Cd274Hi ",length(cd274_high_cells),", Cd274Lo ",length(cd274_low_cells),
       "; displayed cells: ",sum(s2i_heatmap_cells%in%cd274_high_cells)," and ",sum(s2i_heatmap_cells%in%cd274_low_cells),"."))

save_panel("FigS2I_Cd274Hi_vs_Lo_heatmap.pdf", p_s2i_heatmap, width = 6.5, height = 7.0)

# Bidirectional volcano plot.
s2i_volcano <- s2i_de |>
  mutate(
    direction = case_when(
      p_val_adj < 0.05 & log2FC >= 0.25 ~ "Up in Cd274Hi",
      p_val_adj < 0.05 & log2FC <= -0.25 ~ "Up in Cd274Lo",
      TRUE ~ "Not significant"
    )
  )
s2i_labels <- bind_rows(
  s2i_volcano |>
    filter(direction == "Up in Cd274Hi") |>
    arrange(p_val_adj, desc(log2FC)) |>
    slice_head(n = 10),
  s2i_volcano |>
    filter(direction == "Up in Cd274Lo") |>
    arrange(p_val_adj, log2FC) |>
    slice_head(n = 10)
)

p_s2i_volcano <- ggplot(
  s2i_volcano,
  aes(x = log2FC, y = minus_log10_padj, colour = direction)
) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  ggrepel::geom_text_repel(data = s2i_labels, aes(label = gene), size = 2.8) +
  scale_colour_manual(values = c(
    "Up in Cd274Lo" = "#2166AC",
    "Not significant" = "grey75",
    "Up in Cd274Hi" = "#B2182B"
  )) +
  labs(x = "log2 fold change (Cd274Hi / Cd274Lo)", y = "-log10 adjusted P", colour = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top")

save_panel("FigS2I_Cd274Hi_vs_Lo_volcano.pdf", p_s2i_volcano, width = 6.0, height = 5.0)

# Two targeted Slamf7 comparisons form the prespecified Holm family, matching
# the current supplementary panels F and I. Use Seurat raw Wilcoxon p values.
groups <- list(
  FigS2F=list(hi=cell_groups(combined_csf1r,"4","DSS"),lo=cell_groups(combined_csf1r,"2","DSS"),
              labels=c("Cluster 4","Cluster 2")),
  FigS2I=list(hi=cd274_high_cells,lo=cd274_low_cells,labels=c("Cd274Hi","Cd274Lo")))
slamf7_tests <- bind_rows(lapply(names(groups),function(panel) {
  g<-groups[[panel]]
  tt<-FindMarkers(combined_csf1r,ident.1=g$hi,ident.2=g$lo,features="Slamf7",assay="RNA",
                  test.use="wilcox",logfc.threshold=0,min.pct=0)
  tibble(panel=panel,p_raw=tt["Slamf7","p_val"],n_high=length(g$hi),n_low=length(g$lo))
})) |> mutate(p_holm=p.adjust(p_raw,method="holm"))
write_csv(slamf7_tests,file.path(OUTPUT_DIR,"FigS2F_I_Slamf7_Holm_tests.csv"))
for(panel in names(groups)) {
  g<-groups[[panel]]
  dd<-FetchData(combined_csf1r,vars=c("Slamf7","Cd274"),cells=c(g$hi,g$lo),layer="data")
  dd$group<-factor(ifelse(rownames(dd)%in%g$hi,g$labels[1],g$labels[2]),levels=g$labels)
  p<-ggplot(dd,aes(group,Slamf7,fill=group))+geom_violin(scale="width")+
    geom_boxplot(width=.12,outlier.shape=NA,fill="white")+
    stat_summary(fun=mean,geom="point",shape=18,size=2,colour="#F0D000")+
    scale_fill_manual(values=c("#B2182B","#2166AC"),guide="none")+
    labs(x=NULL,y="Slamf7 RNA log1p(counts per 10,000)",
         subtitle=paste0("Holm-adjusted P = ",format.pval(slamf7_tests$p_holm[slamf7_tests$panel==panel],digits=3)))+
    theme_classic(base_size=7)
  save_figure(file.path(OUTPUT_DIR,paste0(panel,"_Slamf7_violin")),p,60,65)
  avg<-dd |> pivot_longer(c(Slamf7,Cd274),names_to="gene",values_to="expression") |>
    group_by(group,gene) |> summarise(mean_expression=mean(expression),percent=100*mean(expression>0),.groups="drop") |>
    group_by(gene) |> mutate(z=scale_or_zero(mean_expression)) |> ungroup()
  p<-ggplot(avg,aes(group,gene,fill=z))+geom_tile(colour="white")+
    geom_text(aes(label=sprintf("%.2f\n%.0f%%+",mean_expression,percent)),size=7/ggplot2::.pt)+
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-1,1),name="Row-scaled mean")+
    labs(x=NULL,y=NULL)+theme_minimal(base_size=7)+theme(panel.grid=element_blank())
  save_figure(file.path(OUTPUT_DIR,paste0(panel,"_Slamf7_Cd274_heatmap")),p,65,45)
  write_csv(avg,file.path(OUTPUT_DIR,paste0(panel,"_Slamf7_Cd274_summary.csv")))
}

# Bidirectional Reactome GSEA.
if (file.exists(MOUSE_REACTOME_GMT)) {
  mouse_reactome <- fgsea::gmtPathways(MOUSE_REACTOME_GMT)
  s2i_ranks <- s2i_de |>
    filter(is.finite(log2FC), !is.na(gene), gene != "") |>
    arrange(desc(abs(log2FC))) |>
    distinct(gene, .keep_all = TRUE) |>
    filter(!str_detect(gene, "^(Igh|Igk|Igl|Jchain|Ms4a1|Cd79a|Cd79b)")) |>
    select(gene, log2FC) |>
    deframe() |>
    sort(decreasing = TRUE)

  s2i_fgsea <- fgsea::fgsea(
    pathways = mouse_reactome,
    stats = s2i_ranks,
    minSize = 10,
    maxSize = 2000,
    eps = 0,
    nPermSimple = 10000
  ) |>
    as_tibble() |>
    arrange(padj, desc(abs(NES)))

  s2i_fgsea_export <- s2i_fgsea |>
    mutate(leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    ))
  write_csv(
    s2i_fgsea_export,
    file.path(OUTPUT_DIR, "FigS2I_Cd274Hi_vs_Lo_Reactome_FGSEA.csv")
  )

  s2i_fgsea_plot <- bind_rows(
    s2i_fgsea |>
      filter(padj < 0.05, NES < 0) |>
      arrange(NES) |>
      slice_head(n = 10),
    s2i_fgsea |>
      filter(padj < 0.05, NES > 0) |>
      arrange(desc(NES)) |>
      slice_head(n = 10)
  ) |>
    mutate(
      enrichment = if_else(NES > 0, "Cd274Hi", "Cd274Lo"),
      pathway_label = tools::toTitleCase(
        gsub("_", " ", sub("^REACTOME_", "", pathway))
      ),
      pathway_label = factor(pathway_label, levels = pathway_label[order(NES)])
    )

  p_s2i_fgsea <- ggplot(
    s2i_fgsea_plot,
    aes(x = NES, y = pathway_label, fill = enrichment)
  ) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("Cd274Lo" = "#2166AC", "Cd274Hi" = "#B2182B")) +
    labs(x = "Normalised enrichment score", y = NULL, fill = NULL) +
    theme_classic(base_size = 9) +
    theme(legend.position = "top")

  save_panel("FigS2I_Cd274Hi_vs_Lo_Reactome_FGSEA.pdf", p_s2i_fgsea, width = 8.4, height = 6.5)
} else {
  message("Fig. S2I FGSEA not run: missing ", MOUSE_REACTOME_GMT)
}

message(
  "Fig. S2E is intentionally not duplicated in R; it was generated by the ",
  "included Figure_S2E_scanpy.py then Figure_S2E_concordance.R."
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)

write_csv(tibble(group=c("Cd274Hi","Cd274Lo"),n=c(length(cd274_high_cells),length(cd274_low_cells)),
                 threshold=c(upper_quartile,lower_quartile)),file.path(OUTPUT_DIR,"FigS2I_group_sizes_and_thresholds.csv"))
# A zero lower quartile means Cd274Lo is undetected; ties can exceed 25%.
write_csv(tibble(cell=c(cd274_high_cells,cd274_low_cells),
                 group=c(rep("Cd274Hi",length(cd274_high_cells)),rep("Cd274Lo",length(cd274_low_cells)))),
          file.path(OUTPUT_DIR,"FigS2I_cell_groups.csv"))
