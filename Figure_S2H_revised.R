# Jones et al. - Supplementary Figure 2H
#
# Re-analysis of the independent Hegarty et al. recovery-colitis scRNA-seq
# dataset. This script contains only the processing and plotting used for
# Fig. S2H: naive/resolution UMAPs, the monocyte-macrophage marker dot plot,
# the resolution cluster-5 versus clusters-2+3 dot plot, and Spearman rank
# correlations between Cxcl9 and Cd274/Ccrl2/Slamf7 in resolution cluster 5.
#
# Run from the repository root:
#   Rscript code/Figure_S2H_revised.R
#
# Expected input layout:
#   data/Figure_S2H/R_control/{raw,filtered}_feature_bc_matrix/
#   data/Figure_S2H/R_resolution/{raw,filtered}_feature_bc_matrix/
#
# All generated files are written to output/Figure_S2H.

suppressPackageStartupMessages({
  library(DoubletFinder)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(SoupX)
  library(tidyr)
})

set.seed(1)
options(future.globals.maxSize = 30 * 1024^3)

# -----------------------------------------------------------------------------
# Inputs, parameters and outputs
# -----------------------------------------------------------------------------

DATA_DIR <- file.path("data", "Figure_S2H")
NAIVE_10X_DIR <- file.path(DATA_DIR, "R_control")
RESOLUTION_10X_DIR <- file.path(DATA_DIR, "R_resolution")
OUTPUT_DIR <- file.path("output", "Figure_S2H")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_FEATURES <- 500

# The supplied analysis used 15%. The manuscript's general mouse scRNA-seq
# methods state 10%; retain 15 here to track the analysis used for the displayed
# S2H panel, but reconcile this value before the public release.
MAX_PERCENT_MT <- 15

EXPECTED_DOUBLET_RATE <- 0.07
DOUBLET_PN <- 0.25
DOUBLET_PK <- 0.09

# Sequential cluster selections and resolutions used in the supplied analysis.
INITIAL_MNP_CLUSTERS <- c("0", "1")
FIRST_NON_DC_CLUSTERS <- c("0", "1", "2")
SECOND_NON_DC_CLUSTERS <- c("0", "1", "2")
FINAL_CLUSTER_LEVELS <- as.character(0:6)

MARKER_FEATURES <- c(
  "Csf1r", "C1qa", "Cd163", "Mrc1", "Apoe", "Itgax", "Cx3cr1",
  "Ccr2", "Ly6c2", "H2-Aa", "Acod1", "Cxcl9", "Cxcl10"
)

IAM_FEATURES <- c("Cd274", "Slamf7", "Ccrl2", "Il1b", "Tnf")
CORRELATION_TARGETS <- c("Cd274", "Ccrl2", "Slamf7")
RESOLUTION_CLUSTER_OF_INTEREST <- "5"
RESOLUTION_COMPARATOR_CLUSTERS <- c("2", "3")

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

check_features <- function(object, features) {
  missing <- setdiff(features, rownames(object))
  if (length(missing) > 0) {
    stop("Genes absent from the object: ", paste(missing, collapse = ", "))
  }
}

numeric_levels <- function(x) {
  values <- unique(as.character(x))
  numeric_values <- suppressWarnings(as.numeric(values))
  if (all(!is.na(numeric_values))) {
    values[order(numeric_values)]
  } else {
    sort(values)
  }
}

scale_or_zero <- function(x) {
  value_sd <- stats::sd(x, na.rm = TRUE)
  if (length(x) < 2 || is.na(value_sd) || value_sd == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

save_panel <- function(stem, plot, width, height) {
  ggsave(
    file.path(OUTPUT_DIR, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      file.path(OUTPUT_DIR, paste0(stem, ".svg")),
      plot,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
  }
}

# -----------------------------------------------------------------------------
# SoupX correction and initial Seurat processing
# -----------------------------------------------------------------------------

read_gene_expression <- function(path) {
  counts <- Seurat::Read10X(path)
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("No 'Gene Expression' matrix found in ", path)
    }
    counts <- counts[["Gene Expression"]]
  }
  counts
}

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
    RunUMAP(dims = 1:30, verbose = FALSE, seed.use = 1) |>
    FindNeighbors(dims = 1:30, verbose = FALSE) |>
    FindClusters(verbose = FALSE, random.seed = 1)

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

  keep_cells <- rownames(object[[]])[
    object$nFeature_RNA > MIN_FEATURES &
      object$percent.mt < MAX_PERCENT_MT
  ]
  subset(object, cells = keep_cells)
}

check_10x_directory(NAIVE_10X_DIR)
check_10x_directory(RESOLUTION_10X_DIR)

naive <- make_sample_object(
  run_soupx(NAIVE_10X_DIR),
  sample_name = "Naive_1",
  timepoint = "Naive"
)
resolution <- make_sample_object(
  run_soupx(RESOLUTION_10X_DIR),
  sample_name = "Resolution_1",
  timepoint = "Resolution"
)

combined_all <- merge(
  naive,
  y = resolution,
  add.cell.ids = c("Naive_1", "Resolution_1"),
  project = "Hegarty_recovery"
)

sample_objects <- SplitObject(combined_all, split.by = "sample")
sample_objects <- lapply(sample_objects, function(object) {
  object |>
    NormalizeData(verbose = FALSE) |>
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = 4000,
      verbose = FALSE
    )
})

# SelectIntegrationFeatures defaults to 2,000 shared features. This retains the
# behaviour of the supplied pipeline after identifying 4,000 features/sample.
integration_features <- SelectIntegrationFeatures(object.list = sample_objects)
VariableFeatures(combined_all) <- integration_features

combined_all <- combined_all |>
  NormalizeData(verbose = FALSE) |>
  ScaleData(features = integration_features, verbose = FALSE) |>
  RunPCA(features = integration_features, verbose = FALSE) |>
  RunHarmony(group.by.vars = "sample", verbose = FALSE) |>
  RunUMAP(
    reduction = "harmony",
    dims = 1:20,
    verbose = FALSE,
    seed.use = 1
  ) |>
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = FALSE) |>
  FindClusters(resolution = 0.1, verbose = FALSE, random.seed = 1)

if (inherits(combined_all[["RNA"]], "Assay5")) {
  combined_all[["RNA"]] <- SeuratObject::JoinLayers(combined_all[["RNA"]])
}

# The merged object has an RNA assay at this point, not an SCT assay; therefore
# DoubletFinder must be run with sct = FALSE. The supplied file used sct = TRUE,
# which is inconsistent with the object created immediately above.
doublet_finder_name <- if (
  "doubletFinder" %in% getNamespaceExports("DoubletFinder")
) {
  "doubletFinder"
} else if (
  "doubletFinder_v3" %in% getNamespaceExports("DoubletFinder")
) {
  "doubletFinder_v3"
} else {
  stop("No supported DoubletFinder function was found.")
}
doublet_finder <- getExportedValue("DoubletFinder", doublet_finder_name)

expected_doublets <- round(ncol(combined_all) * EXPECTED_DOUBLET_RATE)
combined_all <- doublet_finder(
  combined_all,
  PCs = 1:10,
  pN = DOUBLET_PN,
  pK = DOUBLET_PK,
  nExp = expected_doublets,
  sct = FALSE
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

# -----------------------------------------------------------------------------
# Sequential Csf1r+ enrichment and final reclustering used for Fig. S2H
# -----------------------------------------------------------------------------

Idents(combined_all) <- "seurat_clusters"
combined_mnps <- subset(combined_all, idents = INITIAL_MNP_CLUSTERS)

combined_mnps <- SCTransform(
  combined_mnps,
  vars.to.regress = c("percent.mt", "percent.ribo", "percent.Xist"),
  verbose = FALSE,
  return.only.var.genes = FALSE
) |>
  RunPCA(verbose = FALSE) |>
  RunHarmony(
    assay.use = "SCT",
    group.by.vars = "sample",
    verbose = FALSE
  ) |>
  RunUMAP(
    reduction = "harmony",
    dims = 1:30,
    verbose = FALSE,
    seed.use = 1
  ) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindClusters(resolution = 0.3, verbose = FALSE, random.seed = 1)

Idents(combined_mnps) <- "seurat_clusters"
combined_no_dcs <- subset(combined_mnps, idents = FIRST_NON_DC_CLUSTERS)
DefaultAssay(combined_no_dcs) <- "SCT"

combined_no_dcs <- combined_no_dcs |>
  RunPCA(verbose = FALSE) |>
  RunHarmony(
    assay.use = "SCT",
    group.by.vars = "sample",
    verbose = FALSE
  ) |>
  RunUMAP(
    reduction = "harmony",
    dims = 1:30,
    verbose = FALSE,
    seed.use = 1
  ) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindClusters(resolution = 0.3, verbose = FALSE, random.seed = 1)

Idents(combined_no_dcs) <- "seurat_clusters"
figure_object <- subset(combined_no_dcs, idents = SECOND_NON_DC_CLUSTERS)
DefaultAssay(figure_object) <- "SCT"

figure_object <- figure_object |>
  RunPCA(verbose = FALSE) |>
  RunHarmony(
    assay.use = "SCT",
    group.by.vars = "sample",
    verbose = FALSE
  ) |>
  RunUMAP(
    reduction = "harmony",
    dims = 1:30,
    verbose = FALSE,
    seed.use = 1
  ) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindClusters(resolution = 0.5, verbose = FALSE, random.seed = 1)

observed_clusters <- numeric_levels(figure_object$seurat_clusters)
if (!setequal(observed_clusters, FINAL_CLUSTER_LEVELS)) {
  warning(
    "Expected final clusters 0-6, but observed: ",
    paste(observed_clusters, collapse = ", "),
    ". Check package versions, seeds and upstream inputs."
  )
}

DefaultAssay(figure_object) <- "RNA"
if (inherits(figure_object[["RNA"]], "Assay5")) {
  figure_object[["RNA"]] <- SeuratObject::JoinLayers(figure_object[["RNA"]])
}
check_features(
  figure_object,
  unique(c(MARKER_FEATURES, IAM_FEATURES, "Cxcl9", CORRELATION_TARGETS))
)

# -----------------------------------------------------------------------------
# Naive and resolution UMAPs
# -----------------------------------------------------------------------------

umap_coordinates <- as.data.frame(Embeddings(figure_object, "umap"))
colnames(umap_coordinates)[1:2] <- c("UMAP_1", "UMAP_2")

umap_data <- cbind(
  figure_object[[]],
  umap_coordinates[, c("UMAP_1", "UMAP_2"), drop = FALSE]
) |>
  tibble::rownames_to_column("cell") |>
  mutate(
    timepoint = factor(timepoint, levels = c("Naive", "Resolution")),
    cluster = factor(
      as.character(seurat_clusters),
      levels = numeric_levels(seurat_clusters)
    )
  )

cluster_colours <- setNames(
  scales::hue_pal()(length(levels(umap_data$cluster))),
  levels(umap_data$cluster)
)

p_umap <- ggplot(umap_data, aes(UMAP_1, UMAP_2, fill = cluster)) +
  geom_point(shape = 21, size = 0.55, stroke = 0, alpha = 0.9) +
  facet_wrap(~timepoint, nrow = 1) +
  scale_fill_manual(values = cluster_colours, drop = FALSE) +
  coord_equal() +
  labs(x = "UMAP_2", y = "UMAP_1", fill = "Cluster") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    legend.title = element_blank()
  )

write.csv(
  umap_data |>
    select(cell, timepoint, cluster, UMAP_1, UMAP_2),
  file.path(OUTPUT_DIR, "FigS2H_UMAP_coordinates.csv"),
  row.names = FALSE
)
save_panel("FigS2H_UMAPs", p_umap, width = 6.2, height = 3.0)

# -----------------------------------------------------------------------------
# Marker dot plot across final clusters
# -----------------------------------------------------------------------------

cluster_levels <- numeric_levels(figure_object$seurat_clusters)
marker_dot_data <- FetchData(
  figure_object,
  vars = c(MARKER_FEATURES, "seurat_clusters"),
  layer = "data"
) |>
  rename(cluster = seurat_clusters) |>
  pivot_longer(
    cols = all_of(MARKER_FEATURES),
    names_to = "feature",
    values_to = "expression"
  ) |>
  group_by(cluster, feature) |>
  summarise(
    fraction_expressing = 100 * mean(expression > 0, na.rm = TRUE),
    mean_expression = mean(expression, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) |>
  mutate(
    cluster = factor(as.character(cluster), levels = rev(cluster_levels)),
    feature = factor(feature, levels = MARKER_FEATURES)
  )

p_marker_dot <- ggplot(
  marker_dot_data,
  aes(feature, cluster, size = fraction_expressing, fill = mean_expression)
) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.25) +
  scale_size_continuous(
    range = c(1, 7),
    limits = c(0, 100),
    breaks = c(20, 40, 60, 80, 100),
    name = "Fraction of cells\nin group (%)"
  ) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    name = "Mean expression\nin group"
  ) +
  labs(x = NULL, y = "Cluster") +
  theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "italic"
    )
  )

write.csv(
  marker_dot_data,
  file.path(OUTPUT_DIR, "FigS2H_marker_dotplot_data.csv"),
  row.names = FALSE
)
save_panel("FigS2H_marker_dotplot", p_marker_dot, width = 7.2, height = 3.5)

# -----------------------------------------------------------------------------
# Resolution cluster 5 versus pooled clusters 2 + 3
# -----------------------------------------------------------------------------

iam_dot_data <- FetchData(
  figure_object,
  vars = c(IAM_FEATURES, "seurat_clusters", "timepoint"),
  layer = "data"
) |>
  filter(
    timepoint == "Resolution",
    as.character(seurat_clusters) %in% c(
      RESOLUTION_CLUSTER_OF_INTEREST,
      RESOLUTION_COMPARATOR_CLUSTERS
    )
  ) |>
  mutate(
    group = case_when(
      as.character(seurat_clusters) == RESOLUTION_CLUSTER_OF_INTEREST ~
        "Cluster 5",
      as.character(seurat_clusters) %in% RESOLUTION_COMPARATOR_CLUSTERS ~
        "Clusters 2 + 3",
      TRUE ~ NA_character_
    )
  ) |>
  pivot_longer(
    cols = all_of(IAM_FEATURES),
    names_to = "feature",
    values_to = "expression"
  ) |>
  group_by(group, feature) |>
  summarise(
    percent_expressed = 100 * mean(expression > 0, na.rm = TRUE),
    average_expression = mean(expression, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) |>
  group_by(feature) |>
  mutate(average_expression_scaled = scale_or_zero(average_expression)) |>
  ungroup() |>
  mutate(
    group = factor(group, levels = c("Cluster 5", "Clusters 2 + 3")),
    feature = factor(feature, levels = IAM_FEATURES)
  )

expected_iam_groups <- c("Cluster 5", "Clusters 2 + 3")
observed_iam_groups <- unique(as.character(iam_dot_data$group))
if (!setequal(observed_iam_groups, expected_iam_groups)) {
  stop(
    "Missing Fig. S2H comparison group(s): ",
    paste(setdiff(expected_iam_groups, observed_iam_groups), collapse = ", ")
  )
}

p_iam_dot <- ggplot(
  iam_dot_data,
  aes(
    feature,
    group,
    size = percent_expressed,
    fill = average_expression_scaled
  )
) +
  geom_point(shape = 21, colour = "black", stroke = 0.3) +
  scale_size_continuous(
    range = c(1.5, 6),
    limits = c(0, 100),
    name = "Percent expressed"
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Average expression\n(scaled)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")
  )

write.csv(
  iam_dot_data,
  file.path(OUTPUT_DIR, "FigS2H_cluster5_vs_clusters2plus3_dotplot_data.csv"),
  row.names = FALSE
)
save_panel(
  "FigS2H_cluster5_vs_clusters2plus3_dotplot",
  p_iam_dot,
  width = 4.2,
  height = 3.0
)

# -----------------------------------------------------------------------------
# Spearman correlations: Cxcl9 versus Cd274, Ccrl2 and Slamf7
# -----------------------------------------------------------------------------

correlation_data <- FetchData(
  figure_object,
  vars = c(
    "Cxcl9",
    CORRELATION_TARGETS,
    "seurat_clusters",
    "timepoint"
  ),
  layer = "data"
) |>
  filter(
    timepoint == "Resolution",
    as.character(seurat_clusters) == RESOLUTION_CLUSTER_OF_INTEREST
  )

if (nrow(correlation_data) < 3) {
  stop("Fewer than three cells were found in resolution cluster 5.")
}

safe_spearman <- function(x, y) {
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]

  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(data.frame(rho = NA_real_, p_value = NA_real_, n_cells = length(x)))
  }

  test <- suppressWarnings(
    stats::cor.test(x, y, method = "spearman", exact = FALSE)
  )
  data.frame(
    rho = unname(test$estimate),
    p_value = test$p.value,
    n_cells = length(x)
  )
}

correlation_results <- bind_rows(lapply(CORRELATION_TARGETS, function(gene) {
  cbind(
    data.frame(gene_x = "Cxcl9", gene_y = gene),
    safe_spearman(correlation_data$Cxcl9, correlation_data[[gene]])
  )
})) |>
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    timepoint = "Resolution",
    cluster = RESOLUTION_CLUSTER_OF_INTEREST
  )

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 2.2e-16,
      "<2.2e-16",
      ifelse(p < 0.001, formatC(p, format = "e", digits = 2), signif(p, 3))
    )
  )
}

correlation_plot_data <- correlation_results |>
  mutate(
    gene_x = factor(gene_x, levels = "Cxcl9"),
    gene_y = factor(gene_y, levels = rev(CORRELATION_TARGETS)),
    label = paste0(
      "rho = ", ifelse(is.na(rho), "NA", sprintf("%.2f", rho)),
      "\nraw p = ", format_p(p_value)
    )
  )

p_correlation <- ggplot(
  correlation_plot_data,
  aes(gene_x, gene_y, fill = rho)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    limits = c(-0.5, 0.5),
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    oob = scales::squish,
    na.value = "grey90",
    name = expression(rho)
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "italic"),
    axis.text.y = element_text(face = "italic")
  )

write.csv(
  correlation_results,
  file.path(OUTPUT_DIR, "FigS2H_Cxcl9_Spearman_correlations.csv"),
  row.names = FALSE
)
save_panel(
  "FigS2H_Cxcl9_Spearman_heatmap",
  p_correlation,
  width = 3.2,
  height = 3.0
)

# Individual panels are the intended reusable outputs. This combined preview
# mirrors the left-to-right order of Fig. S2H without recreating page layout.
p_combined <- p_umap + p_marker_dot + p_iam_dot + p_correlation +
  patchwork::plot_layout(widths = c(2.0, 1.7, 1.0, 0.9)) +
  patchwork::plot_annotation(title = "Hegarty et al. replication dataset")
save_panel("FigS2H_combined_preview", p_combined, width = 15, height = 3.6)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)
