# Jones et al. - Figure 5 (human scRNA-seq panels only)
#
# Produces the analysis underlying:
#   Fig. 5A (human): Reactome GSEA, IFN-signalling enrichment plot and
#                    leading-edge expression heatmap
#   Fig. 5C (human): selected transcription-factor activities
#   Fig. 5D (human): PROGENy pathway activities
#
# Mouse panels and non-scRNA-seq panels are intentionally excluded.
# Run this script from the repository root. It uses relative paths only.

suppressPackageStartupMessages({
  library(decoupleR)
  library(dplyr)
  library(fgsea)
  library(fmsb)
  library(ggplot2)
  library(pheatmap)
  library(progeny)
  library(readr)
  library(tibble)
  library(tidyr)
})

set.seed(1)

# -----------------------------------------------------------------------------
# Inputs and outputs
# -----------------------------------------------------------------------------

# Gene-by-cell-type specificity/effect-size matrix used in the original script.
# Required columns before renaming:
#   gene, Macrophage, Macrophage.1, Macrophage.2,
#   Mac_intermediate_.2., Monocyte
INPUT_EXPRESSION <- file.path(
  "data", "human_scRNAseq", "fr003_myeloid.celltype_label.esmu.csv.gz"
)

# MSigDB Reactome collection used for the submitted analysis.
REACTOME_GMT <- file.path(
  "data", "gmt", "c2.cp.reactome.v2023.2.Hs.symbols.gmt"
)

OUTPUT_DIR <- file.path("output", "Figure_5_human")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(INPUT_EXPRESSION), file.exists(REACTOME_GMT))

CELL_TYPES <- c(
  "CD163- macrophage",
  "CD163+ macrophage",
  "ITGAX+ macrophage",
  "CXCL9/10+ monocyte",
  "Monocyte"
)

IAM_CELL_TYPE <- "CXCL9/10+ monocyte"
CELL_TYPES_FIGURE <- c(IAM_CELL_TYPE, setdiff(CELL_TYPES, IAM_CELL_TYPE))

CELL_TYPE_COLOURS <- c(
  "CD163- macrophage" = "#AD57DE",
  "CD163+ macrophage" = "#40B043",
  "ITGAX+ macrophage" = "#FF527D",
  "CXCL9/10+ monocyte" = "#5AD1FF",
  "Monocyte" = "#96B200"
)

IG_GENES <- c("IGHA1", "JCHAIN", "IGLC2", "IGLC3", "IGHM", "IGKC")

# -----------------------------------------------------------------------------
# Read and validate the human myeloid expression/specificity matrix
# -----------------------------------------------------------------------------

expression_df <- readr::read_csv(INPUT_EXPRESSION, show_col_types = FALSE)

required_original_columns <- c(
  "gene",
  "Macrophage",
  "Macrophage.1",
  "Macrophage.2",
  "Mac_intermediate_.2.",
  "Monocyte"
)

missing_columns <- setdiff(required_original_columns, colnames(expression_df))
if (length(missing_columns) > 0) {
  stop(
    "The human myeloid matrix is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

expression_df <- expression_df |>
  select(all_of(required_original_columns)) |>
  rename(
    `CD163- macrophage` = Macrophage,
    `CD163+ macrophage` = Macrophage.1,
    `ITGAX+ macrophage` = Macrophage.2,
    `CXCL9/10+ monocyte` = Mac_intermediate_.2.
  ) |>
  filter(!is.na(gene), gene != "") |>
  group_by(gene) |>
  summarise(across(all_of(CELL_TYPES), mean, na.rm = TRUE), .groups = "drop")

expr_mat <- expression_df |>
  column_to_rownames("gene") |>
  as.matrix()

storage.mode(expr_mat) <- "numeric"
stopifnot(identical(colnames(expr_mat), CELL_TYPES))

# -----------------------------------------------------------------------------
# Figure 5A (human): Reactome FGSEA
# -----------------------------------------------------------------------------

# Rank all measured genes by the CXCL9/10+ monocyte score. The historical
# script truncated this list in one of several competing FGSEA blocks; the
# final analysis uses one complete, uniquely named ranking.
human_ranks <- expr_mat[, IAM_CELL_TYPE]
human_ranks <- human_ranks[
  is.finite(human_ranks) &
    human_ranks != 0 &
    !(names(human_ranks) %in% IG_GENES)
]
human_ranks <- sort(human_ranks, decreasing = TRUE)

reactome_pathways <- fgsea::gmtPathways(REACTOME_GMT)

fgsea_human <- fgsea::fgsea(
  pathways = reactome_pathways,
  stats = human_ranks,
  minSize = 50,
  maxSize = 500,
  scoreType = "pos",
  eps = 0,
  nPermSimple = 10000
) |>
  as_tibble() |>
  arrange(padj, desc(NES))

fgsea_export <- fgsea_human |>
  mutate(leadingEdge = vapply(
    leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  ))

write_csv(
  fgsea_export,
  file.path(OUTPUT_DIR, "Fig5A_human_Reactome_FGSEA_all.csv")
)

# The six pathways displayed in the final figure. Only pathways meeting the
# legend-specified FDR < 0.05 criterion are plotted.
display_pathways <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERFERON_SIGNALING",
  "REACTOME_SIGNALING_BY_INTERLEUKINS",
  "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
  "REACTOME_ADAPTIVE_IMMUNE_SYSTEM",
  "REACTOME_NEUTROPHIL_DEGRANULATION"
)

fgsea_bar_data <- fgsea_human |>
  filter(padj < 0.05, pathway %in% display_pathways) |>
  mutate(
    pathway_label = tools::toTitleCase(
      gsub("_", " ", sub("^REACTOME_", "", pathway))
    ),
    pathway_label = factor(
      pathway_label,
      levels = rev(pathway_label[order(NES)])
    )
  )

if (nrow(fgsea_bar_data) == 0) {
  stop("None of the displayed Figure 5A pathways met FDR < 0.05.")
}

p_fig5a_bar <- ggplot(
  fgsea_bar_data,
  aes(x = pathway_label, y = NES)
) +
  geom_col(fill = "#E97878", width = 0.72) +
  coord_flip() +
  labs(x = NULL, y = "Normalised enrichment score") +
  theme_classic(base_size = 10)

ggsave(
  file.path(OUTPUT_DIR, "Fig5A_human_Reactome_barplot.pdf"),
  p_fig5a_bar,
  width = 5.2,
  height = 3.6
)

if (!"REACTOME_INTERFERON_SIGNALING" %in% names(reactome_pathways)) {
  stop("REACTOME_INTERFERON_SIGNALING was not found in the supplied GMT file.")
}

p_fig5a_enrichment <- fgsea::plotEnrichment(
  reactome_pathways[["REACTOME_INTERFERON_SIGNALING"]],
  human_ranks
) +
  labs(title = "IFN signalling", x = "Gene rank", y = "Enrichment score") +
  theme_classic(base_size = 10)

ggsave(
  file.path(OUTPUT_DIR, "Fig5A_human_IFN_signalling_enrichment.pdf"),
  p_fig5a_enrichment,
  width = 5.2,
  height = 3.6
)

# Genes displayed beside the enrichment curve in the final Figure 5A.
leading_edge_genes <- c(
  "IFITM3", "ISG15", "SOCS3", "JAK1", "HSPA1A",
  "IFITM2", "GBP5", "IFI6", "HSPA1B", "ISG20"
)

missing_leading_edge <- setdiff(leading_edge_genes, rownames(expr_mat))
if (length(missing_leading_edge) > 0) {
  stop(
    "Figure 5A leading-edge genes are absent from the input matrix: ",
    paste(missing_leading_edge, collapse = ", ")
  )
}

leading_edge_df <- tibble(
  gene = factor(leading_edge_genes, levels = rev(leading_edge_genes)),
  score = expr_mat[leading_edge_genes, IAM_CELL_TYPE]
)

p_fig5a_leading_edge <- ggplot(
  leading_edge_df,
  aes(x = "Human", y = gene, fill = score)
) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", direction = 1) +
  labs(x = NULL, y = NULL, fill = "Average\nlogFC") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic"))

ggsave(
  file.path(OUTPUT_DIR, "Fig5A_human_IFN_leading_edge_heatmap.pdf"),
  p_fig5a_leading_edge,
  width = 2.5,
  height = 4.0
)

write_csv(
  leading_edge_df |> mutate(gene = as.character(gene)),
  file.path(OUTPUT_DIR, "Fig5A_human_IFN_leading_edge_values.csv")
)

# -----------------------------------------------------------------------------
# Figure 5C (human): transcription-factor activity with decoupleR/CollecTRI
# -----------------------------------------------------------------------------

collectri <- decoupleR::get_collectri(
  organism = "human",
  split_complexes = FALSE
)

ulm_scores <- decoupleR::run_ulm(
  mat = expr_mat,
  network = collectri,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

ulm_mat <- ulm_scores |>
  filter(statistic == "ulm") |>
  select(source, condition, score) |>
  pivot_wider(names_from = condition, values_from = score) |>
  column_to_rownames("source") |>
  as.matrix()

selected_tfs <- c("STAT3", "IRF5", "IRF4", "IRF2", "IRF9", "IRF8")
missing_tfs <- setdiff(selected_tfs, rownames(ulm_mat))
if (length(missing_tfs) > 0) {
  stop("Selected TFs were not inferred: ", paste(missing_tfs, collapse = ", "))
}

tf_subset <- ulm_mat[selected_tfs, CELL_TYPES_FIGURE, drop = FALSE]

write_csv(
  as.data.frame(tf_subset) |> rownames_to_column("TF"),
  file.path(OUTPUT_DIR, "Fig5C_human_TF_activity.csv")
)

tf_breaks <- seq(-5, 8, length.out = 101)
tf_colours <- colorRampPalette(c("blue", "white", "red"))(100)

pdf(
  file.path(OUTPUT_DIR, "Fig5C_human_TF_activity_heatmap.pdf"),
  width = 6.0,
  height = 4.2
)
pheatmap::pheatmap(
  tf_subset,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = 1,
  color = tf_colours,
  breaks = tf_breaks,
  border_color = "white",
  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 45
)
dev.off()

# -----------------------------------------------------------------------------
# Figure 5D (human): PROGENy pathway activity
# -----------------------------------------------------------------------------

# PROGENy expects genes in rows and samples/cell types in columns; its output
# has cell types in rows and pathways in columns.
progeny_scores <- progeny::progeny(
  expr = expr_mat,
  scale = TRUE,
  organism = "Human",
  top = 500,
  perm = 1,
  verbose = FALSE
)

progeny_scores <- progeny_scores[CELL_TYPES, , drop = FALSE]

write_csv(
  as.data.frame(progeny_scores) |> rownames_to_column("cell_type"),
  file.path(OUTPUT_DIR, "Fig5D_human_PROGENy_activity.csv")
)

# Radar axes are ordered by decreasing activity in the IAM population, as in
# the original analysis.
pathway_order <- names(sort(progeny_scores[IAM_CELL_TYPE, ], decreasing = TRUE))
progeny_ordered <- progeny_scores[, pathway_order, drop = FALSE]

radar_data <- as.data.frame(rbind(
  apply(progeny_ordered, 2, max),
  apply(progeny_ordered, 2, min),
  progeny_ordered
))

pdf(
  file.path(OUTPUT_DIR, "Fig5D_human_PROGENy_radar.pdf"),
  width = 6.3,
  height = 6.3
)
fmsb::radarchart(
  radar_data,
  axistype = 0,
  pcol = unname(CELL_TYPE_COLOURS[CELL_TYPES]),
  pfcol = scales::alpha(unname(CELL_TYPE_COLOURS[CELL_TYPES]), 0.12),
  plwd = 2,
  plty = 1,
  cglcol = "grey70",
  cglty = 1,
  cglwd = 0.8,
  axislabcol = NA,
  vlcex = 0.75
)
legend(
  "bottom",
  legend = CELL_TYPES,
  col = unname(CELL_TYPE_COLOURS[CELL_TYPES]),
  lty = 1,
  lwd = 2,
  bty = "n",
  cex = 0.75,
  horiz = TRUE
)
dev.off()

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)
