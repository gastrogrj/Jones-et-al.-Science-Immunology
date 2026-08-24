# ACOD1 gut microbiome analysis
#
# This script imports PacBio abundance and taxonomy data, constructs a
# phyloseq object, calculates alpha diversity, and performs a Bray-Curtis
# principal coordinates analysis (PCoA) with PERMANOVA.

# Package setup -----------------------------------------------------------

# Install missing packages before running this script. Bioconductor packages
# (e.g., phyloseq) should be installed with BiocManager::install().
required_packages <- c(
  "cowplot", "dplyr", "ggExtra", "ggplot2", "microbiome",
  "pairwiseAdonis", "phyloseq", "rlang", "scales", "stringr", "tidyr",
  "vegan"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

# Project paths -----------------------------------------------------------

# Change this path to the project folder.
project_dir <- path.expand(
  ""
)

abundance_file <- file.path(
  project_dir,
  "acod1_output/results/best_tax_merged_freq_tax.tsv"
)
metadata_file <- file.path(project_dir, "input_files/acod1_metadata.tsv")
output_dir <- file.path(project_dir, "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Specify group colours
group_colors <- c(
  KO = "#FF2600",
  WT_HET = "#919191"
)

# Functions ---------------------------------------------------------------

# Plot Bray-Curtis beta diversity and test group differences.
#
# Args:
#   physeq_object: A phyloseq object containing counts and sample metadata.
#   comparator: Name of the sample metadata column used to define groups.
#   colors: Named character vector mapping group names to hex colours.
#
# Returns:
#   A ggExtra marginal plot. PERMANOVA and pairwise PERMANOVA results are
#   printed to the console.
plot_beta_diversity <- function(
    physeq_object,
    comparator,
    colors = group_colors
) {
  metadata <- data.frame(
    phyloseq::sample_data(physeq_object),
    check.names = FALSE
  )

  if (!comparator %in% names(metadata)) {
    stop("Comparator column not found in sample metadata: ", comparator)
  }

  observed_groups <- unique(as.character(metadata[[comparator]]))
  missing_colors <- setdiff(observed_groups, names(colors))
  if (length(missing_colors) > 0L) {
    stop(
      "No colour supplied for: ",
      paste(missing_colors, collapse = ", ")
    )
  }

  # Remove taxa with zero abundance across all retained samples.
  physeq_object <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(physeq_object) > 0,
    physeq_object
  )

  # Calculate Bray-Curtis distances and PCoA coordinates.
  bray_distance <- phyloseq::distance(physeq_object, method = "bray")
  bray_pcoa <- phyloseq::ordinate(
    physeq_object,
    method = "PCoA",
    distance = bray_distance
  )

  variance_axis1 <- bray_pcoa$values$Relative_eig[1] * 100
  variance_axis2 <- bray_pcoa$values$Relative_eig[2] * 100
  axis1_label <- sprintf("PCo1 (%.1f%%)", variance_axis1)
  axis2_label <- sprintf("PCo2 (%.1f%%)", variance_axis2)

  # Test whether community composition differs between comparator groups.
  permanova_formula <- stats::reformulate(comparator, response = "bray_distance")
  permanova_result <- vegan::adonis2(
    permanova_formula,
    data = metadata
  )
  print(permanova_result)

  # Run post-hoc pairwise comparisons between comparator groups.
  pairwise_result <- pairwiseAdonis::pairwise.adonis(
    bray_distance,
    metadata[[comparator]]
  )
  print(pairwise_result)

  p_value <- permanova_result$`Pr(>F)`[1]
  title_text <- sprintf("PERMANOVA p = %.3f", p_value)

  # Generate the PCoA plot, including group confidence ellipses.
  ordination_plot <- phyloseq::plot_ordination(
    physeq_object,
    bray_pcoa,
    color = comparator
  )

  beta_plot <- ordination_plot +
    ggplot2::geom_point(size = 5) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = !!rlang::sym(comparator)),
      linetype = 1
    ) +
    ggplot2::scale_color_manual(values = colors) +
    cowplot::theme_cowplot(font_size = 20) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = axis1_label,
      y = axis2_label,
      title = title_text
    ) +
    ggplot2::theme(
      legend.position = c(0.99, 0.99),
      legend.justification = c(1, 1),
      legend.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black",
        linewidth = 1,
        fill = NA
      ),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7),
        color = "black",
        linewidth = 0.5
      )
    )

  # Add marginal boxplots showing the distributions along each PCoA axis.
  ggExtra::ggMarginal(
    beta_plot,
    type = "boxplot",
    groupFill = TRUE,
    size = 10
  )
}

# Import data -------------------------------------------------------------

pacbio <- utils::read.delim(
  abundance_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "#"
)

acod1_metadata <- utils::read.delim(
  metadata_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "#"
)

# Construct the taxonomy table by splitting semicolon-delimited ranks and
# removing database prefixes such as "k__" and "p__".
taxonomy <- pacbio |>
  dplyr::select(Taxon) |>
  tidyr::separate(
    Taxon,
    into = c(
      "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"
    ),
    sep = ";",
    fill = "right"
  ) |>
  dplyr::mutate(
    dplyr::across(
      Domain:Species,
      ~ stringr::str_trim(stringr::str_remove(.x, "^[a-z]__"))
    )
  ) |>
  as.matrix()

taxonomy_table <- phyloseq::tax_table(taxonomy)

# The remaining sample columns form the ASV abundance matrix.
asv_counts <- pacbio |>
  dplyr::select(-Sequence, -Taxon, -Confidence) |>
  as.matrix()

asv_table <- phyloseq::otu_table(asv_counts, taxa_are_rows = TRUE)
sample_metadata <- phyloseq::sample_data(acod1_metadata)

physeq <- phyloseq::phyloseq(sample_metadata, asv_table, taxonomy_table)

# Relative-abundance filtering -------------------------------------------

# Convert counts to within-sample relative abundances, retain taxa whose
# summed relative abundance across samples exceeds 0.001, and renormalise
# each sample after filtering.
physeq_relative <- phyloseq::transform_sample_counts(
  physeq,
  function(x) x / sum(x)
)

physeq_relative <- phyloseq::filter_taxa(
  physeq_relative,
  function(x) sum(x) > 0.001,
  prune = TRUE
)

physeq_relative <- phyloseq::transform_sample_counts(
  physeq_relative,
  function(x) x / sum(x)
)

# Alpha diversity ---------------------------------------------------------

alpha_diversity <- microbiome::alpha(physeq_relative, index = "all")
alpha_diversity_with_metadata <- merge(
  acod1_metadata,
  alpha_diversity,
  by = "row.names"
)

utils::write.csv(
  alpha_diversity_with_metadata,
  file.path(output_dir, "alpha_diversity_with_metadata.csv"),
  row.names = FALSE
)

# Beta diversity ----------------------------------------------------------

condition_beta <- plot_beta_diversity(physeq_relative, "genotype_2")

ggplot2::ggsave(
  filename = file.path(output_dir, "beta_WTHET_KO.svg"),
  plot = condition_beta,
  width = 7.2,
  height = 7,
  units = "in",
  device = "svg"
)

