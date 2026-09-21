# Shared helpers for the mouse figure scripts. Source from the repository root.
# RNA plots use ln(1 + counts / total counts * 10000); values are dimensionless.
# Cell-level tests describe cells from pooled libraries, not independent mice.

check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse=", "))
}
check_packages(c("Seurat", "SeuratObject", "ggplot2", "patchwork", "svglite", "readr"))

check_features <- function(object, genes) {
  missing <- setdiff(genes, rownames(object[["RNA"]]))
  if (length(missing)) stop("Missing RNA genes: ", paste(missing, collapse=", "))
}

join_rna <- function(object, normalise=FALSE) {
  if (!inherits(object, "Seurat")) stop("Input must be a Seurat object.")
  if (!"RNA" %in% names(object@assays)) stop("RNA assay missing.")
  if (inherits(object[["RNA"]], "Assay5"))
    object[["RNA"]] <- SeuratObject::JoinLayers(object[["RNA"]])
  SeuratObject::DefaultAssay(object) <- "RNA"
  if (normalise) object <- Seurat::NormalizeData(object, assay="RNA",
      normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
  object
}

COLON_MAP <- readr::read_csv("data/annotations/colon_cluster_map.csv", show_col_types=FALSE,
                           col_types=readr::cols(.default=readr::col_character()))
if (anyDuplicated(COLON_MAP$seurat_cluster) || anyDuplicated(COLON_MAP$figure_cluster))
  stop("Cluster mapping must be one-to-one.")
FIGURE_CLUSTER_COLOURS <- setNames(COLON_MAP$colour, COLON_MAP$figure_cluster)
label_colon <- function(object) {
  if (!all(c("seurat_clusters","timepoint") %in% colnames(object[[]])))
    stop("Colon object needs seurat_clusters and timepoint metadata.")
  raw <- as.character(object$seurat_clusters)
  mapped <- COLON_MAP$figure_cluster[match(raw,COLON_MAP$seurat_cluster)]
  if (anyNA(mapped)) stop("Unmapped Seurat cluster: check colon_cluster_map.csv.")
  if ("figure_cluster" %in% colnames(object[[]]) &&
      any(as.character(object$figure_cluster) != mapped, na.rm=TRUE))
    stop("Existing figure_cluster disagrees with colon_cluster_map.csv. Resolve before plotting.")
  object$figure_cluster <- factor(mapped, levels=as.character(1:5))
  SeuratObject::Idents(object) <- "figure_cluster"
  object
}

save_figure <- function(stem, plot, width_mm, height_mm) {
  dir.create(dirname(stem),recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(paste0(stem,".pdf"),plot,width=width_mm,height=height_mm,units="mm",
    device=if(capabilities("cairo")) grDevices::cairo_pdf else "pdf",bg="white")
  ggplot2::ggsave(paste0(stem,".svg"),plot,width=width_mm,height=height_mm,units="mm",
    device=svglite::svglite,bg="white")
}

feature_grid <- function(object, genes, ncol=2, limits=NULL) {
  object <- join_rna(object)
  check_features(object,genes)
  if(is.null(limits)) {
    v <- Seurat::FetchData(object,vars=genes,layer="data")
    limits <- c(0,max(1,ceiling(max(as.matrix(v),na.rm=TRUE))))
  }
  pp <- Seurat::FeaturePlot(object, features=genes,ncol=ncol,order=TRUE,
       raster=TRUE,combine=FALSE,min.cutoff=limits[1],max.cutoff=limits[2])
  pp <- lapply(pp,function(p) p + ggplot2::scale_colour_gradientn(
    colours=c("grey90","#54278F"),limits=limits,oob=scales::squish,
    name="RNA expression\nlog1p(counts per 10,000)") + ggplot2::theme_void(base_size=7) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=7,face="italic"),
                   legend.title=ggplot2::element_text(size=7),legend.text=ggplot2::element_text(size=7)))
  patchwork::wrap_plots(pp,ncol=ncol,guides="collect") &
    ggplot2::theme(legend.position="bottom",legend.direction="horizontal")
}

safe_spearman <- function(x,y) {
  ok <- is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok]
  if(length(x)<3 || stats::sd(x)==0 || stats::sd(y)==0)
    return(data.frame(rho=NA_real_,p_value=NA_real_,n=length(x)))
  a<-suppressWarnings(stats::cor.test(x,y,method="spearman",exact=FALSE,alternative="two.sided"))
  data.frame(rho=unname(a$estimate),p_value=a$p.value,n=length(x))
}
scale_or_zero <- function(x) if(length(x)<2 || !is.finite(stats::sd(x)) || stats::sd(x)==0)
  rep(0,length(x)) else as.numeric(scale(x))

de_table <- function(object,cells1,cells2,min_pct=0,only_positive=FALSE) {
  if(length(cells1)<3 || length(cells2)<3 || length(intersect(cells1,cells2)))
    stop("DE requires two disjoint groups with at least three cells each.")
  a<-Seurat::FindMarkers(object,ident.1=cells1,ident.2=cells2,assay="RNA",
       test.use="wilcox",min.pct=min_pct,logfc.threshold=0,only.pos=only_positive)
  if(!"avg_log2FC" %in% names(a)) stop("Seurat v5 avg_log2FC output is required.")
  a$gene<-rownames(a);a$log2FC<-a$avg_log2FC
  a
}

cell_groups <- function(object,figures,condition=NULL) {
  m<-object[[]]
  keep<-as.character(m$figure_cluster)%in%figures
  if(!is.null(condition)) keep<-keep & as.character(m$timepoint)%in%condition
  rownames(m)[keep]
}

read_gene_expression <- function(path) {
  x<-Seurat::Read10X(path)
  if(is.list(x)) {
    if(!"Gene Expression"%in%names(x)) stop("No Gene Expression matrix: ",path)
    x<-x[["Gene Expression"]]
  };x
}

export_cell_manifest <- function(object,path) {
  m<-object[[]];m$cell<-rownames(m)
  readr::write_csv(m,path)
}
export_session <- function(outdir) {
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
}

# Explicit rectangle fills prevent Affinity from introducing black backgrounds.
acod1_fix_svg <- function(path) {
  doc <- xml2::read_xml(path)
  rects <- xml2::xml_find_all(doc, ".//*[local-name()='rect']")
  for (r in rects) {
    style <- xml2::xml_attr(r, "style")
    has_style_fill <- !is.na(style) && grepl("(^|;)\\s*fill\\s*:", style)
    if (is.na(xml2::xml_attr(r, "fill")) && !has_style_fill) {
      xml2::xml_set_attr(r, "fill", "none")
    }
  }
  xml2::write_xml(doc, path)
  invisible(path)
}

# scores: numeric matrix with cell populations in ROWS and 14 pathways in COLUMNS.
# These should be existing PROGENy scores, NOT the fmsb max/min padding rows.
# This retains the uploaded radar's PER-PATHWAY min/max scaling and centre gap.
# Therefore the rings represent percent of each pathway's observed score range,
# NOT a common z-score, percent of cells, or percentage biological activation.
export_progeny_radar_7pt <- function(
    scores, stem, colours, pathway_order = colnames(scores),
    width_mm = 58, height_mm = 62, font_family = "Arial") {
  check_packages(c("fmsb", "svglite", "xml2"))
  scores <- as.matrix(scores)
  if (!is.numeric(scores) || nrow(scores) < 2L || ncol(scores) != 14L ||
      any(!is.finite(scores)) || is.null(rownames(scores)) ||
      is.null(colnames(scores)) || anyDuplicated(rownames(scores)) ||
      anyDuplicated(colnames(scores))) {
    stop("Supply finite scores: at least two named populations x 14 named pathways.")
  }
  if (!setequal(pathway_order, colnames(scores)) || length(pathway_order) != 14L) {
    stop("pathway_order must contain each of the 14 score columns exactly once.")
  }
  scores <- scores[, pathway_order, drop = FALSE]
  if (is.null(names(colours)) || !all(rownames(scores) %in% names(colours))) {
    stop("Use a named colour vector with a colour for every row of scores.")
  }
  colours <- colours[rownames(scores)]
  hi <- apply(scores, 2, max)
  lo <- apply(scores, 2, min)
  if (any(hi == lo)) {
    stop("No between-population range for: ",
         paste(names(hi)[hi == lo], collapse = ", "),
         ". A min/max radar is undefined for these pathways.")
  }
  radar_data <- as.data.frame(rbind(max = hi, min = lo, scores), check.names = FALSE)
  relative <- 100 * sweep(sweep(scores, 2, lo, "-"), 2, hi - lo, "/")
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(scores, paste0(stem, "_scores.csv"))
  utils::write.csv(relative, paste0(stem, "_within_pathway_range_percent.csv"))
  utils::write.csv(data.frame(pathway = names(lo), minimum = lo, maximum = hi),
                   paste0(stem, "_pathway_bounds.csv"), row.names = FALSE)

  draw <- function() {
    # vlcex is a MULTIPLIER, not a point size. 7 pt device + cex 1 + vlcex 1
    # gives 7 pt pathway labels. No mfrow layout or later scaling is used.
    graphics::par(ps = 7, cex = 1, family = "sans",
                  mar = c(2.5, 2.0, 1.2, 2.0), xpd = NA)
    fmsb::radarchart(
      radar_data, axistype = 1, seg = 4, maxmin = TRUE, centerzero = FALSE,
      vlabels = colnames(scores), vlcex = 1, calcex = 1,
      caxislabels = c("0", "25", "50", "75", "100"),
      pcol = unname(colours), pfcol = grDevices::adjustcolor(colours, alpha.f = 0.15),
      plwd = 1, plty = 1, pty = 32,
      cglcol = "grey80", cglty = 1, cglwd = 0.5, axislabcol = "grey30"
    )
    graphics::mtext("Within-pathway range (%)", side = 1, line = 1.2, cex = 1)
  }
  save_device <- function(file, kind, draw_fun, w = width_mm, h = height_mm) {
    if (kind == "svg") {
      svglite::svglite(file, width = w / 25.4, height = h / 25.4,
                      pointsize = 7, bg = "white", fix_text_size = FALSE,
                      system_fonts = list(sans = font_family))
    } else if (capabilities("cairo")) {
      grDevices::cairo_pdf(file, width = w / 25.4, height = h / 25.4,
                          pointsize = 7, family = font_family, bg = "white")
    } else {
      # Standard PDF fonts may substitute Helvetica for Arial.
      grDevices::pdf(file, width = w / 25.4, height = h / 25.4,
                     pointsize = 7, family = "Helvetica", bg = "white")
    }
    device_id <- grDevices::dev.cur()
    on.exit(grDevices::dev.off(device_id), add = TRUE)
    draw_fun()
  }
  save_device(paste0(stem, ".svg"), "svg", draw)
  acod1_fix_svg(paste0(stem, ".svg"))
  save_device(paste0(stem, ".pdf"), "pdf", draw)

  # Separate key: retain the existing manuscript key if its colours match.
  draw_key <- function() {
    grid::grid.newpage()
    yy <- seq(0.86, 0.14, length.out = nrow(scores))
    grid::grid.points(x = rep(0.04, nrow(scores)), y = yy, pch = 16,
                      size = grid::unit(1.7, "mm"), gp = grid::gpar(col = colours))
    grid::grid.text(rownames(scores), x = 0.09, y = yy, just = "left",
                    gp = grid::gpar(fontsize = 7, fontfamily = "sans"))
  }
  key_height <- max(18, nrow(scores) * 4)
  save_device(paste0(stem, "_key.svg"), "svg", draw_key, 58, key_height)
  acod1_fix_svg(paste0(stem, "_key.svg"))
  save_device(paste0(stem, "_key.pdf"), "pdf", draw_key, 58, key_height)
  invisible(list(scores = scores, plotted_percent = relative, bounds = rbind(lo, hi)))
}
