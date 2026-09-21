# Jones et al. Science Immunology analysis code

Code accompanying **Concurrent induction of proinflammatory and regulatory programmes in inflammation-associated monocyte-derived cells in inflammatory bowel disease**.

This revision covers the mouse single-cell analyses in Figures 2, 3A, 4, 5B/C/D/H and the associated scRNA-seq panels in Figures S2, S3 and S4. It was compared with the supplied v22 manuscript and current figures on 21 September 2026. The starting GitHub commit was `a2f2080bb3e708296078b0c0dee920d725d2c9dc`.

**Validation status:** source/figure review and syntax checks have been completed. The analyses have not been executed end to end against the experimental data. GSE345123 was private at the time of review. Some original figure-generating inputs and settings remain unavailable. See [CODE_REVIEW.md](CODE_REVIEW.md) before declaring this a reproducible manuscript release.


## Panel coverage

| Panel | Script | Scope and provenance |
| --- | --- | --- |
| 2A–F; S2C/D/F/I | `Figure_2_and_S2_revised.R` | Colon processing, clustering, markers, proportions, coexpression, cluster stability and Cd274 groups |
| 3A; S3A/B | `Figure_3_mouse_revised.R` | Added explicit DE and Acod1 plots; original 3A settings and S3 selections require confirmation |
| 4A/B/D–F; S4B/E/F | `Figure_4_and_S4_revised.R` | Blood/colon analysis, fixed original signatures, UCell comparison and Monocle3 trajectories |
| Mouse 5B/C/D/H | `Figure_5_mouse_revised.R` | GSEA, CollecTRI ULM, PROGENy radar and Il1b–Acod1 correlation; 5B is reconstructed from the legend |
| S2H | `Figure_S2H_revised.R` | Independent Hegarty et al. recovery dataset; removed its duplicate, incomplete block from the Figure 2 script |
| S2E | `Export_Figure_S2E_input.R`, `Figure_S2E_scanpy.py`, `Figure_S2E_concordance.R` | Export final colon cells, rerun recovered Scanpy workflow, calculate actual marker overlap and correlations |
| Shared plotting and validation | `mouse_analysis_helpers.R` | RNA units, cluster mapping, PDF/SVG export and 7 pt radar labels |

## Deposited data

| Dataset | Accession or DOI | Access checked for this revision |
| --- | --- | --- |
| Mouse blood and colon scRNA-seq | [GSE345123](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345123) | Private; NCBI displayed release date 27 August 2027 on 21 September 2026 |
| Mouse monocyte NanoString | [GSE345013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345013) | Private; NCBI displayed release date 27 August 2027 |
| Acod1 microbiota 16S | [PRJNA1518240](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1518240) | Accession supplied in manuscript; contents not audited here |
| Mouse Xenium spatial data | [10.5281/zenodo.22829361](https://doi.org/10.5281/zenodo.22829361) | Public record, restricted files |
| Human intestinal scRNA-seq | [10.5281/zenodo.8301000](https://doi.org/10.5281/zenodo.8301000) | Manuscript citation; human data not analysed in this revision |
| Independent recovery scRNA-seq for S2H | Hegarty et al., Mucosal Immunology 19, 1624–1635 (2026) |
| Cadinu et al. MERFISH | Cited separately in the manuscript | Add the exact input file/version and original analysis script |

## Sampling and inference

Five female mice per condition contributed paired blood and colon samples. Cells were pooled by tissue and condition before sequencing, producing four libraries: naive colon, DSS colon, naive blood and DSS blood. The 20,000 cells are input cells per pooled sample, not retained cells after QC.

There is **one pooled library per tissue–condition combination**, not five independent scRNA-seq libraries. Cell-level DE/correlation P values and proportions do not establish mouse-level replication. Cluster-condition expression averages are descriptive profiles, not independent biological pseudobulk replicates. The code does not create mouse identifiers or replicate-level tests that the experiment cannot support.

## Input layout

Run scripts with the repository root as the working directory. The `data/annotations/colon_cluster_map.csv` file is included and must accompany the scripts.

Raw 10x directories, if reprocessing is needed:

| Directory under `data/mouse_scRNAseq/` | Tissue | Condition |
| --- | --- | --- |
| `G1/` | Colon | Naive |
| `G2/` | Colon | DSS day 6 |
| `B1/` | Blood | Naive |
| `B2/` | Blood | DSS day 6 |

Each directory needs both `raw_feature_bc_matrix/` and `filtered_feature_bc_matrix/`, including matrix, barcode and feature files. SoupX needs raw droplet counts as well as the filtered matrix. The public GEO file names/GSM identities could not be checked while the record was private; add that mapping after confirming the deposit.

For reproduction of the final figures, use the **original saved objects**:

| Path | Required contents |
| --- | --- |
| `data/objects/colon_final.rds` | Final Figure 2 Csf1r+ Seurat object, original UMAP, RNA counts/data, `seurat_clusters`, `timepoint`, `sample`; SCT assay and original `SCT_snn` neighbour graph needed for downstream stability work |
| `data/objects/colon_full_postQC.rds` | Full post-QC colon object before Csf1r+ selection, for S2C |
| `data/objects/blood_colon_final.rds` | Original `Combined_Preprocess_7_10` object for Figure 4; RNA/SCT assays, UMAP, `tissue`, `timepoint`, `seurat_clusters` |
| `data/objects/blood_colon_before_cluster10.rds` | Original `Combined_Preprocess_7` object used for S4B; includes the subsequently removed cluster 10 |
| `data/objects/recovery_csf1r.rds` | Final independent recovery object for S2H, with UMAP, RNA assay, `timepoint` (`Naive`/`Resolution`) and clusters 0–6 |
| `data/objects/blood_only_S3A.rds` | Original blood-only embedding for S3A; its upstream script is still required for complete provenance |
| `data/objects/Monocle3_trajectory_objects.rds` | Named list of final CDS objects: `naive_blood_colon`, `dss_blood_colon`, `naive_colon_only`, `dss_colon_only` |

Save these objects from the R sessions in which they were generated. For example:

```r
dir.create("data/objects", recursive = TRUE, showWarnings = FALSE)
# In the original COLON analysis session:
saveRDS(Combined_Preprocess, "data/objects/colon_final.rds", compress = "xz")
# In the original BLOOD/COLON analysis session:
saveRDS(Combined_Preprocess_7_10, "data/objects/blood_colon_final.rds", compress = "xz")
saveRDS(Combined_Preprocess_7, "data/objects/blood_colon_before_cluster10.rds", compress = "xz")
```

The name `Combined_Preprocess` is reused by different historical scripts; verify which dataset the object contains before saving. Do not rerun both historical scripts into the same workspace and assume identically named objects retain their earlier meaning.

If a final object is absent, Figures 2, 4 and S2H attempt the corresponding raw workflow. S2H raw inputs belong under `data/Figure_S2H/R_control/` and `data/Figure_S2H/R_resolution/`, with the same raw/filtered directory layout. Reprocessing requires validating cell calls, marker identities and cluster mapping before using the results as manuscript replacements.

### Additional inputs for exact reproduction

- `data/annotations/Fig2C_final_genes.csv`: optional exact heatmap row order, with `gene` and `cluster` columns. Without it, the script recalculates top ten significant RNA markers per cluster and exports both the recalculated and displayed gene lists.
- `data/annotations/trajectory_roots.csv`: original root-cell barcodes if archived Monocle objects are unavailable. Columns: `trajectory`, `root_cell`; one row per root cell. An example schema is supplied. The script does not invent a new root and label it as the original analysis.
- `data/annotations/Fig4F_S4F_cells.csv`: one `cell` column specifying the exact cells used for the displayed gene-expression trends. The original code contains both all-cell and cluster-11-specific alternatives, so the selection cannot be inferred safely. The current export uses the spline formula in the revised GitHub script; confirm this against the final figure's fit.
- `data/gmt/m2.cp.reactome.v2025.1.Mm.symbols.gmt`: the mouse Reactome file used by the supplied S2I source. Confirm that this was also the final 5B collection; its exact historical recipe was missing. Obtain GMT files under the applicable [MSigDB terms](https://www.gsea-msigdb.org/gsea/msigdb).
- `data/networks/collectri_mouse.csv`: the final CollecTRI network, with `source`, `target`, `mor`. If absent, Figure 5 downloads the current network and exports it. Archive that exported network and record its provenance/version before the final release.

## Cluster identities

The included map retains the **GitHub Figure 2 convention**: Seurat 0/1/2/3/4 become displayed 1/2/3/4/5. In particular, raw cluster 3 is displayed cluster 4 and raw cluster 1 is displayed cluster 2, matching the source S2F comparison. This map must be checked against the final object, not reapplied blindly after reclustering.

An earlier standalone plotting helper inferred a different ordering from PROGENy colours. That inference is not used here. The new scripts use one explicit map throughout Figures 2, 3 and 5, and stop if existing `figure_cluster` metadata conflict with it. Mouse radar keys use cluster numbers; assigning human macrophage names to these mouse colours without checking the annotations is inappropriate.

## Software and execution

R 4.1 or later is required for the native pipe syntax. Use the actual tested Seurat v5/Monocle/PROGENy environment from the final analysis when available. The manuscript names Seurat 5.0, Monocle3 1.3.7 and PROGENy 1.17.3; these versions have not been installed or runtime-tested for this review.

R dependencies: `Seurat`, `SeuratObject`, `Matrix`, `SoupX`, `DoubletFinder`, `harmony`, `monocle3`, `SeuratWrappers`, `SummarizedExperiment`, `UCell`, `igraph`, `fgsea`, `decoupleR`, `progeny`, `fmsb`, `svglite`, `xml2`, `ggplot2`, `ggrepel`, `patchwork`, `scales`, `dplyr`, `tidyr`, `tibble`, `readr`, `stringr`, `purrr`. `splines`, `grid`, `tools` and `utils` are supplied with R. Python S2E dependencies: Scanpy, AnnData, NumPy, pandas, SciPy, harmonypy, matplotlib, Louvain and igraph.

In RStudio, open a project in the repository folder, or use **Session → Set Working Directory → Choose Directory**. Then run each line separately, resolving any missing inputs before continuing:

```r
source("Figure_2_and_S2_revised.R")
source("Figure_3_mouse_revised.R")
source("Figure_5_mouse_revised.R")
source("Figure_4_and_S4_revised.R")
source("Figure_S2H_revised.R")
```

For S2E:

```r
source("Export_Figure_S2E_input.R")
```

Then in a terminal with the analysis Python environment active:

```bash
python Figure_S2E_scanpy.py
```

Finally in R:

```r
source("Figure_S2E_concordance.R")
```

Files go under `output/Figure_2_S2/`, `output/Figure_3_mouse/`, `output/Figure_4_S4/`, `output/Figure_5_mouse/`, `output/Figure_S2H/` and `output/Figure_S2E/`. Scripts export source tables, cell selections and session information alongside figures. Preserve these outputs in a linked data deposit and capture an environment lockfile from the successful final run.

## Analysis details that affect interpretation

- RNA feature/scatter values are `ln(1 + counts / total counts * 10000)`, dimensionless. Heatmaps use per-gene scaled RNA expression, clipped to −2.5 to +2.5. A colour scale is now retained for 2C.
- Seurat `p_val_adj` in the Wilcoxon marker tables is its Bonferroni-adjusted value; it should not be described as BH FDR. FGSEA uses BH-adjusted pathway P values. Fig. 2F displays raw-P significance, matching its current legend, and exports BH values across all six correlations for review. S2H likewise exports both raw and BH P values.
- The two S2F/I targeted Slamf7 tests use Holm correction across those two tests, as shown in the current figures. The S2I heatmap labels are derived from the actual cell grouping. Ties at a zero lower quartile can make the low group larger than 25% and mean that it consists of undetected cells.
- The S2I heatmap shows 15 significant genes in each direction, selected by adjusted P and then fold change, with up to 200 cells per group. Its export caption states both analysed and displayed cell counts. The manuscript's “top 100 by logFC” sentence does not describe this panel.
- Fig. 3A uses DSS cluster 4 versus cluster 2 as the operational monocyte comparator; Fig. 5B uses DSS cluster 4 versus all other Csf1r+ clusters. These are different comparisons. The Results text describing mouse 5B as inflammation versus health needs reconciliation with the legend and original results.
- CollecTRI ULM and PROGENy use cluster-condition **linear-normalised mean expression**, matching Seurat `AverageExpression(layer="data")`; these are not the mean log values used in the custom dot plots. ULM activity is a model t statistic. The displayed six TFs are prespecified, not an unbiased list of all highest-ranked factors.
- PROGENy uses Mouse, top 500 genes, scaled output and `perm=1`. Scaling uses naive and DSS cluster profiles together; DSS is selected afterwards. Radar rings show each pathway's observed range from 0 to 100%, not a common z-score or percentage biological activation. Raw scores and scaling bounds are exported.
- Radar pathway labels are 7 pt at the exported 58 × 62 mm size. Place SVG/PDF at 100% in Affinity; shrinking it also shrinks the type.
- Figure 4B restores the original fixed 16-gene mouse and 14-gene human-derived signatures. S4B uses independently derived top 20/50/100/200 signatures from the pre-cluster-10-removal object. These analyses are not interchangeable.
- The raw Figure 4 workflow retains >200 genes, <15% mitochondrial reads and Harmony by tissue. Figure 2 uses >500 genes, <10% mitochondrial reads and Harmony by pooled sample. S2H uses >500 genes and <15% mitochondrial reads. Report these separately in Methods if they are the final settings.
- The historical pipelines detect doublets after merging distinct libraries. This is discouraged by the [DoubletFinder developers](https://github.com/chris-mcginnis-ucsf/DoubletFinder). The RNA/SCT flag is corrected here, but a per-library sensitivity analysis would be a further scientific reanalysis, not a formatting change. Original processed objects avoid silently replacing the final cell set.

## Release and citation

See [ZENODO_ARCHIVE_GUIDE.md](ZENODO_ARCHIVE_GUIDE.md) for GitHub upload, a tagged release and a permanent Zenodo software DOI. Resolve the remaining provenance items in [CODE_REVIEW.md](CODE_REVIEW.md), choose the code licence with the code authors and add the verified creator metadata before archiving a manuscript release. No licence, article DOI or code DOI has been invented in this package.

Contact: Gareth-Rhys Jones (`gareth.r.jones@glasgow.ac.uk`) and Calum C. Bain (`calum.bain@glasgow.ac.uk`).
