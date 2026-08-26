# Jones et al. - Science Immunology analysis code

This repository contains analysis code accompanying:

> Jones G-R, Drury B, Alegbe T, et al. **Concurrent induction of proinflammatory and regulatory programmes in inflammation-associated monocyte-derived cells in inflammatory bowel disease.** *Science Immunology* (in press).

The study identifies an inflammation-associated monocyte-derived cell population in human inflammatory bowel disease and experimental colitis, and examines the concurrent induction of inflammatory and ACOD1-dependent regulatory programmes.

## Repository contents

| File | Analysis covered |
| --- | --- |
| `Figure_2_and_S2_revised.R` | Mouse colon scRNA-seq preprocessing and analyses for Fig. 2A-F and Fig. S2C/D/F/H/I |
| `Figure_4_and_S4_revised.R` | Combined mouse blood/colon scRNA-seq preprocessing, signature scoring and trajectory analyses for Fig. 4A/B/D-F and Fig. S4B/E/F |
| `Figure_5_human_revised.R` | Human scRNA-seq Reactome GSEA, transcription-factor activity and PROGENy analyses for the human components of Fig. 5A/C/D |
| `S2_python_mouse.ipynb` | Scanpy reanalysis used for Fig. S2E |
| `Figure_6D_16S/` | PacBio HiFi full-length 16S workflow documentation, sample metadata and downstream phyloseq analysis for Fig. 6D |

The existing human scRNA-seq scripts for Fig. 1A/B/D/E, Fig. 3D/E and related supplementary panels are maintained separately and are not duplicated in these files.

Flow-cytometry panels were analysed in FlowJo and/or GraphPad Prism and are not recreated by these R scripts.

## Data availability

| Dataset | Repository and accession | Status at 26 August 2026 |
| --- | --- | --- |
| Human intestinal scRNA-seq | [Zenodo DOI 10.5281/zenodo.8301000](https://doi.org/10.5281/zenodo.8301000) | Available |
| Mouse blood and colon scRNA-seq | [NCBI GEO: GSE345123](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345123) | Accession assigned; record currently private |
| Mouse colonic monocyte NanoString profiling | [NCBI GEO: GSE345013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345013) | Private until 27 August 2027 |
| Mouse faecal microbiota full-length 16S rRNA sequencing | [NCBI BioProject: PRJNA1518240](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1518240) | Private until 1 September 2027; 12 BioSamples and 12 SRA records |

Accessions under embargo may not display a public record until their scheduled release. NCBI release dates can be brought forward when the article is published.

Large raw sequencing files and processed expression matrices are not stored directly in this GitHub repository. After downloading the data, arrange the scRNA-seq files as described below. The Figure 6D subfolder contains its own processing notes and metadata, but not the deposited FASTQ files.

## Expected local data structure

Run all scripts from the repository root. The expected local structure is:

```text
Jones-et-al.-Science-Immunology/
├── README.md
├── Figure_2_and_S2_revised.R
├── Figure_4_and_S4_revised.R
├── Figure_5_human_revised.R
├── S2_python_mouse.ipynb
├── Figure_6D_16S/
│   ├── README.md
│   ├── phyloseq_analysis_acod1.R
│   ├── azurebatch_hifi.config
│   ├── input_files/
│   │   ├── acod1_samplesheet.tsv
│   │   └── acod1_metadata.tsv
│   └── acod1_output/
│       └── results/
│           └── best_tax_merged_freq_tax.tsv
├── data/
│   ├── mouse_scRNAseq/
│   │   ├── G1/
│   │   │   ├── raw_feature_bc_matrix/
│   │   │   └── filtered_feature_bc_matrix/
│   │   ├── G2/
│   │   │   ├── raw_feature_bc_matrix/
│   │   │   └── filtered_feature_bc_matrix/
│   │   ├── B1/
│   │   │   ├── raw_feature_bc_matrix/
│   │   │   └── filtered_feature_bc_matrix/
│   │   └── B2/
│   │       ├── raw_feature_bc_matrix/
│   │       └── filtered_feature_bc_matrix/
│   ├── human_scRNAseq/
│   │   └── fr003_myeloid.celltype_label.esmu.csv.gz
│   ├── objects/
│   │   └── recovery_csf1r.rds
│   └── gmt/
│       ├── m2.cp.reactome.v2025.1.Mm.symbols.gmt
│       └── c2.cp.reactome.v2023.2.Hs.symbols.gmt
└── output/
```

The downstream Figure 6D script additionally expects the HiFi-16S workflow output `Figure_6D_16S/acod1_output/results/best_tax_merged_freq_tax.tsv`. This derived abundance/taxonomy table should be generated from the deposited 16S reads using the workflow described in `Figure_6D_16S/README.md`. Set `project_dir` near the top of `phyloseq_analysis_acod1.R` to the local `Figure_6D_16S` directory before running it. Before committing `azurebatch_hifi.config`, remove all cloud account names, keys, SAS tokens, passwords and other credentials.

Dataset labels used by the mouse scripts are:

| Directory | Tissue | Condition |
| --- | --- | --- |
| `G1` | Colon | Naive |
| `G2` | Colon | Acute DSS colitis |
| `B1` | Blood | Naive |
| `B2` | Blood | Acute DSS colitis |

`recovery_csf1r.rds` is the processed Seurat object used for the independent recovery/replication analysis in Fig. S2H. Its upstream preprocessing was not part of the historical Figure 2 script.

The Reactome GMT files can be obtained from the [Molecular Signatures Database](https://www.gsea-msigdb.org/gsea/msigdb). Users should comply with the applicable MSigDB terms when downloading and using these files.

## Software requirements

The scripts were prepared for R 4.x and use Seurat v5-style objects. Required packages include:

- `Seurat`
- `SoupX`
- `DoubletFinder`
- `harmony`
- `monocle3`
- `SeuratWrappers`
- `UCell`
- `fgsea`
- `decoupleR`
- `progeny`
- `fmsb`
- `patchwork`
- `dplyr`, `tidyr`, `tibble`, `readr`, `stringr`, `purrr`
- `ggplot2`, `ggrepel`, `scales`, `pheatmap`

The Figure 6D analysis additionally uses `phyloseq`, `microbiome`, `vegan`, `pairwiseAdonis`, `ggExtra`, `cowplot` and `rlang`. Reprocessing the raw 16S reads requires the [PacBio HiFi-16S workflow v0.9](https://github.com/PacificBiosciences/HiFi-16S-workflow), Nextflow v25.04.6 and a compatible Singularity environment; see the subfolder README for the exact command and configuration.

The manuscript specifies Seurat v5.0, Monocle3 v1.3.7 and PROGENy v1.17.3 for the corresponding analyses. Each of the three scRNA-seq R scripts writes `sessionInfo.txt` to its output directory so that the executed package versions can be recorded.

## Running the analyses

After installing the required packages and downloading the input data, run:

```bash
Rscript Figure_2_and_S2_revised.R
Rscript Figure_4_and_S4_revised.R
Rscript Figure_5_human_revised.R
Rscript Figure_6D_16S/phyloseq_analysis_acod1.R
```

Outputs are written to:

```text
output/Figure_2_S2/
output/Figure_4_S4/
output/Figure_5_human/
Figure_6D_16S/figures/
```

The scRNA-seq scripts export panel-level plots, intermediate result tables, processed Seurat/Monocle objects where relevant, and session information. Input checks stop execution with an informative message when a required file, column, cluster or gene is missing. The 16S script exports alpha-diversity results and the Figure 6D beta-diversity plot.

## Panel provenance

### Figure 2 and Supplementary Figure 2

- Fig. 2A-F: generated by `Figure_2_and_S2_revised.R`.
- Fig. 2G-I: flow cytometry; not generated in R.
- Fig. S2A/B/G: flow cytometry; not generated in R.
- Fig. S2C/D/F/I: generated by `Figure_2_and_S2_revised.R`.
- Fig. S2E: generated using the separate Scanpy notebook.
- Fig. S2H: generated from the processed recovery/replication Seurat object.

### Figure 4 and Supplementary Figure 4

- Fig. 4A/B/D-F: generated by `Figure_4_and_S4_revised.R`.
- Fig. 4C/G: flow cytometry; not generated in R.
- Fig. S4B/E/F: generated by `Figure_4_and_S4_revised.R`.
- Fig. S4A/C/D/G: flow cytometry; not generated in R.

### Figure 5 human scRNA-seq analyses

- Human components of Fig. 5A/C/D: generated by `Figure_5_human_revised.R`.
- Mouse and experimental components of Fig. 5 are outside the scope of this script.

### Figure 6D full-length 16S rRNA analysis

- Raw PacBio HiFi reads were processed with the PacBio HiFi-16S Nextflow workflow.
- Fig. 6D beta-diversity ordination and PERMANOVA are generated by `Figure_6D_16S/phyloseq_analysis_acod1.R` from the merged ASV abundance/taxonomy table and sample metadata.
- The same script also exports alpha-diversity results as a CSV file.

## Reproducibility notes

- Figure 2 cluster numbers displayed in the manuscript are one greater than the original zero-indexed Seurat identities; this mapping is explicit in the script.
- The Figure 2 correlation panel uses cell-level Spearman rank correlations based on normalised RNA expression.
- The DoubletFinder `sct = TRUE` setting is retained from the supplied historical pipelines and is declared explicitly in both mouse scripts.
- The Figure 4 workflow retains Harmony integration using `tissue`, matching the supplied analysis code.
- The human-derived Figure 4 signature contains the 14 mouse orthologues retained in the supplied analysis script; these are exported as a separate CSV when the analysis runs.
- Monocle3 trajectory roots are selected reproducibly from the highest-*Ly6c2*-expressing blood cells, or the highest-*Ly6c2*-expressing colon cells for colon-only analyses.
- Figure 5 human Reactome results are filtered at FDR < 0.05, consistent with the figure legend.

## Citation

If you use this code or the associated datasets, please cite the final published article. Replace the placeholder below when the DOI is available:

```text
Jones G-R, Drury B, Alegbe T, et al. Concurrent induction of proinflammatory
and regulatory programmes in inflammation-associated monocyte-derived cells
in inflammatory bowel disease. Science Immunology. 2026. DOI: [add DOI]
```

## Contact

For questions about the analysis, please contact:

- Gareth-Rhys Jones - `gareth.r.jones@glasgow.ac.uk`
- Calum C. Bain - `calum.bain@glasgow.ac.uk`
