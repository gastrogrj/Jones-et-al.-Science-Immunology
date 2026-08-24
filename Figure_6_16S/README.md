# ACOD1 full-length 16S rRNA gene analysis

## Overview

PacBio HiFi full-length 16S rRNA gene sequencing data were processed using
the Pacific Biosciences
[HiFi-16S workflow](https://github.com/PacificBiosciences/HiFi-16S-workflow).
This Nextflow workflow uses QIIME 2 and DADA2 to generate amplicon sequence
variants (ASVs), assign taxonomy, and produce quality-control summaries and
analysis-ready abundance tables.

The workflow was run using Microsoft Azure Batch with Singularity containers.
Default workflow parameters were retained for data processing unless an
alternative value is shown in the command below. The DADA2 and VSEARCH CPU
allocations were set to 32 and 16, respectively. Azure-specific execution and
storage settings were defined in `azurebatch_hifi.config`.

## Input files

The workflow used the following input files in Azure Blob Storage:

- `acod1_samplesheet.tsv` — sample
  manifest containing the locations of the input sequencing files.
- `acod1_metadata.tsv` — sample
  metadata used for downstream analysis.

## Workflow execution

The repository was cloned and the following command was run from its root
directory, where `main.nf` was located:

```bash
nextflow run main.nf \
    --input az://nextflow/HiFi-16S-workflow/input_files/acod1_samplesheet.tsv \
    --metadata az://nextflow/HiFi-16S-workflow/input_files/acod1_metadata.tsv \
    --dada2_cpu 32 \
    --vsearch_cpu 16 \
    -profile singularity \
    -c azurebatch_hifi.config \
    --outdir az://nextflow/HiFi-16S-workflow/acod1_output
```

The `singularity` profile provided the containerized software environment, and
`azurebatch_hifi.config` specified the Azure Batch executor and associated
cloud-storage configuration. Workflow outputs were written to:

```text
az://nextflow/HiFi-16S-workflow/acod1_output
```

Apart from the resource allocations and execution/storage options specified
above, the workflow's default parameters were used.

## Downstream analysis in R

The merged ASV abundance and taxonomy output from the HiFi-16S workflow was
used for downstream ecological analysis in R with the `phyloseq` package. The
sample metadata file supplied to the Nextflow workflow was also imported into
R and combined with the ASV abundance and taxonomy tables to construct a
`phyloseq` object.

The downstream R script:

1. imports the merged ASV abundance and taxonomic assignments;
2. constructs a `phyloseq` object using the ASV table, taxonomy table, and
   sample metadata;
3. converts counts to relative abundance and applies the specified abundance
   filter;
4. calculates alpha-diversity measures;
5. calculates Bray-Curtis dissimilarities;
6. performs principal coordinates analysis (PCoA) and PERMANOVA; and
7. generates the beta-diversity figure and associated output tables.

The analysis code is provided in `phyloseq_analysis_acod1.R`.

## Reproducibility information

For full computational reproducibility, the following files should accompany
this README and the R analysis script:

- `azurebatch_hifi.config`, with credentials removed;
- `acod1_samplesheet.tsv`;
- `acod1_metadata.tsv`

## Software

- [Pacific Biosciences HiFi-16S workflow (v0.9)](https://github.com/PacificBiosciences/HiFi-16S-workflow)
- [Nextflow (v25.04.6)](https://www.nextflow.io/)
- [phyloseq](https://joey711.github.io/phyloseq/)
