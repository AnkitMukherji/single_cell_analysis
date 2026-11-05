# Single-cell Analysis in Python and R

A practical collection of single-cell analysis scripts, notebooks, and Quarto notes spanning **scverse** and **Bioconductor**.

## Analysis pathway

| Area | Material |
| --- | --- |
| Technologies and quantification | [Technology basics](scverse/01_basics_single_cell_technology.py), [alevin/simpleaf example](scverse/02_alevin_simpleaf_example_processing.sh) |
| Data representations | [Core structures](scverse/03_fundamental_data_structures_and_frameworks/), [multimodal and spatial structures](scverse/04_multimodal_and_spatial_data_structures/) |
| Preprocessing | [QC, normalization, feature selection, and dimensionality reduction](scverse/05_scRNA_analysis/01_preprocessing_and_visualization/) |
| Cellular structure | [Clustering and annotation](scverse/05_scRNA_analysis/02_identifying_cellular_structure/) |
| R ecosystem | [Bioconductor/OSCA Quarto notes](bioconductorverse/scripts/) |

## Environments and use

Python environments are recorded per topic rather than as one shared lockfile. For example:

```bash
conda env create -f scverse/05_scRNA_analysis/01_preprocessing_and_visualization/preprocessing.yml
```

Activate the environment name declared in the selected YAML file, then open the relevant notebook. Clustering and annotation have separate environment files. Inspect data download requirements and local paths before execution.

The R side has its own [renv.lock](bioconductorverse/renv.lock). From the repository root in R:

```r
renv::restore(project = "bioconductorverse")
```

Open the `.qmd` sources in Quarto-compatible tooling. These are learning and workflow-development materials; availability of an environment file does not imply a fresh-machine execution has been verified.

## Interpretation focus

Use QC to understand compromised cells and technical variation before clustering. Evaluate normalization, feature selection, and embedding choices together. A visually separated cluster is not sufficient evidence for a distinct cell type: inspect markers, donor/batch composition, and biological context.

Spatial and multimodal data-structure examples establish representation and access patterns; their presence does not imply a completed spatial research study.

Related: [single-cell paper replication](https://github.com/AnkitMukherji/Paper_Replication) · [statistical genetics](https://github.com/AnkitMukherji/statistical_genetics).
