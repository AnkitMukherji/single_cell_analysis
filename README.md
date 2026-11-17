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

## Cell-type annotation: markers, classifiers, and reference mapping

The [annotation notebook](scverse/05_scRNA_analysis/02_identifying_cellular_structure/02_annotation.ipynb) compares three ways to assign cell identities to an example from the Single-cell Best Practices teaching collection. It starts with already processed single-cell counts and shows how each method produces labels, where the labels disagree, and how to inspect uncertain assignments.

### Setup and inputs

From the repository root:

```bash
conda env create -f scverse/05_scRNA_analysis/02_identifying_cellular_structure/annotation.yml
conda activate annotation
```

Open the notebook using the `annotation` Python environment as its Jupyter kernel. Set the notebook working directory to its own folder to reproduce the relative output paths. The [environment specification](scverse/05_scRNA_analysis/02_identifying_cellular_structure/annotation.yml) records Python 3.14.7, Scanpy 1.12.3, scvi-tools 1.5.0.post1, CellTypist 1.6.3, LaminDB 2.10.0 and CellMapper 0.2.6.

The notebook connects to the public teaching instance `theislab/sc-best-practices` in LaminDB. It retrieves `cellular_structure/s4d8_clustered.h5ad`, the reference feature table, a pretrained bone marrow scVI model, and the corresponding labeled reference embedding. The example query object supplies raw counts, scran-normalized expression, gene identifiers and highly deviant gene flags. Network access is required for these artifacts and the CellTypist model downloads. The artifact lookups request the latest versions. Record the resolved artifact versions and model files when comparing runs.

### Run the notebook in order

| Method | What the notebook does | What to inspect |
| --- | --- | --- |
| Manual annotation | Starts with scran-normalized expression, builds a PCA neighbor graph and Leiden clusters, then inspects marker expression in UMAPs and dot plots. It also ranks cluster markers with a Wilcoxon test and filters them by expression inside and outside the cluster. | Marker combinations, expression outside the target cluster, and provisional or missing labels in `manual_celltype_annotation`. |
| CellTypist | Copies raw counts, normalizes to 10,000 counts per cell, applies log1p, and converts the matrix to dense format. It downloads coarse and fine immune models and applies majority voting. | The coarse and fine label columns and their confidence scores in `adata.obs`. The dense copy can require substantial memory on larger datasets. |
| scArches-style mapping with scVI and CellMapper | Uses raw counts and matches features by gene ID. Missing reference features are filled with zeros, and columns are ordered to match the reference model. It adapts the pretrained scVI model to the query and uses the latent representation for mapping. | Gene overlap and order, reference/query coverage, marker agreement, and label uncertainty. The notebook includes a compatibility shim for loading the serialized reference model with newer pandas. |
| Label transfer and review | CellMapper maps labels from 15 reference neighbors with the scArches kernel. The notebook stores uncertainty as one minus the returned confidence and marks scores above 0.2 as `Unknown`. It checks selected assignments against marker expression. | Compare `transf_cell_type`, `transf_cell_type_unc` and `transf_cell_type_certain`. The 0.2 cutoff is an example threshold for this teaching workflow. |

The final cell writes `../intermediate_data/annotation.h5ad` relative to the notebook working directory. The file includes annotation columns in `obs`, the query scVI embedding in `obsm["X_scVI"]`, and marker ranking information in `uns`. The notebook also displays marker dot plots, UMAPs, and a box plot of label-transfer uncertainty. Inspect the notebook for the example outputs. This review checked the notebook source and package specification without rerunning the training.

### Interpreting the three approaches

Read marker evidence and reference coverage together. Marker-based labels depend on the cluster resolution and marker specificity. Rankings derived from the same data used to define clusters are exploratory. Classifier labels depend on the training data and its label definitions. Reference mapping can assign a familiar label to a state missing from the reference, so check uncertainty and marker agreement before interpreting a subtype. UMAP mixing is a useful visual diagnostic for the mapping and needs supporting checks of biological plausibility.

## Interpretation focus

Use QC to understand compromised cells and technical variation before clustering. Evaluate normalization, feature selection, and embedding choices together. A visually separated cluster is not sufficient evidence for a distinct cell type: inspect markers, donor/batch composition, and biological context.

Spatial and multimodal data-structure examples establish representation and access patterns; their presence does not imply a completed spatial research study.

Related: [single-cell paper replication](https://github.com/AnkitMukherji/Paper_Replication) · [statistical genetics](https://github.com/AnkitMukherji/statistical_genetics).
