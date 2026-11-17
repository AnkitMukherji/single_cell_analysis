# Cell clustering and annotation

Open the notebooks in order with the environment YAML in this directory. Inspect input artifact access, paths, and the data layers that each notebook uses. The root README describes annotation model versions, reference data and output fields.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`01_clustering.ipynb`](01_clustering.ipynb) | Uses graph clustering to identify groups of transcriptomically similar cells and inspects the partition in an embedding. |
| [`02_annotation.ipynb`](02_annotation.ipynb) | Compares manual marker labels, CellTypist models and scVI/CellMapper reference transfer. Writes an annotated H5AD and shows how to inspect uncertainty and mark low-confidence cells. |

[Repository overview](../../../README.md)
