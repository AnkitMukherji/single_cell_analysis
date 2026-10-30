import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Preprocessing of scRNA-seq data

    Two essential properties of scRNA-seq datasets:
    - Contains drop-outs i.e., there is an excessive number of zeros in the data due to the limitation of mRNA.
    - The potential for correcting the data and performing quality control might be limited as the data may be confounded with biology.

    The main steps in preprocessing involves:
    - **Doublet detection**
    - **Quality control**
    - **Normalization**
    - **Feature selection**
    - **Dimensionality reduction**

    ## Example analysis on single-cell multiomics data from bone marrow mononuclear cells of 12 healthy human donors measured at four different sites
    """)
    return


@app.cell
def _():
    from typing import Literal

    import anndata as ad
    import lamindb as ln
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import seaborn as sns
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.conversion import localconverter
    from scipy.sparse import csc_matrix
    from scipy.special import logit
    from scipy.stats import median_abs_deviation

    # Suppress verbose logging from Scanpy
    sc.settings.verbosity = 0

    # Set figure parameters for clean, minimal plots
    sc.set_figure_params(dpi=80, facecolor="white", frameon=False)
    return Literal, ad, ln, median_abs_deviation, np, pd, sc, sns


@app.cell
def _(ln):
    # Loading the dataset
    af = ln.Artifact.connect("theislab/sc-best-practices").get(
        key="preprocessing_visualization/quality_control_adata.h5ad", is_latest=True
    )
    adata = af.load()
    adata
    return (adata,)


@app.cell
def _(adata):
    # A warning shows that not all variable names (genes) are unique
    adata.var_names_make_unique() # Makes the variable names unique by appending a number string to each duplicate index element: '1', '2', etc
    adata
    return


@app.cell(hide_code=True)
def _(mo):
    mo.vstack([
        mo.md("## Necessity of quality control"),
        mo.image("../../images/quality_control.jpeg")
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Filtering low quality cells
    - Signs of low quality cells
        - Low number of detected genes
        - A low count depth per cell
        - A high fraction of mitochondrial counts (Reflects cells with broken membranes whose cytoplasmic mRNA has leaked out, and therefore only the mRNA in the mitochondria is still present)
    - Cell QC is performed on the following three QC covariates:
        - The number of counts per barcode (count depth)
        - The number of genes per barcode
        - The fraction of counts from mitochondrial genes per barcode

    - Special cases:
        - Cells with a relatively high fraction of mitochondrial counts might be involved in respiratory processes and should not be filtered out.
        - Cells with low or high counts might correspond to quiescent cell populations or cells larger in size.

    **It is therefore preferred to consider multiple covariates when thresholding decisions on single covariates are made. In general, it is advised to exclude fewer cells and be as permissive as possible to avoid filtering out viable cell populations or small sub-populations.**

    - QC on small datasets can be performed manually by looking at the distribution of different QC covariates and identifying outliers.
    - Large datasets, however, needs automatic thresholding via **MAD (median absolute deviations)**. We can mark cells as outliers if they differ by 5 MADs which is a permissive filtering strategy.
    - It might be reasonable to re-assess the filtering after the annotation of cells.
    - Other QCs:
        - A high fraction of hemoglobin counts points at contamination with red blood cells.
        - The ribosomal fraction varies with cell state and can help explain what separates clusters based on it.
    """)
    return


@app.cell
def _(adata):
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith("RPS", "RPL")
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
    return


@app.cell
def _(adata, sc):
    # Calculating the QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    adata
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Important columns added by the "calculate qc metrics" function
    - **n_genes_by_counts** in **.obs** is the number of genes with positive counts in a cells
    - **total_counts** is the total number of counts for a cell, this might also be known as library size
    - **pct_counts_mt** is the proportion of total counts for a cell which are mitochondrial
    - The boolean **mt**, **ribo** and **hb** flags in **.var** define which genes are summed into **total_count_mt**, **pct_counts_mt** and the corresponding ribosomal and hemoglobin columns.
    """)
    return


@app.cell
def _(adata, sc, sns):
    # Plotting the three QC covariates
    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    p2 = sc.pl.violin(adata, "pct_counts_mt")
    p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - The p2 plot indicates that some cells have a relatively high percentage of mitochondrial counts.
    - Since the number of counts per cell is sufficiently high and the percentage of mitochondrial counts is for most cells below 20%, we can still process the data.
    """)
    return


@app.cell
def _(Literal, ad, median_abs_deviation, np, pd):
    def mad_thresholds(values: pd.Series, nmads:float) -> tuple[float, float]:
        """Thresholds `nmads` median absolute deviations below and above the median of `values`."""
        median, mad = np.median(values), median_abs_deviation(values)
        return median - nmads * mad, median + nmads * mad

    def is_outlier(
        adata: ad.AnnData,
        metric: str,
        nmads: float,
        side: Literal["both","upper"] = "both",
    ) -> pd.Series:
        """Flag the cells beyond the MAD thresholds of `metric`, on both sides or only above,"""
        values = adata.obs[metric]
        lower, upper = mad_thresholds(values, nmads)
        return values > upper if side == "upper" else (values < lower) | (values > upper)

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
