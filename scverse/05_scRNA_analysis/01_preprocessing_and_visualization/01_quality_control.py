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
    return Literal, ad, ln, logit, median_abs_deviation, np, pd, plt, sc, sns


@app.cell
def _(ln):
    # Loading the dataset
    af = ln.Artifact.connect("theislab/sc-best-practices").get(
        key="preprocessing_visualization/quality_control_adata.h5ad", is_latest=True
    )
    adata = af.load()
    adata
    return adata, af


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

    return is_outlier, mad_thresholds


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - The bounded and skewed mitochondrial fraction is logit-transformed before computing MADs.
    - A pseudocount keeps cells without mitochondrial counts finite.
    """)
    return


@app.cell
def _(adata, logit):
    mt_fraction = (adata.obs["total_counts_mt"] + 1) / (adata.obs["total_counts"] + 2)
    adata.obs["logit_mt_fraction"] = logit(mt_fraction)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - Using 5 MADs for total_counts, n_genes_by_counts and pct_counts
    - For mitochondrial fraction, flag cells that exceed both 3 MADs and a fixed cutoff of 8% which is not universal since mitochondrial content differs between cell types and tissues
    """)
    return


@app.cell
def _(adata, is_outlier):
    qc_filters = {
        "log1p_total_counts": (5, "both"),
        "log1p_n_genes_by_counts": (5, "both"),
        "pct_counts_in_top_20_genes": (5, "both"),
        "logit_mt_fraction":(3, "upper"),
    }

    flags = {
        metric: is_outlier(adata, metric, nmads, side)
        for metric, (nmads, side) in qc_filters.items()
    }
    adata.obs["outlier"] = (
        flags["log1p_total_counts"]
        | flags["log1p_n_genes_by_counts"]
        | flags["pct_counts_in_top_20_genes"]
    )
    adata.obs["mt_outlier"] = flags["logit_mt_fraction"] & (adata.obs["pct_counts_mt"] > 8)
    adata.obs[["outlier", "mt_outlier"]].value_counts()
    return flags, qc_filters


@app.cell
def _(flags):
    flags
    return


@app.cell
def _(adata, logit, mad_thresholds, plt, qc_filters, sns):
    # Plotting the thresholds
    fig, axes = plt.subplots(1, len(qc_filters), figsize=(16, 3.5))

    for ax, (metric, (nmads, side)) in zip(axes, qc_filters.items(), strict=True):
        sns.histplot(adata.obs[metric], bins=100, ax=ax)
        lower, upper = mad_thresholds(adata.obs[metric], nmads)
        for threshold in [upper] if side == "upper" else [lower, upper]:
            ax.axvline(threshold, color="red")
    axes[-1].axvline(logit(0.08), color="red", linestyle="--") # 8 % cutoff
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(adata):
    # Filtering the AnnData object
    print(f"Total number of cells: {adata.n_obs}")
    adata_filtered = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    print(f"Number of cells after filtering of low quality cells: {adata_filtered.n_obs}")
    return (adata_filtered,)


@app.cell
def _(adata_filtered, sc):
    pl = sc.pl.scatter(adata_filtered, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Correction of ambient RNA
    - Only applies to **droplet-based protocols**.
    - A certain amount of background mRNAs (cell-free mRNAs) is present in the dilution that gets distributed into the droplets with cells and sequenced along with them.
    - Methods for removal of ambient mRNA like **SoupX** and **DecontX** aim to estimate the composition of the soup (term for contamination of cell-free mRNA in the input solution) and correct the count table with respect to the soup-based expression.

    N.B.: SoupX works better if a basic clustering is provided.
    """)
    return


@app.cell
def _(adata_filtered, sc):
    # Normalize and lop1p transform the data
    adata_pp = adata_filtered.copy()
    sc.pp.normalize_total(adata_pp, target_sum=1e4)
    sc.pp.log1p(adata_pp)
    return (adata_pp,)


@app.cell
def _(adata_filtered, adata_pp, sc):
    # Computing the principal components of the data
    # which in turn is used to generate a neighborhood graph
    # and then run leiden clustering on the KNN-graph
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(
        adata_pp, key_added="soupx_groups", flavor="igraph", n_iterations=2, directed=False
    )

    # Preprocess variables for SoupX
    adata_filtered.obs["soupx_groups"] = adata_pp.obs["soupx_groups"]
    return


@app.cell
def _(adata_filtered):
    adata_filtered
    return


@app.cell
def _(adata_filtered):
    cells = adata_filtered.obs_names
    genes = adata_filtered.var_names

    # SoupX requires a matrix of shape features x barcodes, so we have to transpose .X
    data = adata_filtered.X.T
    return


@app.cell
def _(af):
    # SoupX also requires the raw gene by cells matrix typically called raw_feature_bc_matrix.h5 in the cellranger output
    adata_raw = af.load()
    adata_raw.var_names_make_unique()

    genes_raw = adata_raw.var_names
    cells_raw = adata_raw.obs_names

    data_tod = adata_raw.X.T
    return


@app.cell
def _():
    # Converting the data into the right structure to be used in R

    return


if __name__ == "__main__":
    app.run()
