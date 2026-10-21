import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.vstack([
        mo.md("# Fundamental Data Structures"),
        mo.image("../images/fundamental_data_str.png")
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Storing unimodal data with AnnData

    - AnnData supports sparse matrices and partial reading.
    - At the core, AnnData stores **X**, a sparse or dense matrix (count matrix in scRNA-Seq). The matrix has the dimensions of **obs_names x var_names**, where the obs (=observations) correspond to the cell's barcodes and the var (=variables) correspond to the gene identifiers. X in different forms (e.g., raw and normalized count matrix) can be stored in layers.
    - The matrix X is surrounded by Pandas DataFrames obs and var, which save annotations of cells and genes, respectively.
    - AnnData saves whole matrices of calculations for the observations in **obsm** or variables in **varm** with their corresponding dimensions.
    - Graph-like structures that associate cells with cells or genes with genes are usually saved in **obsp** and **varp**.
    - Any other unstructured data which does not fit any other slot is saved as unstructured data in **uns**.
    """)
    return


@app.cell
def _(mo):
    mo.vstack([
        mo.md("# AnnData structure"),
        mo.image("../images/anndata.jpg")
    ])
    return


@app.cell
def _():
    import anndata as ad
    import lamindb as ln
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    ln.connect("theislab/sc-best-practices")

    ln.track("zmuES9DakbmX")
    return ad, csr_matrix, np


@app.cell
def _(ad, csr_matrix, np):
    # Initializing an AnnData object with random Poisson distribution
    counts = csr_matrix(np.random.default_rng().poisson(1, size=(100, 2000)), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata
    return (adata,)


@app.cell
def _(adata):
    # 100 cells x 2000 genes
    adata.X
    return


@app.cell
def _(adata):
    # Providing the index to both obs and var axes using .obs.names and .var.names
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var.names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    print(adata.obs_names[:10])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
