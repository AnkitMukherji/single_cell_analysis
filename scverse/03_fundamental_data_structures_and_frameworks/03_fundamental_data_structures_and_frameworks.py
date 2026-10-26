import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", auto_download=["html"])


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

    # ln.track("zmuES9DakbmX")
    return ad, csr_matrix, ln, np, pd


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
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    print(adata.obs_names[:10])
    return


@app.cell
def _(adata, np, pd):
    # Adding aligned metadata
    # Observational (in .obs) or variable level (in .var)
    ct = np.random.default_rng().choice(["B", "T", "Monocyte"], size=adata.n_obs)
    adata.obs["cell_type"] = pd.Categorical(ct) # Categoricals are preferred for efficiency
    adata.obs
    return


@app.cell
def _(adata):
    adata
    return


@app.cell
def _(adata):
    # Subsetting using metadata
    bdata = adata[adata.obs.cell_type == "B"]
    bdata
    return


@app.cell
def _(adata, np):
    # Observation/Variable-level matrices stored in .obsm/.varm
    # .obsm matrices must have a length equal to the number of observations as .n_obs
    # .varm matrices must have a length equal to the number of variables as .n_vars

    # Example of storing UMAP embedding of the data
    adata.obsm["X_umap"] = np.random.default_rng().normal(0, 1, size=(adata.n_obs, 2))

    # Adding random gene-level metadata
    adata.varm["gene_stuff"] = np.random.default_rng().normal(0, 1, size=(adata.n_vars, 2))

    adata.obsm
    return


@app.cell
def _(adata):
    adata
    return


@app.cell
def _(adata):
    # Unstructured metadata
    adata.uns["random"] = [1, 2, 3]
    adata.uns
    return


@app.cell
def _(adata, np):
    # Layers
    adata.layers["log_transformed"] = np.log1p(adata.X)
    adata
    return


@app.cell
def _(adata):
    # Comparing original matrix X with log-transformed X
    # .nnz returns the number of non-zero elements in the boolean matrix
    (adata.X != adata.layers["log_transformed"]).nnz == 0 # The inner comparison counts the number of true instances which is returned by .nnz
    return


@app.cell
def _(adata):
    # Conversion to DataFrames
    adata.to_df(layer="log_transformed")
    return


@app.cell
def _(adata):
    # Readin and writing of AnnData objects
    # Can be saved in HDF5-based file format .h5ad
    # AnnData autotransforms string data into categorical variables
    adata.write("my_results.h5ad", compression="gzip")
    return


@app.cell
def _(ad):
    adata_new = ad.read_h5ad("my_results.h5ad")
    adata_new
    return


@app.cell
def _(adata, np, pd):
    # Efficient Data Access
    # Observations come from instruments characterizing 10 readouts in a multi-year study with samples taken from different subjects at different sites
    obs_meta = pd.DataFrame(
        {
            "time_yr": np.random.default_rng().choice([0, 2, 4, 8], adata.n_obs),
            "subject_id": np.random.default_rng().choice(
                ["subject 1", "subject 2", "subject 4", "subject 8"], adata.n_obs
            ),
            "instrument_type": np.random.default_rng().choice(
                ["type a", "type b"], adata.n_obs
            ),
            "site": np.random.default_rng().choice(["site x", "site y"], adata.n_obs)
        },
        index=adata.obs.index
    )
    return (obs_meta,)


@app.cell
def _(obs_meta):
    obs_meta
    return


@app.cell
def _(ad, adata, obs_meta):
    adata_2 = ad.AnnData(adata.X, obs=obs_meta, var=adata.var)
    adata_2
    return (adata_2,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Subsetting AnnData objects always returns views, which has two advantages:
    - No new memory is allocated.
    - It is possible to modify the underlying AnnData object
    """)
    return


@app.cell
def _(adata_2):
    # Indexing into AnnData will assume that integer arguments to [] behave like .iloc in pandas
    # whereas string arguments behave like .loc

    # Example of view
    adata_view = adata_2[:5, ["Gene_1", "Gene_3"]]
    adata_view
    return


@app.cell
def _(adata_2):
    adata_2
    return


@app.cell
def _(adata_2):
    # For a view, we can also set the first three elements of a column
    print(adata_2[:3, "Gene_1"].X.toarray().tolist())
    adata_2[:3, "Gene_1"].X = [0, 0, 0]
    print(adata_2[:3, "Gene_1"].X.toarray().tolist())
    return


@app.cell
def _(adata_2):
    adata_subset = adata_2[:3, ["Gene_1", "Gene_2"]]
    adata_subset
    return (adata_subset,)


@app.cell
def _(adata_subset):
    adata_subset.obs["foo"] = range(3)
    adata_subset
    return


@app.cell
def _(adata_2):
    # Using pandas to slice with sequences or boolean indices
    adata_2[adata_2.obs.time_yr.isin([2, 4])].obs.head()
    return


@app.cell
def _(ad):
    # Partial reading of large data
    adata_partial = ad.read_h5ad("my_results.h5ad", backed="r")
    return (adata_partial,)


@app.cell
def _(adata_partial):
    adata_partial.isbacked
    return


@app.cell
def _(adata_partial):
    # The AnnData object has an open connection to the file used for reading
    adata_partial.filename
    return


@app.cell
def _(adata_partial):
    # Always explicitly close the connection
    adata_partial.file.close()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Scanpy API design

    - All preprocessing functions are available in the **scanpy.preprocessing** module (**scanpy.pp**).
    - All transformations of a data matrix that are not preprocessing are available in **scanpy.tools** (**scanpy.tl**).
    - All visualizations are available in **scanpy.plot** (**scanpy.pl**).
    - **scanpy.datasets** stores various datasets.
    """)
    return


@app.cell
def _(mo):
    mo.vstack([
        mo.md("# Scanpy API overview"),
        mo.image("../images/scanpy_api.png")
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Example workflow of an analysis using scanpy
    """)
    return


@app.cell
def _():
    import scanpy as sc

    sc.set_figure_params(dpi=80, facecolor="white")
    sc.logging.print_header()
    return (sc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Dataset

    2700 peripheral blood mononuclear cells of a healthy donor sequenced on the Illumina NextSeq 500
    """)
    return


@app.cell
def _(ln):
    adata_example = ln.Artifact.get(
        key="introduction/fundamental_data_structures_and_frameworks.h5ad", is_latest=True
    ).load()
    adata_example
    return (adata_example,)


@app.cell
def _(adata_example):
    # The var slow contains the gene IDs
    adata_example.var
    return


@app.cell
def _(adata_example, sc):
    # Genes that yield the highest fraction of counts in each single cell, across all cells
    sc.pl.highest_expr_genes(adata_example, n_top=20)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    *MALAT1* is the most expressed gene which is frequently detected in poly-A captured scRNA-Seq data, independent of protocol. This gene has been shown to have an inverse correlation with cell health. Especially dead/dying cells have a higher expression of *MALAT1*.
    """)
    return


@app.cell
def _(adata_example):
    adata_example
    return


@app.cell
def _(adata_example, sc):
    # Filter cells with less than 200 detected genes
    # Filter genes found in less than 3 cells
    sc.pp.filter_cells(adata_example, min_genes=200)
    sc.pp.filter_genes(adata_example, min_cells=3)
    return


@app.cell
def _(adata_example):
    adata_example
    return


@app.cell
def _(adata_example, sc):
    # PCA
    sc.tl.pca(adata_example, svd_solver="arpack")
    return


@app.cell
def _(adata_example):
    adata_example
    return


@app.cell
def _(adata_example, sc):
    # Plotting a specific gene
    sc.pl.pca(adata_example, color="CST3")
    return


@app.cell
def _(adata_example, sc):
    # Calculating neighborhood graph using PCA
    # Used by other tools such as calculation of a UMAP
    sc.pp.neighbors(adata_example, n_neighbors=10, n_pcs=40)
    return


@app.cell
def _(adata_example):
    adata_example
    return


@app.cell
def _(adata_example, sc):
    sc.tl.umap(adata_example)
    return


@app.cell
def _(adata_example):
    adata_example
    return


@app.cell
def _(adata_example, sc):
    sc.pl.umap(adata_example, color=["CST3", "NKG7", "PPBP"])
    return


@app.cell
def _(adata_example):
    adata_example
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Therefore,

    ### The fundamental data structure for single-cell analysis in the scverse is **AnnData**.

    ### The fundamental framework for single-cell analysis in the scverse is **Scanpy**.
    """)
    return


if __name__ == "__main__":
    app.run()
