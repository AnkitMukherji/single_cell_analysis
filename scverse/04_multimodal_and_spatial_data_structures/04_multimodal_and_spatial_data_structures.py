import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", auto_download=["html"])


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.vstack([
        mo.md("# Multimodal and Spatial Data Structures"),
        mo.image("../images/advanced_data_struct.png")
    ])
    return


@app.cell
def _(mo):
    mo.vstack([
        mo.md("# Using MuData to store multimodal data"),
        mo.md("""
        - MuData builds on top of AnnData to store and manipulate multimodal data.
        - MuData object contains references to single AnnData objects of the unimodal data where each AnnData object is a dictionary stored in **.mod** (=modality) attribute
            """),
        mo.image("../images/mudata_overview.png")
    ])
    return


@app.cell
def _():
    import lamindb as ln
    import mudata as md
    import numpy as np

    ln.connect("theislab/sc-best-practices")
    return ln, md, np


@app.cell
def _(ln):
    adata = ln.Artifact.get(
        key="introduction/advanced_data_structures_and_frameworks_mudata_1.h5ad",
    ).load()
    adata
    return (adata,)


@app.cell
def _(ln):
    adata2 = ln.Artifact.get(
        key="introduction/advanced_data_structures_and_frameworks_mudata_2.h5ad",
    ).load()
    adata2
    return (adata2,)


@app.cell
def _(adata, adata2, md):
    # The two AnnData objects (two "modalities") can be wrapped into a single MuData object
    # Naming modality one A and modality two B
    mdata = md.MuData({"A": adata, "B": adata2})
    return (mdata,)


@app.cell
def _(mdata):
    # Observations with identical names in different modalities are considered to be the same observation
    # This indicates that variable names should be unique
    # Therefore, in mdata there are 1000 observations and 150 = 100 + 50 variables
    mdata.shape
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## MuData attributes

    ### .mod
    """)
    return


@app.cell
def _(mdata):
    # Modalities are stored in a collection accessible via the .mod attribute of the MuData object with names of modalities as keys and AnnData objects as values
    list(mdata.mod.keys())
    return


@app.cell
def _(mdata):
    # Accessing individual modalities
    print(mdata.mod["A"])
    print(mdata["A"])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### .obs and .var

    - Sample (cells) annotation is accessible via the **.obs** attribute and includes **.obs** data from individual modalities.
    - Annotations of variables (features) are stored in **.var**.
    - Observatons column from individual modalities contain modality name as their prefix e.g., rna:n_genes. This is also true for variables column.
    - If there are columns with identical names in **.var** of multiple modalities, these columns are merged across modalities and no prefix is added. In such cases, when new columns are added or samples (cells) are filtered out, the changes have to be fetched with the **.update()** method.
    """)
    return


@app.cell
def _(mdata):
    mdata.var_names
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - Multidimensional annotations of sample (cells) are accessible in the **.obsm** attribute.
    - If the shape of a modality is changed, **MuData.update()** has to be run to bring the respective updates to the MuData object.
    """)
    return


@app.cell
def _(adata2, mdata):
    adata2.var_names = ["var_ad2_" + e.split("_")[1] for e in adata2.var_names]
    print("Outdated variables names: ...,", ", ".join(mdata.var_names[-3:]))
    mdata.update()
    print("Updated variables names: ...,", ", ".join(mdata.var_names[-3:]))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - Individual modalities are stored as references to the original objects. Hence, if the original AnnData object is changed the change will also be reflected in the MuData object
    """)
    return


@app.cell
def _(adata):
    # Adding some unstructured data to the original object
    adata.uns["misc"] = {"adata": True}
    return


@app.cell
def _(mdata):
    # Access modality A via the .mod attribute
    mdata.mod["A"]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Variable mappings
    """)
    return


@app.cell
def _(mdata):
    # Since all the observations are the same across modalities in mdata, all the values in the observations mappings are set to True
    mdata.obsm["A"]
    return


@app.cell
def _(mdata, np):
    np.sum(mdata.obsm["A"]) == np.sum(mdata.obsm["B"]) == 1000
    return


@app.cell
def _(mdata):
    # For variables, there are 150-long vectors where the A modality has 100 True values followed by 50 False values
    mdata.varm["A"]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## MuData views
    """)
    return


@app.cell
def _(mdata):
    # Analogous to the behaviour of AnnData objects, slicing MuData objects returns views of the original data
    view = mdata[:100, :1000]
    print(view)
    return (view,)


@app.cell
def _(view):
    # Subsetting MuData objects slices them across modalities
    # A copy of MuData object is not a view and has no dependence on the original object
    mdata_sub = view.copy()
    mdata_sub.is_view
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading and Writing of MuData objects
    """)
    return


@app.cell
def _(md, mdata):
    # Similar to AnnData objects, MuData objects are also designed to be serialized into HDF5 based .h5mu files
    mdata.write("my_mudata.h5mu")
    mdata_r = md.read("my_mudata.h5mu", backed=True)
    mdata_r
    return (mdata_r,)


@app.cell
def _(mdata_r):
    # Individual modalities are backed as well inside the .h5mu file
    mdata_r["A"].isbacked
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Multimodal methods

    - For running multimodal methods such as dimensionality reduction, the most simple approach is to concatenate matrices from multiple modalities by horizontally stacking them.
    - This is possible because the number of observations matches across modalities (the number of features per modality does not have to match).
    """)
    return


@app.cell
def _(np):
    # Running a PCA on a concatenated matrix
    def simple_pca(mdata):
        from sklearn import decomposition

        x = np.hstack([m.X for m in mdata.mod.values()])

        pca = decomposition.PCA(n_components=2)
        components = pca.fit_transform(x)

        # By default, methods operate in-place and embeddings are stored in the .obsm slot
        mdata.obsm["X_pca"] = components

    return (simple_pca,)


@app.cell
def _(mdata, simple_pca):
    simple_pca(mdata)
    print(mdata)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.vstack([
        mo.md("""
        # Multimodal data analysis with muon
        - Analogous to Scanpy, general and multimodal functions are grouped in preprocessing (**muon.pp**), tools (**muon.tl**) and plots (**muon.pl**).
        - Unimodal tools are available from the corresponding muon, which are again separated into preprocessing, tools and plots. For example, all ATAC preprocessing functions are grouped into **muon.atac.pp** while CITE-Seq preprocessing functions in **muon.prot.pp**.
        """),
        mo.image("../images/muon_api.png")
    ])

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Demo analysis on 10x Genomics Multiome dataset for human PBMCs
    """)
    return


@app.cell
def _(ln):
    import muon as mu

    mdata_10x = ln.Artifact.get(
        key="introduction/advanced_data_structures_and_frameworks_muon.h5mu", is_latest=True
    ).load()

    fragments = ln.Artifact.get(
        key="introduction/advanced_data_structures_and_frameworks_atac_fragments.tsv.gz",
        is_latest=True,
    ).cache()

    fragments_tbi = ln.Artifact.get(
        key="introduction/advanced_data_structures_and_frameworks_atac_fragments.tsv.gz.tbi",
        is_latest=True,
    ).cache()
    return (mdata_10x,)


@app.cell
def _(mdata_10x):
    # Subsetting to the ATAC modality
    atac = mdata_10x.mod["atac"]
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
