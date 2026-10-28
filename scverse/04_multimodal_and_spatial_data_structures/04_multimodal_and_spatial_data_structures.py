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

    # MuData object with n_obs × n_vars = 11909 × 144978
    #   var:	'gene_ids', 'feature_types', 'genome', 'interval'
    #   2 modalities
    #     rna:	11909 x 36601
    #       var:	'gene_ids', 'feature_types', 'genome', 'interval'
    #     atac:	11909 x 108377
    #       var:	'gene_ids', 'feature_types', 'genome', 'interval'
    #       uns:	'atac', 'files'
    return fragments, mdata_10x, mu


@app.cell
def _(mdata_10x):
    # Subsetting to the ATAC modality
    atac = mdata_10x.mod["atac"]
    return (atac,)


@app.cell
def _(atac):
    # Filtering out cells with too few peaks and peaks detected in too few cells
    import scanpy as sc

    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
    return (sc,)


@app.cell
def _(atac):
    atac
    return


@app.cell
def _(atac, sc):
    sc.pl.violin(atac, ["total_counts", "n_genes_by_counts"], jitter=0.4, multi_panel=True)
    return


@app.cell
def _(atac, mu):
    # Filter peaks whose expression is not detected
    mu.pp.filter_var(atac, "n_cells_by_counts", lambda x: x >= 10)
    # This is analogous to sc.pp.filter_genes(rna, min_cells=10)
    # but does in-place filtering and avoids copying the object
    return


@app.cell
def _(mdata_10x):
    mdata_10x.mod["atac"]
    return


@app.cell
def _(atac):
    atac
    return


@app.cell
def _(atac, mu):
    # Filtering cells
    mu.pp.filter_obs(atac, "n_genes_by_counts", lambda x: (x >= 2000) & (x <= 15000))
    # This is analogous to
    # sc.pp.filter_cells(atac, max_genes=15000)
    # sc.pp.filter_cells(atac, min_genes=2000)
    # but does in-place filtering avoiding copying the object

    mu.pp.filter_obs(atac, "total_counts", lambda x: (x >= 4000) & (x <= 40000))
    return


@app.cell
def _(atac, sc):
    sc.pl.violin(atac, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True)
    return


@app.cell
def _(atac, mu):
    # Same plots using Muon
    mu.pl.histogram(atac, ["n_genes_by_counts", "total_counts"])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## ATAC specific quality control with muon
    """)
    return


@app.cell
def _():
    from muon import atac as ac

    return (ac,)


@app.cell
def _(ac, atac, fragments, mdata_10x):
    atac.obs["NS"] = 1
    ac.tl.locate_fragments(mdata_10x, fragments=str(fragments))
    ac.tl.locate_fragments(atac, fragments=str(fragments))
    ac.tl.nucleosome_signal(atac, n=1e6)
    ac.tl.get_gene_annotation_from_rna(mdata_10x["rna"]).head(3)
    tss = ac.tl.tss_enrichment(mdata_10x, n_tss=1000)
    ac.pl.tss_enrichment(tss)
    return


@app.cell
def _(atac, sc):
    # Save original counts, normalize data and select highly variable genes with scanpy
    atac.layers["counts"] = atac.X
    sc.pp.normalize_per_cell(atac, counts_per_cell_after=1e4)
    sc.pp.log1p(atac)
    sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=0.5)
    atac.raw = atac
    return


@app.cell
def _(atac):
    atac
    return


@app.cell
def _(ac, atac, sc):
    # Athough PCA is used for ATAC data, latent semantic indexing (LSI) is another option
    ac.tl.lsi(atac)
    sc.pp.neighbors(atac, use_rep="X_lsi", n_neighbors=10, n_pcs=30)
    sc.tl.leiden(atac, resolution=0.5)
    sc.tl.umap(atac, spread=1.5, min_dist=0.5, random_state=20)
    sc.pl.umap(atac, color=["leiden", "n_genes_by_counts"], legend_loc="on data")
    return


@app.cell
def _(ac, atac):
    # Gene specific UMAP
    ac.pl.umap(atac, color=["KLF4"], average="peak_type")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Using SpatialData to store multimodal and spatial data
    - Spatial omics measure molecular information (gene expression, chromatin accessibility, proteins) while preserving their physical locations inside a tissue.
    - Single-cell methods lose spatial context during tissue dissociation, whereas spatial omics retains both molecular profiles and spatial positions of individual cells.
    - The datasets combine images, coordinates and multiple molecular measurements making them large, complex and requiring specialized spatially aware data structures.

    ## SpatialData object
    - SpatialData is a container for various **Elements**. An **Element** is either a **SpatialElement** (Images, Labels, Points, Shapes) or a Table.
        - **Images:** H&E, staining images
        - **Labels:** pixel-level segmentation
        - **Points:** transcripts locations with gene information, landmark points
        - **Shapes:** cell/nucleus boundaries, subcellular structures, anatomical annotations, regions of interest (ROIs)
        - **Tables:** sparse/dense matrices annotating the SpatialElements or storing arbitrary (non-spatial) metadata. They do not contain spatial coordinates.
    - **SpatialElements** can be segregated into two broad types:
        - **Rasters:** Data made up of pixels: including **Images** and **Labels**.
        - **Vectors:** Data made up of points and lines. Polygons are also vectors, since they are simply a list of connected points. **Points** and **Shapes** are elements of this type.
    """)
    return


@app.cell
def _(mo):
    mo.vstack([
        mo.md("## Elements of a SpatialData object"),
        mo.image("../images/SpatialData_elements.png")
    ])
    return


@app.cell
def _():
    import spatialdata as sd
    import spatialdata_plot  # noqa: F401

    return (sd,)


@app.cell
def _(ln):
    sdata = ln.Artifact.get(
        key="introduction/advanced_data_structures_and_frameworks_spatialdata.zarr"
    ).load()
    sdata
    return (sdata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## **Elements** of the **SpatialData** object and their representation
    - **Images:** **xarray.DataArray** or **xarray.DataTree** objects respectively for single-scale or multi-scale images.
    - **Labels:** **xarray.DataArray** or **xarray.DataTree** objects containing integer codes for different labels.
    - **Points:** **dask.DataFrame** objects containing point coordinates, (lazy version of pandas.DataFrame objects).
    - **Shapes:** **geopandas.GeoDataFrame** objects containing geometric objects like polygons.
    - **Tables:** **anndata.AnnData** objects for tabular data with annotations.

    Images can have multiple scales, which are stored in a pyramid for downsampled representations useful for more efficient visualization and analysis. We can access a single scale using the **get_pyramid_levels** function, where 0 is the highest resolution, 1 is downsampled by a factor of (usually 2), 2 is downsampled by another factor (usually 2, with respect to the previous scale).
    """)
    return


@app.cell
def _(sdata):
    # The image has two scales
    sdata["raw_image"]
    return


@app.cell
def _(sd, sdata):
    # Accessing the highest resolution image
    sd.get_pyramid_levels(sdata["raw_image"], n=1)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - To view image properties such as axes and channel names, **SpatialData** provides some helper functions.
    - The axes **x** and **y** represent the two-dimensional plane of the image, while **c** corresponds to the different channels (e.g., wavelength) within the image.
    """)
    return


@app.cell
def _(sd, sdata):
    sd.models.get_axes_names(sdata["raw_image"])
    return


@app.cell
def _(sd, sdata):
    # The image has one channel for DAPI staining (highlights nucleus)
    sd.models.get_channel_names(sdata["raw_image"])
    return


@app.cell
def _(sdata):
    # Using spatialdata-plot to show the image (DAPI staining for nuclei)
    sdata.pl.render_images("raw_image", cmap="gray").pl.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
