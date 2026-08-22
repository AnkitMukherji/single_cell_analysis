import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.vstack([
        mo.md("# Transcript Quantification\n\n### Transcription takes place in stochastic bursts and leads to a total of 1,00,000 to 10,00,000 mRNA molecules in a typical mammalian cell"),
        mo.image(src="images/quantifying_gene_expression.png", width=600)
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    There are two major approaches:
    - **Full-length**
        - Covers the whole transcript
        - Can be captured by *plate-based* protocols only
        - Helps in the detection of splice variants
    - **Tag-based**
        - Captures only the 5' or 3' ends of transcripts
        - Allows for the usage of *UMIs (Unique Molecular Identifiers)* which helps in resolving biases in the transcript amplification process

    # Single-cell sequencing protocols

    There are two major protocols:
    - **Separation in droplets**
        - High-throughput
        - Examples include *inDrop*, *Drop-seq*, *10X Genomics Chromium*
        - Droplet contains Cell + Bead
        - There are on-bead primers comprising of:
            - PCR handle
            - Cell barcode
            - 4-8 bp long UMI
            - poly-T tail
        - Cell lysis inside the droplet ----> mRNA captured by bead ----> Droplets broken ----> Release STAMPs (single-cell transcriptomes attached to microparticles) ----> PCR and Reverse transcription ----> Tagmentation (adaptor attachment) ----> Sequencing
        - *<u>Proportion of reads originating from valid barcodes</u>*
            - 10X Genomics -> *75%*
            - InDrop -> *25%*
            - Drop-seq -> *30%*
        - *<u>Number of transcripts captured</u>*
            - 10X Genomics -> *17000 transcripts from 3000 genes*
            - InDrop -> *8000 transcripts from 2500 genes*
            - Drop-seq -> *2700 transcripts from 1250 genes*
        - *<u>Based on the numbers above, in terms of technical noise:</u>*
              <br>**InDrop >> Drop-seq >> 10X Genomics**</br>
        - *<u>Biasedness</u>*
            - 10X Genomics
                - Amplification of shorter genes
                - Genes with high GC content
            - Drop-seq
                - Genes with lower GC content
    - **Separation in physical compartments (wells)**
        - Examples include *Plate based protocols*, *Fluidigm C1*
        - *Plate-based protocols*
            - Cell-sorting by FACS (sorted according to specific-cell surface markers / micropipetting) ----> placed into individual wells containing cell lysis buffers ----> Reverse transcription
            - Approximately 5k - 10k genes captured per cell
            - Examples include SMART-seq2, MARS-seq, QUARTZ-seq, SRCB-seq
            - QUARTZ-seq capture more genes than the other three
        - *Fluidigm C1*
            - Examples include CEL-seq2, SMART-seq1
            - Uses microfluidic chips that loads and separates cells into small reaction chambers


    # Selection of protocol
    - *Deep-characterization of a specific cell type population* -> **Plate-based methods**
    - *Capturing an overall heterogenous mixture of cells* -> **Droplet-based assays**

    # Single-cell vs Single-nuclei sequencing
    Single-cell sequencing  | Single-nuclei sequencing
    ------------- | -------------
    Requires fresh tissue samples | Can be captured from frozen tissues
    Sensitive to tissue dissociating enzymes | Resistant to mechanical force
    """)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
