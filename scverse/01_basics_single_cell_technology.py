import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium", auto_download=["html"])


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # ----------------Basics of Single-cell RNA sequencing----------------
    """)
    return


@app.cell
def _(mo):
    mo.vstack([
        # Section 1: Transcript Quantification
        mo.md("# Transcript Quantification\n\n### Transcription takes place in stochastic bursts and leads to a total of 1,00,000 to 10,00,000 mRNA molecules in a typical mammalian cell"),
        mo.image(src="images/quantifying_gene_expression.png", width=600),

        mo.md("---"),

        # Section 2: Non-coding RNA
        mo.md("# [Overview of non-coding RNA](https://www.bio-rad.com/de-de/applications-technologies/coding-non-coding-rna?ID=Q1070M70KWE7)"),
        mo.image(src="images/ncRNA_Overview-Fig1a.jpg", width=600),
        mo.image(src="images/ncRNA_Overview-Fig1b.jpg", width=600)
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Why we need negative binomial to model gene transcription?
    - Gene transcription occurs in stochastic bursts - short, irregular periods of activity during which a gene might suddenly produce multiple mRNA transcripts before returning to silence.
    - Counting events over a fixed interval is a Poisson distribution, but a Poisson has a single parameter, which forces its variance to equal its mean.
    - Negative binomial is ideal as it models event counts (mRNA) while capturing overdispersion (variance exceeding the mean) caused by transcriptional bursts.

    # There are two major approaches for transcript quantification:
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


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # -----------------------Raw Data Processing-----------------------
    """)
    return


@app.cell
def _(mo):
    mo.vstack([
        mo.md("# Overview of raw data processing"),
        mo.image(src="images/overview_raw_data_processing.jpg")
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Errors that can happen during lane demultiplexing
    ## Index Hopping
    Reads from a RNA-seq library are assigned to a wrong sample.
    - Caused by a partial library fragment annealing to a fragment in another nanowell (in flow cell), creating a chimeric molecule with a swapped sample index.
    <!---->
    - May result from unligated adaptors in solution, which can hybridize to the complementary sequence of the adaptor on the unbound end of a properly hybridized fragment and extend, adding the incorrect sample index to that sequence and creating a hybrid library in that nanowell. If the index-hopped strand then denatures and seeds another nanowell, it will generate a false signal, incorrectly assigning that read to a different sample during demultiplexing.
    <!---->
    - Long stretches of low-diversity in cDNA sequence that allow one fragment to bind to another in a separate nanowell.
    <!---->
    - DNA breakage during amplification.

    Ways to deal with index hopping:
    - Cleaning unligated adaptors from library.
    - Dual indexing eliminates sample misalignment due to sequencing errors.
    """)
    return


@app.cell
def _(mo):
    mo.vstack([
        mo.md("# A representative Chromium Single Cell 3’ dual index library schematic with sample index i5 and i7."),
        mo.image(src="images/dual_indexing.png"),
        mo.md("----"),
        mo.md("# Removal of index-hopped reads during demultiplexing."),
        mo.image(src="images/lane_demux_dual_index.png")
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Important points to remember while QC:
    - **Per base sequence content:** For single-cell data, it is common to observe fluctuations at the start of the reads. This occurs because the initial bases represent the sequence of the priming sites, which is often not perfectly random.
    - **Per base GC content:**
        - GC-content of transcriptome may appear wider or narrower than the genome's expected GC distribution.
        - Interpreting GC content in transcriptomes can be challenging because the expected GC distribution depends not only on the sequence composition of the transcriptome but also on gene expression levels in the sample (unknown beforehand).
    - **Sequence length distribution:** Trimming performed for quality assessment leads to variation in read lengths.
    - **Sequence duplication levels:**
        - Although there is a step of deduplication, FastQC is not UMI-aware. Single-cell platforms perform multiple rounds of PCR therefore, highly expressed genes naturally produce large number of transcripts.
        - It is common for small subset of sequence to show high duplication levels.
    - **Overrepresented sequence:**
        - Identifies read sequence that constitute more than 0.1% of the total reads.
        - In scRNA-seq, some overrepresented sequence may arise from highly expressed genes amplified during PCR. However, the majority of sequence should not be overrepresented.

    # Alignment and Mapping

    Raw sequence files include:
    - Cell Barcode (CB): Unique identifiers for individual cells.
    - Unique Molecular Identifiers (UMIs): Tags that distinguish individual molecules to account for amplification bias.
    - Raw cDNA sequences: The actual read sequence generated from the molecules.

    Alignment approaches can be categorized along two axes:
    1. **The type of mapping they perform**
        - Alignment-based methods (Uses Dynamic Programming algorithm)
            - Spliced alignment
            - Contiguous alignment
        - Lightweight mapping-based approaches [*Note: Do not yet support spliced mapping of reads against a reference genome*]
    3. **The type of reference sequence against which they map reads**
        - Full reference genome (annotated)
            - Tools using this approach are zUMIs, CellRanger and STARsolo. All are splice-aware aligners.
            - Suitable for both single-cell and single-nuclei data
        - Annotated transcriptome
            - Length of transcriptome << Length of Genome therefore, reducing the computational load of mapping.
            - Used in contiguous alignment and lightweight mapping.
            - Fails to capture reads arising outside the spliced transcriptome. Therefore, not suitable for processing single-nuclei data.
        - Augmented transcriptome
            - Spliced transcript sequences are augmented with additional reference sequence such as full-length unspliced transcripts or excised intronic sequence.
            - Suitable for both single-cell and single-nuclei data and methods like RNA velocity.

    # Cell Barcode Correction

    In droplet-based methods, RNA of each captured cell is extracted within a water-based droplet encapsulation along with a barcoded bead. The barcoded bead tags the RNA with **unique cell barcodes (CBs)**.

    ## Errors in barcoding

    Observed cell barcodes can differ significantly from the number of originally encapsulated cells. Key sources of error are:
    - Doublets/Multiplets: Single barcode -> Multiple cells
    - Empty droplets with no encapsulated cells
    - Sequencing errors

    ## Strategies for cell barcode identification and correction
    1. **Correction against a known list of potential barcodes (called whitelist)**
        - Any barcode not in list is corrected by finding its closest match.
        - Limitations: Corrupted barcode closely resembles multiple barcodes in the permit list therefore, its correction becomes ambiguous. It can be resolved by:
            - Corrected barcode should appear in sample in a high frequency threshold.
            - Base quality of the corrected position.
            - Complete removal of reads with ambiguous barcodes.
    2. **Knee or elbow-based methods**
        - Assumes high quality barcodes are those associated with the highest number of reads in the sample.
        - Cumulative frequency graph can be plot where barcodes are sorted in descending order based on the number of distinct reads or UMIs they are associated with. It will contain a “knee” or “elbow” – an inflection point that can be used to characterize frequently occurring barcodes from infrequent (and therefore likely erroneous) barcodes.
        - The set of barcodes that appear “above” the knee can be treated as a permit list against which the rest of the barcodes may be corrected, as in the first method list above.
    3. **Filtering and correction based on an expected cell count**
        - When barcode frequency distributions lack a clear knee or show bimodal patterns due to technical artifacts, barcode correction can be guided by a user-provided expected cell count.
        - The cell barcodes having less frequency are corrected against the valid list [barcodes with the expected cell count] by attempting to correct uniquely to one of the valid barcodes based on sequence similarity.
    4. **Filtering based on a forced number of valid cells**
        - The simplest approach, although potentially problematic, is for the user to manually specify the number of valid barcodes.

    # UMI Resolution
    """)
    return


if __name__ == "__main__":
    app.run()
