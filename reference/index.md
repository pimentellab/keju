<div id="main" class="col-md-9" role="main">

# Package index

<div class="section level2">

## All functions

</div>

<div class="section level2">

-   `filter()` : Filters architectures without at least an average of
    barcode_threshold barcodes per batch with DNA counts \>
    dna_threshold and RNA counts \> rna_threshold.
-   `keju_fit()` : Fits keju model to data.
-   `keju_from_counts()` : Create keju object from count data. Counts
    and metadata should be matched by barcode and batch.
-   `pretty_summarize()` : Processes parameter estimates using existing
    metadata from the keju object. Calls significance for effect sizes
    and writes to output folder.
-   `use_covariate_slope_intercept()` : Process data to fit keju model
    with motif-level shrinkage and covariate-level effects on
    transcription rate.
-   `use_motif_shrinkage()` : Process counts to fit keju model with
    motif-level shrinkage.
-   `use_no_motif()` : Process data to fit keju model without
    motif-level shrinkage. Prefers to run on filtered_counts, but will
    run on raw counts.

</div>

</div>
