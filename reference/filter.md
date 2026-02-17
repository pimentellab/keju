<div id="main" class="col-md-9" role="main">

# Filters architectures without at least an average of barcode_threshold barcodes per batch with DNA counts \> dna_threshold and RNA counts \> rna_threshold.

<div class="ref-description section level2">

Filters architectures without at least an average of barcode_threshold
barcodes per batch with DNA counts \> dna_threshold and RNA counts \>
rna_threshold.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
filter(keju, dna_threshold = 5, rna_threshold = 5, barcode_threshold = 10)
```

</div>

</div>

<div class="section level2">

## Arguments

-   keju:

    a keju object containing counts

-   dna_threshold:

    minimum DNA counts per barcode

-   rna_threshold:

    minimum RNA counts per barcode

-   barcode_threshold:

    minimum average barcodes per batch

</div>

<div class="section level2">

## Value

a keju object with the filtered_counts object

</div>

</div>
