<div id="main" class="col-md-9" role="main">

# Create keju object from count data. Counts and metadata should be matched by barcode and batch.

<div class="ref-description section level2">

Create keju object from count data. Counts and metadata should be
matched by barcode and batch.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
keju_from_counts(
  barcode,
  R,
  D,
  batch,
  dna_batch,
  architecture,
  treatment,
  is_control_treatment = FALSE,
  is_control_architecture = FALSE,
  covariates = "correction",
  motif = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   barcode:

    a length N vector of barcodes

-   R:

    a length N vector of RNA counts

-   D:

    a length N vector of DNA counts

-   batch:

    a length N vector denoting RNA batches

-   dna_batch:

    a length N vector denoting DNA batches

-   architecture:

    a length N vector of enhancer architectures or names. Transcription
    rate and effect size estimates will be fit on these.

-   treatment:

    a length N vector detailing the treatments.

-   is_control_treatment:

    a length N boolean vector denoting the control treatment.

-   is_control_architecture:

    a length N boolean vector denoting architectures that are negative
    controls.

-   covariates:

    a length N vector denoting covariates to be used in effect size
    correction and potentially slope/intercept fitting. If none, will
    default to a vector of ones of length N, to perform global effect
    size correction.

-   motif:

    a length N vector denoting motifs shared across architectures. Used
    to shrink architecture-level estimates to a shared motif-level
    estimate. Necessary for motif_shrinkage and
    covariate_motif_slope_intercept models.

</div>

<div class="section level2">

## Value

an initialized keju object with a counts object

</div>

</div>
