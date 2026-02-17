<div id="main" class="col-md-9" role="main">

# Process counts to fit keju model with motif-level shrinkage.

<div class="ref-description section level2">

Process counts to fit keju model with motif-level shrinkage.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
use_motif_shrinkage(keju, G = 50, infer_differential_activity = FALSE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   keju:

    a keju object containing counts or filtered_counts object.

-   G:

    number of enhancers per dispersion parameter

-   infer_differential_activity:

    boolean to use effect size estimation with a control and alternate
    treatment

</div>

<div class="section level2">

## Value

a keju object processed for the motif_shrinkage model

</div>

</div>
