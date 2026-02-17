<div id="main" class="col-md-9" role="main">

# Process data to fit keju model without motif-level shrinkage. Prefers to run on filtered_counts, but will run on raw counts.

<div class="ref-description section level2">

Process data to fit keju model without motif-level shrinkage. Prefers to
run on filtered_counts, but will run on raw counts.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
use_no_motif(keju, G = 50, infer_differential_activity = FALSE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   keju:

    a keju object containing counts

-   G:

    number of enhancers per dispersion parameter

-   infer_differential_activity:

    boolean to use effect size estimation with a control and alternate
    treatment

</div>

<div class="section level2">

## Value

a keju object processed for the no_motif model

</div>

</div>
