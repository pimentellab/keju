<div id="main" class="col-md-9" role="main">

# Fits keju model to data.

<div class="ref-description section level2">

Fits keju model to data.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
keju_fit(
  keju,
  output_dir,
  model = "no_motif",
  infer_differential_activity = FALSE,
  seed = 1,
  chains = 4,
  parallel_chains = 4
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   keju:

    a filtered and processed keju object

-   output_dir:

    directory for output files

-   model:

    model type. Must be one of 'no_motif', 'motif_shrinkage', or
    'covariate_motif_slope_intercept.

-   infer_differential_activity:

    boolean to use effect size estimation with a control and alternate
    treatment

-   seed:

    random seed for MCMC sampling. See
    https://mc-stan.org/cmdstanr/reference/model-method-sample.html

-   chains:

    number of Markov chains to use. See
    https://mc-stan.org/cmdstanr/reference/model-method-sample.html

-   parallel_chains:

    maximum number of MCMC chains to run in parallel. See
    https://mc-stan.org/cmdstanr/reference/model-method-sample.html

</div>

<div class="section level2">

## Value

a fitted keju object

</div>

</div>
