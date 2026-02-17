<div id="main" class="col-md-9" role="main">

# Processes parameter estimates using existing metadata from the keju object. Calls significance for effect sizes and writes to output folder.

<div class="ref-description section level2">

Processes parameter estimates using existing metadata from the keju
object. Calls significance for effect sizes and writes to output folder.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
pretty_summarize(keju, model, output_folder, infer_differential_activity)
```

</div>

</div>

<div class="section level2">

## Arguments

-   keju:

    a fitted keju object

-   model:

    the type of model fit by keju

-   output_folder:

    folder for output estimate csv files

-   infer_differential_activity:

    boolean to use effect size estimation with a control and alternate
    treatment

</div>

</div>
