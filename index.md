<div id="main" class="col-md-9" role="main">

# Keju

<div class="section level1">

**keju** is an R package for statistical analysis of Massively Parallel
Reporter Assay (MPRA) data, and outputs transcription rate estimates and
differential activity estimates with batch-specific uncertainty
quantification.

<div class="section level2">

## Installation

<div class="section level3">

### R package installation

**keju** relies on [cmdstanr](https://mc-stan.org/cmdstanr/) and
[basilisk](https://www.bioconductor.org/packages/release/bioc/html/basilisk.html).
To install cmdstanr, run

<div id="cb1" class="sourceCode">

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 2)
```

</div>

The compiler requirements can be seen at
[stan-dev](https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms#supported-cpp-versions-and-compilers).
If you run into issues with installation, please ensure your gcc version
is \> 5. Additionally, if you run into issues with the TBB library, try
downgrading CmdStan to version 2.33.1.

**keju** performs some sparse matrix processing steps that require
python. As a result, the R API is a thin wrapper around these python
calls, which require python==3.13.3, numpy==2.2.5, pandas==2.2.3, and
[formulaic](https://matthewwardrop.github.io/formulaic/latest/)==1.1.1.
We use basilisk to install a new python environment with these packages.
See
[basilisk](https://www.bioconductor.org/packages/release/bioc/html/basilisk.html),
specifically `basilisk::configureBasiliskEnv()`, for more information
and for details about installation location. To install basilisk, run

<div id="cb2" class="sourceCode">

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("basilisk")
```

</div>

After installing cmdstanr and basilisk, install **keju** by running

<div id="cb3" class="sourceCode">

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("pimentellab/keju")
library(keju)
```

</div>

</div>

</div>

<div class="section level2">

## Using keju

An introductory vignette can be found
[here](https://github.com/asxue/pimentellab/vignettes/introduction.html).
If you run into problems, please submit and issue on github or email
<albertsxue@gmail.com>.

<div class="section level3">

### Choosing the correct keju model

**keju** is not a singular model, but a suite of models for different
use cases that fit within each other, a little like Russian nesting
dolls. As a result, choosing the correct model for your use case can be
confusing.

-   the simplest modeling option is `no_motif`, and that will get you
    most of the benefits of using keju. If you’re not sure which model
    to use, you should probably at least start with `no_motif`.
-   `motif_shrinkage` is slightly more specialized version of
    `no_motif`. If architectures have some kind of shared structure
    (i.e., multiple architectures test the same transcription factor
    binding motif), `motif_shrinkage` has some statistical niceties that
    may be of interest and slightly better performance given this
    motif-level metadata.
-   `covariate_motif_slope_intercept` is a slightly more specialized
    version of `motif_shrinkage`. If you think your covariates are
    affecting the transcription rate of your architectures (i.e. minimal
    promoter choice (see our paper)), `covariate_motif_slope_intercept`
    can quantify those effects for you given motif-level and
    covariate-level metadata.

If you have questions, the
[vignette](https://github.com/pimentellab/keju/vignettes/introduction.html)
may be helpful.

</div>

</div>

</div>

</div>
