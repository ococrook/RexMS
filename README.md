
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ReX

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/ReX.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/ReX)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/ReX.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ReX)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/ReX.svg)](http://bioconductor.org/packages/stats/bioc/ReX/)
[![Bioc
support](https://bioconductor.org/shields/posts/ReX.svg)](https://support.bioconductor.org/tag/ReX)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/ReX.svg)](https://bioconductor.org/packages/release/bioc/html/ReX.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/ReX.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/ReX/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/ReX.svg)](https://bioconductor.org/packages/release/bioc/html/ReX.html#since)
[![R-CMD-check-bioc](https://github.com/ococrook/ReX/actions/workflows/R-CMD-check-bioc.yaml/badge.svg)](https://github.com/ococrook/ReX/actions/workflows/R-CMD-check-bioc.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ococrook/ReX/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/ococrook/ReX?branch=devel)
<!-- badges: end -->

The goal of `ReX` is to analyse hydrogen deuterium exchange mass
spectrometry (HDX-MS) data at the residue level. The package is designed
to be used with data from various platforms that have already been
process at Spectrum-level. `ReX` takes the processed peptide-level data
and infers the residue level deuterium uptake. The underlying model is a
Bayesian non-parametric model that recasts HDX-MS data analysis as a
(latent) change-point detection problem. The unique benefits of this
model are the following: (1) It can provide statistical confidence at
the level of residues (e.g. a probability of change) (2) Borrow
statistical power from overlapping peptides (3) Infer uptakes patterns
that are hidden at the peptide-level because of averaging (4) Provide
global and pre-residue resolution metrics (5) It can perform single
protein analysis, differential analysis and confrontational signature
analysis (many states/compounds) (6) You can build predictive models
with the inferred uptakes using partial least squares discriminant
analysis (PLS-DA). (7) You can costumize the model to your specific
needs.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `ReX` from
[Bioconductor](http://bioconductor.org/) using the following code:

Note that you cannot currently install from Bioconductor

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("ReX")
```

To install the development version of `ReX`, it’s easiest to use the
`remotes` package.

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("ococrook/ReX")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("ReX")
## basic example code
```

## Citation

Below is the citation output from using `citation('ReX')` in R. Please
run this yourself to check for any updates on how to cite **ReX**.

``` r
print(citation("ReX"), bibtex = TRUE)
#> To cite package 'ReX' in publications use:
#> 
#>   ococrook (2024). _Inferring residue level hydrogen deuterium exchange
#>   with ReX_. doi:10.18129/B9.bioc.ReX
#>   <https://doi.org/10.18129/B9.bioc.ReX>,
#>   https://github.com/ococrook/ReX/ReX - R package version 0.99.0,
#>   <http://www.bioconductor.org/packages/ReX>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Inferring residue level hydrogen deuterium exchange with ReX},
#>     author = {{ococrook}},
#>     year = {2024},
#>     url = {http://www.bioconductor.org/packages/ReX},
#>     note = {https://github.com/ococrook/ReX/ReX - R package version 0.99.0},
#>     doi = {10.18129/B9.bioc.ReX},
#>   }
#> 
#>   ococrook (2024). "Inferring residue level hydrogen deuterium exchange
#>   with ReX." _bioRxiv_. doi:10.1101/TODO
#>   <https://doi.org/10.1101/TODO>,
#>   <https://www.biorxiv.org/content/10.1101/TODO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Inferring residue level hydrogen deuterium exchange with ReX},
#>     author = {{ococrook}},
#>     year = {2024},
#>     journal = {bioRxiv},
#>     doi = {10.1101/TODO},
#>     url = {https://www.biorxiv.org/content/10.1101/TODO},
#>   }
```

Please note that the `ReX` was only made possible thanks to many other R
and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `ReX` project is released with a [Contributor Code
of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
