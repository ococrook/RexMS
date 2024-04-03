
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RexMS

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/RexMS.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/RexMS)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/RexMS.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/RexMS)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/RexMS.svg)](http://bioconductor.org/packages/stats/bioc/RexMS/)
[![Bioc
support](https://bioconductor.org/shields/posts/RexMS.svg)](https://support.bioconductor.org/tag/RexMS)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/RexMS.svg)](https://bioconductor.org/packages/release/bioc/html/RexMS.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/RexMS.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/RexMS/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/RexMS.svg)](https://bioconductor.org/packages/release/bioc/html/RexMS.html#since)
[![R-CMD-check-bioc](https://github.com/ococrook/RexMS/actions/workflows/R-CMD-check-bioc.yaml/badge.svg)](https://github.com/ococrook/RexMS/actions/workflows/R-CMD-check-bioc.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ococrook/RexMS/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/ococrook/RexMS?branch=devel)
<!-- badges: end -->

The goal of `RexMS` is to analyse hydrogen deuterium exchange mass
spectrometry (HDX-MS) data at the residue level. The package is designed
to be used with data from various platforms that have already been
process at Spectrum-level. `RexMS` takes the processed peptide-level
data and infers the residue level deuterium uptake. The underlying model
is a Bayesian non-parametric model that recasts HDX-MS data analysis as
a (latent) change-point detection problem. The unique benefits of this
model are the following:

1)  It can provide statistical confidence at the level of residues (e.g.
    a probability of change)

2)  Borrow statistical power from overlapping peptides

3)  Infer uptake patterns that are hidden at the peptide-level because
    of averaging

4)  Provide global and per-residue resolution metrics

5)  It can perform single protein analysis, differential analysis and
    conformational signature analysis (many states/compounds), linking
    results to downstream functional outcomes

6)  You can build predictive models with the inferred uptakes using
    partial least squares discriminant analysis (PLS-DA).

7)  You can costumize the model to your specific needs.

## Documentation

The official documentation is available at
[ococrook.github.io/RexMS](https://ococrook.github.io/RexMS/). We
encourage use to look there and in particular the vignettes. We suggest
reading through all of the vignettes first so you understand which part
of `RexMS` you need to use.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `RexMS` from
[Bioconductor](http://bioconductor.org/) using the following code:

Note that you cannot currently install from Bioconductor

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("RexMS")
```

To install the development version of `RexMS`, it’s easiest to use the
`remotes` package.

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("ococrook/RexMS")
```

## Example

RexMS has structure visualisation at part of its core functionality.
Below is an example of how to use the `hdx_to_pdb_colours` function to
map HDX-MS data onto a PDB file.

``` r
library("RexMS")
library(NGLVieweR)

# generate random HDX data
v <- matrix(rnorm(n = 477), nrow = 1)
colnames(v) <- seq.int(ncol(v)) # residue numbering

v2 <- v[, seq.int(344, 477), drop = FALSE]
colnames(v2) <- seq.int(ncol(v2))

pdb_filepath <- system.file("extdata", "test_BRD4.pdb", mustWork = TRUE, package = "RexMS")

# generate a protection-deprotection colour mapping
mycolor_parameters <- hdx_to_pdb_colours(v2, pdb = pdb_filepath, cmap_name = "ProtDeprot")

# Note this will open in a view
view_structure(pdb_filepath = pdb_filepath, color_parameters = mycolor_parameters)
```

The expected output is a 3D structure of the protein with the HDX-MS
data mapped onto it. The colouring is based on the
protection-deprotection scale.The following image is indicative and may
not be the same as the output you will see.

<figure>
<img src="man/figures/pdb_hdx.png"
alt="An example figure generated from the RexMS package" />
<figcaption aria-hidden="true">An example figure generated from the
RexMS package</figcaption>
</figure>

## Citation

Below is the citation output from using `citation('RexMS')` in R. Please
run this yourself to check for any updates on how to cite **RexMS**.

``` r
print(citation("RexMS"), bibtex = TRUE)
#> To cite package 'RexMS' in publications use:
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

Please note that the `RexMS` was only made possible thanks to many other
R and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `RexMS` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.18/BiocCheck)*.
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
*[biocthis](https://bioconductor.org/packages/3.18/biocthis)*.
