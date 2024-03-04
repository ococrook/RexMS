
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

The goal of `ReX` is to …

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `ReX` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("ReX")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("ReX")
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!

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
