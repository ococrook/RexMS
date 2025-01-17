% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rex-conformationalSignatures.R
\name{UnsupervisedCSA}
\alias{UnsupervisedCSA}
\title{Funcition for performing unsupervised conformational signature analysis.
This performs PCA on the TRE values of the RexDifferentialList and
returns the PCA object and the TRE values in wide format.}
\usage{
UnsupervisedCSA(
  RexDifferentialList,
  quantity = "TRE",
  states,
  whichTimepoint = 600,
  pca_params = list(scale = FALSE, center = TRUE)
)
}
\arguments{
\item{RexDifferentialList}{A list of RexDifferential objects}

\item{quantity}{The quantity to use for the analysis. Default is "TRE"}

\item{states}{The state name to use for the analysis/ e.g. ligand used in
differential analysis}

\item{whichTimepoint}{The timepoint to use for the analysis. Default is 600}

\item{pca_params}{The parameters to use for the PCA.
Default is list(scale = FALSE, center = TRUE)}
}
\value{
A list containing the PCA object and the TRE values in wide format
}
\description{
Funcition for performing unsupervised conformational signature analysis.
This performs PCA on the TRE values of the RexDifferentialList and
returns the PCA object and the TRE values in wide format.
}
\examples{
library("RexMS")
data("out_lxr_compound_proccessed")
data("LXRalpha_compounds")

states <- names(LXRalpha_compounds)


ucsa <- UnsupervisedCSA(out_lxr_compound_proccessed,
                       quantity = "TRE",
                       states = states,
                       whichTimepoint = 600,
                       pca_params = list(scale = FALSE,
                                         center = TRUE))
                       


}
