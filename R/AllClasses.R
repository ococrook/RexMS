##' @slot chains `list()` containing the individual full MCMC chain
##'     results in an `RexChains` instance. Each element must be a
##'     valid `RexChain` instance.
##' @md
##' @rdname RexParams
.RexChains <- setClass("RexChains",
                          slots = c(chains = "list"),
                          validity = function(object) {
                            msg <- validMsg(NULL, NULL)
                            cls <- vapply(object@chains,
                                          function(x) inherits(x, "RexChain"),
                                          logical(1))
                            if (!all(cls))
                              msg <- validMsg(msg, "Not all items are RexChains.")
                            if (is.null(msg)) TRUE
                            else msg
                          })

##' @slot posteriorEstimates A `DataFrame` documenting the posteriors
##'  in an `RexSummary` instance
##' @slot diagnostics A `matrix` of dimensions 1 by 2 containing the
##'     `RexSummary` diagnostics.
##' @slot Rex.joint A `matrix` of dimensions N (samples) by R storing the samples
##'   from the joint posterior distribution of the parameters in an `RexSummary`
##' 
##' @md
##' @rdname RexParams
.RexSummary <- setClass("RexSummary",
                           slots = c(posteriorEstimate = "DFrame",
                                     diagnostics = "matrix",
                                     Rex.joint = "DFrame",
                                     Rex.quantiles = "DFrame",
                                     Rex.globals = "DFrame",
                                     Rex.resolution = "DFrame",
                                     Rex.peptideError = "DFrame"))

##' The `RexParams` infrastructure is used to store and process MCMC results for
##' Rex model from Crook et al 2024.
##' 
##'
##' Objects of the `RexParams` class are created with the `Rex()` function
##' These objects store the *priors* for the model and the results of the MCMC
##' chains, which themselves are stored as an instance of class `RexChains` and
##' can be accessed with the `chains()` function. A summary of the `RexChains`
##' (or class `RexSummary`) can be further computed with the `RexProcess`
##' function.
##' 
##' see the *Rex-MS* vignette for examples
##' @title Infrastructure to to store and process MCMC results
##' 
##' @slot method A `character()` storing the Rex method name
##' @slot priors A `list()` with the priors for the parameters
##' @slot seed An `integer()` with the random number generation seed.
##' @slot summary Object of class `RexSummary` the summarised MCMC results
##' available in the `RexParams` instance.
##' @slot chains Object of class `RexChains` containing the full MCMC results
##' in the `RexParams` instance
##' @slot interval A `numeric()` with the interval for the residues cover
##' @return An object of class `RexParams` which stores the main results
##' for the analysis when using Rex
##' 
##' @md
##' @rdname RexParams
.RexParams <- setClass("RexParams",
                          slots = c(method = "character",
                                    priors = "list",
                                    seed = "integer",
                                    summary = "RexSummaries",
                                    chains = "RexChains",
                                    interval = "numeric"))

##' @title Container for a single Rex chain results
##'
##' @slot dataset `character` indicating the dataset usaully control or treatment
##' @slot replicate `integer` indicating the number of dataset replicate
##' @slot n `integer(1)` indicating the number of MCMC interactions.
##'     Stored in an `RexChain` instance.
##' @slot K `integer(1)` indicating the number of components. Stored
##'     in an `RexChain` instance.
##' @slot N `integer(1)` indicating the number of proteins. Stored in
##'     an `RexChain` instance.
##' @slot niche `matrix(N, n)` component allocation results of an
##'     `RexChain` instance.
##' @slot nicheProb `matrix(N, n, K)` component allocation
##'     probabilities of an `RexChain` instance.
##' @slot outlier `matrix(N, n)` outlier allocation results.
##' @slot outlierProb `matrix(N, n, 2)` outlier allocation
##'     probabilities of an `RexChain` instance.
##' @md
##' @rdname RexParams
.RexChain <- setClass("RexChain",
                         slots = c(dataset = "character",
                                   blong = "matrix",
                                   pilong = "matrix",
                                   qlong = "matrix",
                                   dlong = "matrix",
                                   loglikelihood = "numeric",
                                   UptakeGuess = "matrix",
                                   numBreak = "numeric",
                                   Sigma = "numeric",
                                   R = "integer",
                                   timepoints = "numeric",
                                   numIter = "integer",
                                   numPeptides = "integer",
                                   numTimepoints = "integer"
                                   ),
                         validity = function(object) {
                           msg <- validMsg(NULL, NULL)
                           R <- object@R
                           numIter <- object@numIter
                           numPeptides <- object@numPeptides
                           numTimepoints <- object@numTimepoints
                           
                           if (!identical(nrow(object@blong), R))
                             msg <- validMsg(msg, "Wrong number of residues in blong")
                           if (!identical(nrow(object@pilong), R))
                             msg <- validMsg(msg, "Wrong number of residues in pilong")
                           if (!identical(nrow(object@qlong), R))
                             msg <- validMsg(msg, "Wrong number of residues in qlong")
                           if (!identical(nrow(object@dlong), R))
                             msg <- validMsg(msg, "Wrong number of residues in dlong")
                           if (!identical(nrow(object@UptakeGuess), R))
                             msg <- validMsg(msg, "Wrong number of residues in UptakeGuess")
                           if (!identical(length(object@Sigma), numIter))
                             msg <- validMsg(msg, "Wrong number of iterations in Sigma")
                           if (!identical(length(object@numBreak), numIter))
                             msg <- validMsg(msg, "Wrong number of iterations in numBreak")
                           if (!identical(length(object@loglikelihood), numIter))
                             msg <- validMsg(msg, "Wrong number of iterations in loglikelihood")
                           if (!identical(ncol(object@blong), numIter))
                             msg <- validMsg(msg, "Wrong number of iterations in blong")
                           if (!identical(ncol(object@pilong), numIter))
                             msg <- validMsg(msg, "Wrong number of iterations in pilong")
                           if (!identical(ncol(object@qlong), numIter))
                             msg <- validMsg(msg, "Wrong number of iterations in qlong")
                           if (!identical(ncol(object@dlong), numIter))
                             msg <- validMsg(msg, "Wrong number of iterations in dlong")
                           if (!identical(ncol(object@UptakeGuess), numTimepoints))
                             msg <- validMsg(msg, "Wrong number of timepoints in UptakeGuess")
                           if (is.null(msg)) TRUE
                           else msg
                         })
