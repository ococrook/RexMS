##' @slot chains `list()` containing the individual full MCMC chain
##'     results in an `RexChains` instance. Each element must be a
##'     valid `RexChain` instance.
##' @md
##' @rdname RexParams
.RexChains <- setClass("RexChains",
    slots = c(chains = "list"),
    validity = function(object) {
        msg <- validMsg(NULL, NULL)
        cls <- vapply(
            object@chains,
            function(x) inherits(x, "RexChain"),
            logical(1)
        )
        if (!all(cls)) {
            msg <- validMsg(msg, "Not all items are RexChains.")
        }
        if (is.null(msg)) {
            TRUE
        } else {
            msg
        }
    }
)

##' @slot posteriorEstimates A `DataFrame` documenting the posteriors
##'  in an `RexSummary` instance
##' @slot diagnostics A `matrix` of dimensions 1 by 2 containing the
##'     `RexSummary` diagnostics.
##' @slot Rex.quantiles A `DataFrame` of with R rows indicated quantiles for key
##' parameters in the Rex model
##' @slot Rex.globals A `DataFrame` of with R rows indicated global parameters
##' at the moment the standard deviation sigma and the model likelihoods
##' @slot Rex.resolution A `DataFrame` of with R rows indicated resolution metrics
##' @slot Rex.peptideError A `DataFrame` indicated peptide error
##'
##' @md
##' @rdname RexParams
.RexSummary <- setClass("RexSummary",
    slots = c(
        posteriorEstimates = "DFrame",
        diagnostics = "matrix",
        Rex.quantiles = "DFrame",
        Rex.globals = "DFrame",
        Rex.resolution = "DFrame",
        Rex.peptideError = "DFrame"
    )
)

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
    slots = c(
        method = "character",
        priors = "list",
        seed = "integer",
        summary = "RexSummary",
        chains = "RexChains",
        interval = "numeric"
    )
)

##' @title Container for a single Rex chain results
##'
##' @slot dataset `character` indicating the dataset
##' @slot blong `matrix` of dimensions R by numIter containing the
##'    blong values for the residues
##' @slot pilong `matrix` of dimensions R by numIter containing the
##'   pilong values for the residues
##' @slot qlong `matrix` of dimensions R by numIter containing the
##'  qlong values for the residues
##' @slot dlong `matrix` of dimensions R by numIter containing the
##'  dlong values for the residues
##' @slot loglikelihood `numeric` of length numIter containing the
##'  loglikelihood values
##' @slot phi `numeric` of length 1 containing the phi values for the model
##' @slot UptakeGuess `list` each with a matrix of dimensions R by numTimepoints
##'  containing the UptakeGuess values for the residues
##' @slot numBreak `numeric` of length numIter containing the number of breakpoints
##' @slot Sigma `numeric` of length numIter containing the Sigma values
##' @slot R `integer` of length 1 containing the number of residues
##' @slot timepoints `numeric` of length numTimepoints containing the timepoints
##' @slot numIter `integer` of length 1 containing the number of mcmc iterations
##' @slot numPeptides `integer` of length 1 containing the number of peptides
##' @slot numTimepoints `integer` of length 1 containing the number of timepoints
##'
##' @md
##' @rdname RexParams
.RexChain <- setClass("RexChain",
    slots = c(
        dataset = "character",
        blong = "matrix",
        pilong = "matrix",
        qlong = "matrix",
        dlong = "matrix",
        loglikelihood = "numeric",
        phi = "numeric",
        UptakeGuess = "list",
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

        if (!identical(nrow(object@blong), R)) {
            msg <- validMsg(msg, "Wrong number of residues in blong")
        }
        if (!identical(nrow(object@pilong), R)) {
            msg <- validMsg(msg, "Wrong number of residues in pilong")
        }
        if (!identical(nrow(object@qlong), R)) {
            msg <- validMsg(msg, "Wrong number of residues in qlong")
        }
        if (!identical(nrow(object@dlong), R)) {
            msg <- validMsg(msg, "Wrong number of residues in dlong")
        }
        if (!identical(nrow(object@UptakeGuess[[1]]), R)) {
            msg <- validMsg(msg, "Wrong number of residues in UptakeGuess")
        }
        if (!identical(length(object@Sigma), numIter)) {
            msg <- validMsg(msg, "Wrong number of iterations in Sigma")
        }
        if (!identical(length(object@numBreak), numIter)) {
            msg <- validMsg(msg, "Wrong number of iterations in numBreak")
        }
        if (!identical(length(object@loglikelihood), numIter)) {
            msg <- validMsg(msg, "Wrong number of iterations in loglikelihood")
        }
        if (!identical(ncol(object@blong), numIter)) {
            msg <- validMsg(msg, "Wrong number of iterations in blong")
        }
        if (!identical(ncol(object@pilong), numIter)) {
            msg <- validMsg(msg, "Wrong number of iterations in pilong")
        }
        if (!identical(ncol(object@qlong), numIter)) {
            msg <- validMsg(msg, "Wrong number of iterations in qlong")
        }
        if (!identical(ncol(object@dlong), numIter)) {
            msg <- validMsg(msg, "Wrong number of iterations in dlong")
        }
        if (!identical(ncol(object@UptakeGuess[[1]]), numTimepoints)) {
            msg <- validMsg(msg, "Wrong number of timepoints in UptakeGuess")
        }
        if (is.null(msg)) {
            TRUE
        } else {
            msg
        }
    }
)
##' Container for rex differential results
##'
##' @slot dataset `character` indicating the dataset
##' @slot Rex.predictionError `DFrame` containing the prediction error
##' @slot Rex.estimates `DFrame` containing the estimates ARE, TRE, signed ARE
##' @slot Rex.probs `DFrame` containing the probabilities and total probabilities
##' @slot Rex.eFDR `DFrame` containing the eFDR
##'
##' @md
##' @rdname RexParams
.RexDifferential <- setClass("RexDifferential",
    slots = c(
        dataset = "character",
        Rex.predictionError = "DFrame",
        Rex.estimates = "DFrame",
        Rex.probs = "DFrame",
        Rex.eFDR = "DFrame"
    ),
    validity = function(object) {
        msg <- validMsg(NULL, NULL)
        R <- nrow(object@Rex.estimates)
        numPeptides <- nrow(object@Rex.predictionError)
        numTimepoints <- ncol(object@Rex.predictionError) - 1

        if (!identical(nrow(object@Rex.estimates), R)) {
            msg <- validMsg(msg, "Wrong number of residues in estimates")
        }
        if (!identical(nrow(object@Rex.probs), R)) {
            msg <- validMsg(msg, "Wrong number of residues in probs")
        }
        if (!identical(ncol(object@Rex.eFDR), as.integer(numTimepoints))) {
            msg <- validMsg(msg, "Wrong number of columns in eFDR")
        }
        if (is.null(msg)) {
            TRUE
        } else {
            msg
        }
    }
)
