##' @param object An instance of appropriate class.
##' @rdname RexParams
chains <- function(object) {
    stopifnot(inherits(object, "RexParams"))
    object@chains
}

##' @md
##' @rdname RexParams
setMethod(
    "show", "RexParams",
    function(object) {
        cat("Object of class \"", class(object), "\"\n", sep = "")
        cat("Method:", object@method, "\n")
        cat("Number of chains:", length(object@chains), "\n")
        invisible(NULL)
    }
)

##' @rdname RexParams
setMethod(
    "show", "RexChain",
    function(object) {
        cat(" Object of class \"", class(object), "\"\n", sep = "")
        cat(" Number of residues:", object@R, "\n")
        cat(" Number of timepoints:", object@numTimepoints, "\n")
        cat(" Number of peptides:", object@numPeptides, "\n")
        cat(" Number of iterations:", object@numIter, "\n")
        invisible(NULL)
    }
)

##' @rdname RexParams
setMethod(
    "show", "RexSummary",
    function(object) {
        cat(" Object of class \"", class(object), "\"\n", sep = "")
        cat(" Number of Residues:", nrow(object@posteriorEstimates), "\n")
        cat(" Number of Timepoints:", ncol(object@Rex.peptideError) - 1, "\n")
        cat(" Number of Peptides:", nrow(object@Rex.peptideError), "\n")
        cat(" Number of chains:", nrow(object@Rex.globals), "\n")
        cat(" Global resolutions:", object@Rex.globals$sigma, "\n")
        invisible(NULL)
    }
)

##' @rdname RexParams
setMethod(
    "length", "RexChains",
    function(x) length(x@chains)
)

##' @rdname RexParams
setMethod(
    "length", "RexParams",
    function(x) length(chains(x))
)


##' @param object An instance of appropriate class.
##' @rdname RexParams
posteriorEstimates <- function(object) {
    stopifnot(inherits(object, "RexSummary"))
    object@posteriorEstimates
}

##' @rdname RexParams
setMethod(
    "posteriorEstimates", "RexSummary",
    function(object) object@posteriorEstimates
)

##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.quantiles <- function(object) {
    stopifnot(inherits(object, "RexSummary"))
    object@Rex.quantiles
}

##' @rdname RexParams
setMethod(
    "Rex.quantiles", "RexSummary",
    function(object) object@Rex.quantiles
)

##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.resolution <- function(object) {
    stopifnot(inherits(object, "RexSummary"))
    object@Rex.resolution
}

##' @rdname RexParams
setMethod(
    "Rex.resolution", "RexSummary",
    function(object) object@Rex.resolution
)

##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.peptideError <- function(object) {
    stopifnot(inherits(object, "RexSummary"))
    object@Rex.peptideError
}

##' @rdname RexParams
setMethod(
    "Rex.peptideError", "RexSummary",
    function(object) object@Rex.peptideError
)

##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.globals <- function(object) {
    stopifnot(inherits(object, "RexSummary"))
    object@Rex.globals
}

##' @rdname RexParams
setMethod(
    "Rex.globals", "RexSummary",
    function(object) object@Rex.globals
)


##' @rdname RexParams
setMethod(
    "posteriorEstimates", "RexParams",
    function(object) posteriorEstimates(object@summary)
)


##' @rdname RexParams
setMethod(
    "Rex.quantiles", "RexParams",
    function(object) Rex.quantiles(object@summary)
)


##' @rdname RexParams
setMethod(
    "Rex.resolution", "RexParams",
    function(object) Rex.resolution(object@summary)
)

##' @rdname RexParams
setMethod(
    "Rex.peptideError", "RexParams",
    function(object) Rex.peptideError(object@summary)
)


##' @rdname RexParams
setMethod(
    "Rex.globals", "RexParams",
    function(object) Rex.globals(object@summary)
)

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##'
##' @md
##' @rdname RexParams
setMethod(
    "[[", "RexChains",
    function(x, i, j = "missing", drop = "missing") x@chains[[i]]
)

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##'
##' @md
##' @rdname RexParams
setMethod(
    "[", "RexChains",
    function(x, i, j = "missing", drop = "missing") {
        if (any(i > length(x))) {
            stop("Index out of bounds. Only ", length(x), " chain(s) available.")
        }
        x@chains <- x@chains[i]
        x
    }
)
##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##'
##' @md
##' @rdname RexParams
setMethod(
    "[", "RexParams",
    function(x, i, j = "missing", drop = "missing") {
        if (any(i > length(x))) {
            stop("Index out of bounds. Only ", length(x), " chain(s) available.")
        }
        x@chains <- chains(x)[i]
        x
    }
)
##' @param object of RexChains
##' @rdname RexParams
setMethod(
    "show", "RexChains",
    function(object) {
        cat(" Object of class \"", class(object), "\"\n", sep = "")
        cat(" Number of chains:", length(object), "\n")
        invisible(NULL)
    }
)

##' @param object of class RexDifferential
##' @rdname RexParams
setMethod(
    "show", "RexDifferential",
    function(object) {
        cat(" Object of class \"", class(object), "\"\n", sep = "")
        invisible(NULL)
    }
)
##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.predictionError <- function(object) {
    stopifnot(inherits(object, "RexDifferential"))
    object@Rex.predictionError
}

##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.estimates <- function(object) {
    stopifnot(inherits(object, "RexDifferential"))
    object@Rex.estimates
}

##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.probs <- function(object) {
    stopifnot(inherits(object, "RexDifferential"))
    object@Rex.probs
}

##' @param object An instance of appropriate class.
##' @rdname RexParams
Rex.eFDR <- function(object) {
    stopifnot(inherits(object, "RexDifferential"))
    object@Rex.eFDR
}

##' @rdname RexParams
setMethod(
    "Rex.predictionError", "RexDifferential",
    function(object) object@Rex.predictionError
)

##' @rdname RexParams
setMethod(
    "Rex.estimates", "RexDifferential",
    function(object) object@Rex.estimates
)

##' @rdname RexParams
setMethod(
    "Rex.probs", "RexDifferential",
    function(object) object@Rex.probs
)

##' @rdname RexParams
setMethod(
    "Rex.eFDR", "RexDifferential",
    function(object) object@Rex.eFDR
)
