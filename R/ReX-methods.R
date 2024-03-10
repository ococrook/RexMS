##'@param object An instance of appropriate class.
##'@rdname RexParams
chains <- function(object) {
  stopifnot(inherits(object, "RexParams"))
  object@chains
}

##' @md
##' @rdname RexParams
setMethod("show", "RexParams",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat("Method:", object@method, "\n")
            cat("Number of chains:", length(object@chains), "\n")
            invisible(NULL)
          })

##' @rdname RexParams
setMethod("show", "RexChain",
          function(object) {
            cat(" Object of class \"", class(object), "\"\n", sep = "")
            cat(" Number of residues:", object@R, "\n")
            cat(" Number of timepoints:", object@numTimepoints, "\n")
            cat(" Number of peptides:", object@numPeptides, "\n")
            cat(" Number of iterations:", object@numIter, "\n")
            invisible(NULL)
          })

##' @rdname RexParams
setMethod("length", "RexChains",
          function(x) length(x@chains))

##' @rdname RexParams
setMethod("length", "RexParams",
          function(x) length(chains(x)))

##' @rdname RexParams
setMethod("length", "RexSummaries",
          function(x) length(summaries(x)))

##'@param object An instance of appropriate class.
##'@rdname RexParams
posteriorEstimates <- function(object) {
  stopifnot(inherits(object, "RexSummary"))
  object@posteriorEstimates
}

##' @rdname RexParams
setMethod("posteriorEstimates", "RexSummary",
          function(object) object@posteriorEstimates)

##'@param object An instance of appropriate class.
##'@rdname RexParams
RexJoint <- function(object) {
  stopifnot(inherits(object, "RexSummary"))
  object@Rex.joint
}

##'@rdname RexParams
setMethod("RexJoint", "RexSummary",
          function(object) object@Rex.joint)


##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname RexParams
setMethod("[[", "RexChains",
          function(x, i, j = "missing", drop = "missing") x@chains[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname RexParams
setMethod("[[", "RexParams",
          function(x, i, j = "missing", drop = "missing") params(x)[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname RexParams
setMethod("[", "RexChains",
          function(x, i, j = "missing", drop = "missing") {
            if (any(i > length(x)))
              stop("Index out of bounds. Only ", length(x), " chain(s) available.")
            x@chains <- x@chains[i]
            x
          })
##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname RexParams
setMethod("[", "RexParams",
          function(x, i, j = "missing", drop = "missing") {
            if (any(i > length(x)))
              stop("Index out of bounds. Only ", length(x), " chain(s) available.")
            x@chains <- chains(x)[i]
            x
          })
##' @param object of RexChains
##' @rdname RexParams
setMethod("show", "RexChains",
          function(object) {
            cat(" Object of class \"", class(object), "\"\n", sep = "")
            cat(" Number of chains:", length(object), "\n")
            invisible(NULL)
          })


