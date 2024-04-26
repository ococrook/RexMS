##' Function to help find good initial parameters
##'
##' @title Initialiser
##' @param uptake_guess A matrix of uptake data
##' @param timepoints A vector of timepoints
##' @param numRep The number of replicates
##' @param numtimepoints The number of timepoints
##' @param R The number of residues
##' @param phi The maximum uptake possible. Default is 1.
##' @return A matrix of initial parameters
##'
##' @md
initialiser <- function(uptake_guess,
                        timepoints,
                        numRep,
                        numtimepoints,
                        R,
                        phi) {
    params <- matrix(NA, ncol = 4, nrow = R)
    Err <- vector(mode = "numeric", length = R)
    for (j in seq_len(R)) {
        data <- cbind(c(sapply(uptake_guess, function(X) X[j, ])), rep(unique(timepoints), times = numRep))

        if (sum(is.na(data)) == 0) {
            colnames(data) <- c("y", "t")

            if (numtimepoints > 4) {
                out <- try(
                    minpack.lm::nlsLM(
                        data = as.data.frame(data),
                        formula = y ~ phi * ((1 - pi) * (1 - exp(-b * t^p)) + pi * (1 - exp(-d * t))),
                        start = list(b = 0.1, p = 0.1, d = 0.001, pi = 0.5),
                        lower = c(0, 0, 0, 0),
                        upper = c(2, 0.99, 1, 1),
                        algorithm = "LM",
                        control = nls.control(maxiter = 1000)
                    ),
                    silent = TRUE
                )
                
                if (inherits(out, "try-error")) {
                    out <- try(
                        minpack.lm::nlsLM(
                            data = as.data.frame(data),
                            formula = y ~ phi * ((1 - exp(-b * t^p))),
                            start = list(b = 0.1, p = 0.1),
                            lower = c(0, 0),
                            upper = c(2, 0.99),
                            algorithm = "LM",
                            control = nls.control(maxiter = 1000)
                        ),
                        silent = TRUE
                    )
                }
                  
                
            } else {
                # reduce initialising model if not many data points
                out <- try(
                    minpack.lm::nlsLM(
                        data = as.data.frame(data),
                        formula = y ~ phi * ((1 - exp(-b * t^p))),
                        start = list(b = 0.1, p = 0.1),
                        lower = c(0, 0),
                        upper = c(2, 0.99),
                        algorithm = "LM",
                        control = nls.control(maxiter = 1000)
                    ),
                    silent = TRUE
                )
                # force init param to be b
                init_param <- "b"
            }

            if (inherits(out, "try-error")) {
                if (j > 1) {
                    params[j, ] <- params[j - 1, ]
                    Err[j] <- NA
                } else {
                    params[j, ] <- rep(NA, 4)
                }
            } else {
                Err[j] <- sum((predict(out) - data[, 1])^2)

                if ((numtimepoints > 4)) {
                  if (length(coef(out)) == 4) {
                    params[j, ] <- coef(out)
                  } else {
                    params[j, 1:2] <- coef(out)
                    params[j, 3] <- 1
                    params[j, 4] <- 0.001
                  }
                  
                } else{
                    params[j, 1:2] <- coef(out)
                    params[j, 3] <- 1
                    params[j, 4] <- 0.001
                }
            }
        } else {
            params[j, ] <- rep(NA, 4)
        }
    }

    return(params)
}

##' Function to help find good initial parameters by guessing the initial uptake
##' using a generalised inverse appraoach
##' @title uptakeGuess
##' @param res An object of class DataFrame containing the results of the
##'   hdx-ms  experiment
##' @param numRep The number of replicates
##' @param numPeptides The number of peptides
##' @param numtimepoints The number of timepoints
##' @param R The number of residues
##' @param C A matrix indicating the position of the residues in the peptides.
##' @param phi The maximum uptake possible. Default is 1.
##' @return A matrix of intial uptake guesses
##'
##' @examples
##' require(RexMS)
##' require(dplyr)
##' data(BRD4_apo)
##'
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##' BRD4_apo <- DataFrame(BRD4_apo)
##' C <- coverageHeatmap(res = BRD4_apo, plot = FALSE)
##'
##' uptakeGuess(BRD4_apo,
##'             numRep = 3,
##'             numPeptides =  numPeptides,
##'             numtimepoints =  numTimepoints,
##'             R =  99,
##'             C = C,
##'             phi = 0.92)
##'
##' @export
uptakeGuess <- function(res,
    numRep,
    numPeptides,
    numtimepoints,
    R,
    C,
    phi = 1) {
    uptake <- lapply(
        seq.int(numRep),
        function(j) {
            matrix(res$Uptake[res$replicate == j],
                nrow = numPeptides,
                ncol = numtimepoints, byrow = TRUE
            )
        }
    )

    uptake_guess <- lapply(uptake, function(X) MASS::ginv(t(C)) %*% X)
    uptake_guess <- lapply(uptake_guess, function(X) {
        X[X > phi] <- phi
        return(X)
    })
    uptake_guess <- lapply(uptake_guess, function(X) {
        X[X < 0] <- 0
        return(X)
    })

    return(uptake_guess)
}
