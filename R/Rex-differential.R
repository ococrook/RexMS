##' Functions for differential analysis of HDX-MS data
##'
##' @param HdxData A data frame containing the HDX-MS data for which you want
##'  to perform the differential analysis. Typically this dataframe contains the
##'  apo data.
##' @param params An object of class `RexParams` containing a fitted ReX model
##'  typically to the differential of interest (e.g. ligand binding)
##' @param whichChain A numeric value indicating the chain to use. Default is 1. 
##' @param num_montecarlo A numeric value indicating the number of montecarlo
##'  samples to use for the error analysis. Default is 5000.
##'
##' @return An object of class `RexDifferential` containing the results of the
##' differential analysis
##'
##' @examples
##'
##' require(RexMS)
##' require(dplyr)
##' data(BRD4_apo)
##' data(BRD4_ibet)
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' BRD4_ibet <- BRD4_ibet %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##' 
##' 
##'
##' rex_example <- rex(HdxData = DataFrame(BRD4_ibet),
##'                  numIter = 100,
##'                  R = max(BRD4_apo$End),
##'                  numtimepoints = numTimepoints,
##'                  timepoints = Timepoints,
##'                  seed = 1L,
##'                  tCoef = c(0, rep(1, 5)),
##'                  BPPARAM = SerialParam())
##'
##' rex_example <- RexProcess(HdxData = DataFrame(BRD4_ibet),
##'                          params = rex_example,
##'                          range = 50:100,
##'                          thin = 1,
##'                          whichChains = c(1,2))
##'
##' # Note change of dataset to compare with apo
##' rex_diff <- processDifferential(HdxData = DataFrame(BRD4_apo),
##'                                params = rex_example,
##'                                whichChain = 1,
##'                                num_montecarlo = 5000)
##'
##' @md
##'
##' @export
processDifferential <- function(HdxData,
                                params,
                                whichChain = 1,
                                num_montecarlo = 5000) {
    # check for comptability
    stopifnot(
        "Residues are incompatible between expreiments" =
            max(HdxData$End) <= params@chains[[whichChain]]@R
    )

    # global quantities of itnerest
    R <- params@chains[[whichChain]]@R
    phi <- params@chains@chains[[whichChain]]@phi
    numPeptides <- length(unique(HdxData$Sequence))
    numTimepoints <- length(unique(HdxData$Exposure))
    timepoints <- unique(HdxData$Exposure)
    index <- prepareIndexes(res = HdxData)
    interval <- params@interval
    Residues <- seq.int(interval[1], interval[2], by = 1)
    C <- matrix(NA, nrow = R, ncol = numPeptides)

    # storage
    ARE <- signedARE <- TRE <- matrix(nrow = numTimepoints - 1, ncol = R)
    probs <- matrix(nrow = numTimepoints - 1, ncol = R)
    eFDR <- matrix(nrow = numTimepoints - 1, ncol = 100)
    threshold <- rev(seq.int(1, 100) / 100)

    # These are computed from the params data
    blong <- params@summary@posteriorEstimates$blong
    pilong <- params@summary@posteriorEstimates$pilong
    qlong <- params@summary@posteriorEstimates$qlong
    dlong <- params@summary@posteriorEstimates$dlong
    sigma <- params@summary@Rex.globals$sigma[1]
    
    names(blong) <- names(pilong) <- names(qlong) <- names(dlong) <- Residues

    err <- error_prediction(
        res = HdxData,
        blong = blong,
        phi = phi,
        qlong = qlong,
        pilong = pilong,
        dlong = dlong
    )

    for (i in seq_len(numTimepoints)[-1]) {
        for (j in seq_len(ncol(C))) {
            C[index[[j]], j] <- err[j, i]
        }
        C2 <- (t((t(C) / sapply(index, length))))
        redundancy <- colSums(apply(C, 1, function(x) !is.na(x)))
        lplace <- sapply(seq_along(redundancy), function(x) {
            abs(rowSums(matrix(
                LaplacesDemon::rlaplace(
                    n = num_montecarlo * redundancy[x],
                    scale = sqrt(sigma)
                ),
                ncol = redundancy[x]
            )))
        })


        signedARE[i - 1, ] <- rowMeans(C2, na.rm = TRUE)
        ARE[i - 1, ] <- rowMeans(abs(C2), na.rm = TRUE)
        TRE[i - 1, ] <- rowSums(C2, na.rm = TRUE)
        probs[i - 1, ] <- sapply(seq_along(redundancy), function(x) {
            sum(abs(rowSums(C2, na.rm = TRUE))[x] > abs(lplace[[x]])) / num_montecarlo
        })
        eFDR[i - 1, ] <- sapply(threshold, function(x) EFDR(probs[i - 1, ], threshold = x))
    }
    totalprobs <- colProds(probs)

    Rex.predictionError <- DataFrame(unique(HdxData$Sequence), err)
    colnames(Rex.predictionError) <- c("Peptide", timepoints)

    Rex.estimates <- DataFrame(
        Residues = Residues,
        signedARE = t(signedARE)[Residues,],
        ARE = t(ARE)[Residues,],
        TRE = t(TRE)[Residues,]
    )

    colnames(Rex.estimates) <- c(
        "Resdiues",
        paste0(
            rep(c("signedARE_", "ARE_", "TRE_"),
                each = numTimepoints - 1
            ),
            timepoints[-1]
        )
    )

    Rex.probs <- DataFrame(
        Residues = Residues,
        probs = t(probs)[Residues,],
        totalprobs = totalprobs[Residues]
    )

    colnames(Rex.probs) <- c(
        "Residues",
        paste0(
            rep(c("probs_"),
                each = numTimepoints - 1
            ),
            timepoints[-1]
        ), "totalprobs"
    )

    Rex.eFDR <- DataFrame(
        Threshold = threshold,
        eFDR = t(eFDR)
    )

    colnames(Rex.eFDR) <- c(
        "Threshold",
        paste0(
            rep(c("eFDR_"),
                each = numTimepoints - 1
            ),
            timepoints[-1]
        )
    )




    .out <- .RexDifferential(
        dataset = "Differential Rex Experiment",
        Rex.predictionError = Rex.predictionError,
        Rex.estimates = Rex.estimates,
        Rex.probs = Rex.probs,
        Rex.eFDR = Rex.eFDR
    )

    return(.out)
}

##' Function to calculate the corresponding expected False discovery rate
##' (eFDR) for a given threshold
##'
##' @param prob A numeric vector of probabilities
##' @param threshold A numeric value between 0 and 1
##' @return A numeric value indicating the expected FDR
##'
##' @examples
##' EFDR(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 0.9)
##'
##' @md
##' @export
EFDR <- function(prob, threshold = 0.90) {
    stopifnot("prob must be numeric" = is(prob, "numeric"))
    stopifnot("prob must be probabilities" = max(prob) <= 1)
    stopifnot("prob must be probabilities" = min(prob) >= 0)
    stopifnot("threshold must be numeric" = is(threshold, "numeric"))
    stopifnot("threshold must be a single values" = length(threshold) == 1)

    .out <- sum((1 - prob) * I(prob >= threshold)) / sum(I(prob >= threshold))
    return(.out)
}
##' Function to obtain uncertainty estimates for the total relative error (TRE)
##'
##' @param HdxData A data frame containing the HDX-MS data for which you want
##' to perform the differential analysis. Typically this dataframe contains the
##' apo data.
##' @param params An object of class `RexParams` containing a fitted ReX model
##' typically to the differential of interest (e.g. ligand binding)
##' @param num_montecarlo A numeric value indicating the number of montecarlo
##' samples to use for the error analysis. Default is 5000.
##' @param whichChain A numeric value indicating the chain to use. Default is 1.
##' @param whichSamples A numeric vector indicating which samples to use. Default
##' is seq.int(50).
##' @return An array of the total relative error (TRE) for each residue at each
##' timepoint represnting uncertainty in the distribution
##'
##' @examples
##' require(RexMS)
##' require(dplyr)
##' data(BRD4_apo)
##' data(BRD4_ibet)
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' BRD4_ibet <- BRD4_ibet %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##' 
##' rex_example <- rex(HdxData = DataFrame(BRD4_ibet),
##'                 numIter = 100,
##'                 R = max(BRD4_apo$End),
##'                 numtimepoints = numTimepoints,
##'                 timepoints = Timepoints,
##'                 seed = 1L,
##'                 tCoef = c(0, rep(1, 5)),
##'                 BPPARAM = SerialParam())
##'
##' rex_example <- RexProcess(HdxData = DataFrame(BRD4_ibet),
##'                         params = rex_example,
##'                         range = 50:100,
##'                         thin = 1,
##'                         whichChains = c(1,2))
##'
##' rex_TRE_uncertainty <- processTREuncertainty(HdxData = DataFrame(BRD4_apo),
##'                                            params = rex_example,
##'                                            whichChain = 1,
##'                                            num_montecarlo = 5000)
##'
##' @md
##'
##' @export
processTREuncertainty <- function(HdxData,
                                  params,
                                  whichChain = 1,
                                  whichSamples = seq.int(50),
                                  num_montecarlo = 5000) {
    # global quantities of itnerest
    R <- params@chains[[whichChain]]@R
    phi <- params@chains@chains[[whichChain]]@phi
    numPeptides <- length(unique(HdxData$Sequence))
    numTimepoints <- length(unique(HdxData$Exposure))
    timepoints <- unique(HdxData$Exposure)
    index <- prepareIndexes(res = HdxData)
    interval <- params@interval
    Residues <- seq.int(interval[1], interval[2], by = 1)
    C <- matrix(NA, nrow = R, ncol = numPeptides)
    TRE <- array(NA, c(length(whichSamples), numTimepoints - 1, R))


    for (i in whichSamples) {
        blong <- params@chains@chains[[whichChain]]@blong[, i]
        pilong <- params@chains@chains[[whichChain]]@pilong[, i]
        qlong <- params@chains@chains[[whichChain]]@qlong[, i]
        dlong <- params@chains@chains[[whichChain]]@dlong[, i]

        err <- error_prediction(
            res = HdxData,
            blong = blong,
            phi = phi,
            qlong = qlong,
            pilong = pilong,
            dlong = dlong
        )

        for (k in seq_len(numTimepoints)[-1]) {
            for (j in seq_len(ncol(C))) {
                C[index[[j]], j] <- err[j, k]
            }
            C2 <- (t((t(C) / sapply(index, length))))
            redundancy <- colSums(apply(C, 1, function(x) !is.na(x)))
            lplace <- sapply(seq_along(redundancy), function(x) {
                abs(rowSums(matrix(
                    LaplacesDemon::rlaplace(
                        n = num_montecarlo * redundancy[x],
                        scale = sqrt(params@chains@chains[[whichChain]]@Sigma[i])
                    ),
                    ncol = redundancy[x]
                )))
            })

            TRE[which(whichSamples == i), k - 1, ] <- rowSums(C2, na.rm = TRUE)
        }
    }

    .out <- TRE

    return(.out)
}
