##' Function that computes errors between the predicted peptide values and
##' the observed peptide values
##' @param res An object of class DataFrame containing the results of the
##'  hdx-ms  experiment
##' @param blong A matrix of dimensions R by 1 containing the blong values
##' @param pilong A matrix of dimensions R by 1 containing the pilong values
##' @param qlong A matrix of dimensions R by 1 containing the qlong values
##' @param dlong A matrix of dimensions R by 1 containing the dlong values
##' @param phi A numeric value containing the phi value
##' @return Returns exchange indexes as a list
##'
error_prediction <- function(res, blong, pilong, qlong, dlong, phi) {
    timepoints <- unique(res$Exposure)
    out_res <- sapply(seq.int(length(blong)), function(z) {
        doublelogisticFunction(timepoints,
            b = blong[z],
            a = phi, q = qlong[z],
            pi = pilong[z], d = dlong[z]
        )
    })
    numtimepoints <- length(unique(res$Exposure))
    numRep <- table(res$Sequence)[1] / numtimepoints


    diff_coupling <- matrix(NA, nrow = length(unique(res$Sequence)), ncol = numtimepoints)

    reduce_rep <- c(rep(FALSE, numRep - 1), TRUE)
    for (j in seq_len(length(unique(res$Sequence)))) {
        start <- res$Start[res$Sequence == unique(res$Sequence)[j]][1]
        end <- res$End[res$Sequence == unique(res$Sequence)[j]][1]
        numExch <- res$MaxUptake[res$Sequence == unique(res$Sequence)[j]][1]
        seq <- strsplit(unique(res$Sequence)[j], "")[[1]]
        seq <- seq[-c(1:2)]
        index <- which(seq != "P") + start + 1
        mu <- rowSums(out_res[, index])

        diff_coupling[j, ] <- mu - colMeans(matrix(res$Uptake[res$Sequence == unique(res$Sequence)[j]], nrow = numRep))
    }

    return(diff_coupling)
}
