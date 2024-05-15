##' These functions are used to implement the ReX model from Crook et al 2024.
##' using HDX-MS data using RJ-MCMC for inference.
##'
##' The 'ReX' function generates the sample from the joint posterior distribution
##' (object or class 'Rexparams') based on annotated HDX-MS data (object of class
##' [DataFrame]). Both are then passed to the `RexProcess` function to generate
##' summary data for further analysis and visualisation
##'
##' @title residue-resolved HDX-MS data analysis using Rex
##' @param HdxData An object of class DataFrame containing the HDX-MS data
##' @param numIter The number of iterations for the MCMC. Default is 1000
##'  but recommend at least 5000 in practice
##' @param R The number of residues in the protein
##' @param numtimepoints The number of timepoints in the HDX-MS experiment
##' @param multivariate A logical indicating whether to use a multivariate
##'  model. Default is FALSE.
##' @param timepoints A numeric vector of timepoints in the HDX-MS experiment.
##'  Default is c(0, 30, 300).
##' @param tCoef A numeric value for the timpeoint standard deviation
##' coefficient. Default is 1. This
##'  allows for the timepoint standard deviation to vary by a fixed factor.
##'   1 indicates the same for
##'  each timepoint.
##' @param density A character string indicating the density function to use.
##'  Default is "Gaussian" but
##'  typically "laplace" is also used.
##' @param R_lower A numeric value indicating the lower bound for the residue
##'  number. Default is 1. This
##'  may change if you want to subset the protien or avoid a his-tag.
##' @param R_upper A numeric value indicating the upper bound for the residue
##'  number. Default is R. This
##' may change if you want to subset the protien.
##' @param priors A list of prior parameters for the model
##' @param phi A numeric value indicating the maximum uptake possible.
##'  Default is 0.92. Typically the
##' deuterium buffer concentration of the HDX-MS experiment.
##' @param init_param A character string indicating the initial parameter to
##'  use. Default is "d" but "b" is also possible. This is used to help find
##'  good initial parameters by guessing the initial uptake and is used
##'   in the `uptakeGuess` function. We do not recommend changing
##'    this unless you are an expert.
##' @param numChains A numeric value indicating the number of chains to run. Default
##'  is 2. This is used to run the chains in parallel.
##' @param seed An integer value indicating the random number generation seed.
##'  Default is NULL. This is used to set the seed for reproducibility.
##' @param BPPARAM A BiocParallelParam object indicating the parallel backend to use.
##'  Default is BiocParallel::bpparam(). This is used to run the chains in parallel.
##' @return `Rex` returns an object of class `RexParams` containing the results
##' @md
##' @examples
##' require(RexMS)
##' require(dplyr)
##' data("BRD4_apo")
##' BRD4_apo <- BRD4_apo %>% filter(End < 40)
##'
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
## numPeptides <- length(unique(BRD4_apo$Sequence))
##'
##' rex_example <- rex(HdxData = DataFrame(BRD4_apo),
##'                 numIter = 4, # typically much larger
##'                 R = max(BRD4_apo$End),
##'                 numtimepoints = numTimepoints,
##'                 timepoints = Timepoints,
##'                 seed = 1L,
##'                 numChains = 1L,
##'                 tCoef = c(0, rep(1, 5)),
##'                 BPPARAM = SerialParam())
##'
##'
##' @rdname rex
##' @export
rex <- function(HdxData,
                numIter = 1000,
                R = 379,
                numtimepoints = 3,
                multivariate = FALSE,
                timepoints = c(0, 30, 300),
                tCoef = 1,
                density = "Gaussian",
                R_lower = 1,
                R_upper = R,
                priors = list(
                    lambda = 100 / (R - 1),
                    meanlog = -3,
                    sdlog = 1,
                    rho = 0.5,
                    shape1 = 1,
                    shape2 = 5,
                    shape = 1,
                    b_alpha = 1,
                    b_beta = 200,
                    dshape = 1,
                    d_alpha = 1,
                    d_beta = 1,
                    sigma_sd = 0.5,
                    pishape1 = 1,
                    pishape2 = 10
                ),
                phi = 0.92,
                init_param = "d",
                numChains = 2L,
                seed = NULL,
                BPPARAM = BiocParallel::bpparam()) {
    # checks
    stopifnot(exprs = {
        "HdxData must be a DataFrame" <- is(HdxData, "DFrame")
        "numIter must be a numeric" <- is.numeric(numIter)
        "R must be a numeric" <- is.numeric(R)
        "numtimepoints must be a numeric" <- is.numeric(numtimepoints)
        "multivariate must be a logical" <- is.logical(multivariate)
        "timepoints must be a numeric" <- is.numeric(timepoints)
        "tCoef must be a numeric" <- is.numeric(tCoef)
        "density must be a character" <- is.character(density)
        "R_lower must be a numeric" <- is.numeric(R_lower)
        "R_upper must be a numeric" <- is.numeric(R_upper)
        "priors must be a list" <- is.list(priors)
        "phi must be a numeric" <- is.numeric(phi)
        "init_param must be a character" <- is.character(init_param)
    })

    # check the priors
    stopifnot(exprs = {
        "lambda must be a numeric" <- is.numeric(priors$lambda)
        "meanlog must be a numeric" <- is.numeric(priors$meanlog)
        "sdlog must be a numeric" <- is.numeric(priors$sdlog)
        "rho must be a numeric" <- is.numeric(priors$rho)
        "shape1 must be a numeric" <- is.numeric(priors$shape1)
        "shape2 must be a numeric" <- is.numeric(priors$shape2)
        "shape must be a numeric" <- is.numeric(priors$shape)
        "b_alpha must be a numeric" <- is.numeric(priors$b_alpha)
        "b_beta must be a numeric" <- is.numeric(priors$b_beta)
        "dshape must be a numeric" <- is.numeric(priors$dshape)
        "d_alpha must be a numeric" <- is.numeric(priors$d_alpha)
        "d_beta must be a numeric" <- is.numeric(priors$d_beta)
        "sigma_sd must be a numeric" <- is.numeric(priors$sigma_sd)
        "pishape1 must be a numeric" <- is.numeric(priors$pishape1)
        "pishape2 must be a numeric" <- is.numeric(priors$pishape2)
    })

    stopifnot({
        "multivariate not currently implemented" <- isFALSE(multivariate)
    })

    # check the data
    HdxData <- cleanHDX(res = HdxData, clean = FALSE)

    # check missing values
    if (sum(is.na(HdxData$Uptake)) > 0) {
        stop("The dataset contains missing values. Please examine the dataset first")
    }

    ## chains run in parallel, repeating number of iterations
    .res <- BiocParallel::bplapply(rep(numIter, numChains),
        FUN = resolver,
        res = HdxData,
        # numIter = numIter, #parrallelised
        R = R,
        numtimepoints = numtimepoints,
        multivariate = multivariate,
        timepoints = timepoints,
        tCoef = tCoef,
        density = density,
        R_lower = R_lower,
        R_upper = R_upper,
        priors = priors,
        phi = phi,
        init_param = init_param,
        seed = seed
    )


    # construct the RexChains
    .chains <- .RexChains(chains = .res)

    ## construct the RexParams
    .out <- .RexParams(
        method = "ReX",
        priors = priors,
        seed = as.integer(seed),
        summary = .RexSummary(),
        chains = .chains,
        interval = c(R_lower, R_upper)
    )



    return(.out)
}

##' @title proccess Rex results
##' @param HdxData An object of class `DataFrame` containing the HDX-MS data
##' @param params An object of class `RexParams`
##' @param range The iterations to keep. Default is seq.int(4000, 5000, by = thin)
##' @param thin The thinning interval. Default is 1
##' @param whichChains The chains to keep. Default is c(1,2)
##' @return An object of class `RexParams` with its summary slot populated
##'
##' @examples
##' require(RexMS)
##' require(dplyr)
##' data("BRD4_apo")
##' BRD4_apo <- BRD4_apo %>% filter(End < 40)
##'
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
## numPeptides <- length(unique(BRD4_apo$Sequence))
##'
##' rex_example <- rex(HdxData = DataFrame(BRD4_apo),
##'                 numIter = 4, # typically much larger
##'                 R = max(BRD4_apo$End),
##'                 numtimepoints = numTimepoints,
##'                 timepoints = Timepoints,
##'                 seed = 1L,
##'                 numChains = 1L,
##'                 tCoef = c(0, rep(1, 5)),
##'                 BPPARAM = SerialParam())
##' rex_example <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                           params = rex_example,
##'                           thin = 1,
##'                           range = 1:4,
##'                           whichChains = c(1))
##' @md
##' @rdname rex-process
##' @export
RexProcess <- function(HdxData,
                       params,
                       thin = 1,
                       range = seq.int(4000, 5000, by = thin),
                       whichChains = c(1, 2)) {
    ## get required slots
    numChains <- length(params@chains)
    R <- params@chains[[1]]@R
    numIter <- params@chains[[1]]@numIter
    phi <- params@chains[[1]]@phi
    numPeptides <- params@chains[[1]]@numPeptides
    numTimepoints <- params@chains[[1]]@numTimepoints
    interval <- params@interval
    Residues <- seq.int(interval[1], interval[2], by = 1)
    timepoints <- params@chains[[1]]@timepoints


    stopifnot("range must be compartible with numIter" = max(range) <= numIter)
    stopifnot("whichChains must be a numeric" = is.numeric(whichChains))
    stopifnot("whichChains must be compartible with numChains" = max(whichChains) <= numChains)

    # checks
    stopifnot(exprs = {
        "params must be a RexParams" <- is(params, "RexParams")
        "thin must be a numeric" <- is.numeric(thin)
    })

    # check the chains
    stopifnot(exprs = {
        "chains must be a RexChains" <- is(params@chains, "RexChains")
        "summary must be a RexSummary" <- is(params@summary, "RexSummary")
        "thin must be a numeric" <- is.numeric(thin)
    })

    # check the range
    stopifnot(exprs = {
        "range must be a numeric" <- is.numeric(range)
    })

    ## storage
    mean_blong <- matrix(nrow = length(whichChains), ncol = R)
    mean_pilong <- matrix(nrow = length(whichChains), ncol = R)
    mean_qlong <- matrix(nrow = length(whichChains), ncol = R)
    mean_dlong <- matrix(nrow = length(whichChains), ncol = R)
    C <- matrix(NA, nrow = R, ncol = numPeptides)


    ## sigma
    sigma <- vector(mode = "numeric", length = length(whichChains))
    likelihood <- vector(mode = "numeric", length = length(whichChains))

    ## average over mcmc interations
    for (i in seq_along(whichChains)) {
        j <- whichChains[i]
        mean_blong[i, ] <- rowMeans(params@chains[[j]]@blong[, range])
        mean_pilong[i, ] <- rowMeans(params@chains[[j]]@pilong[, range])
        mean_qlong[i, ] <- rowMeans(params@chains[[j]]@qlong[, range])
        mean_dlong[i, ] <- rowMeans(params@chains[[j]]@dlong[, range])
        sigma[i] <- mean(params@chains[[j]]@Sigma[range])
        likelihood[i] <- mean(params@chains[[j]]@loglikelihood[range])
    }

    # store quantiles
    for (i in seq_along(whichChains)) {
        k <- whichChains[i]
        .blong_quantiles <- vapply(seq_len(R),
            function(j) {
                quantile(
                    params@chains[[k]]@blong[j, range],
                    c(0.025, 0.5, 0.975)
                )
            },
            FUN.VALUE = numeric(3)
        )
        .pilong_quantiles <- vapply(seq_len(R),
            function(j) {
                quantile(
                    params@chains[[k]]@pilong[j, range],
                    c(0.025, 0.5, 0.975)
                )
            },
            FUN.VALUE = numeric(3)
        )
        .qlong_quantiles <- vapply(seq_len(R),
            function(j) {
                quantile(
                    params@chains[[k]]@qlong[j, range],
                    c(0.025, 0.5, 0.975)
                )
            },
            FUN.VALUE = numeric(3)
        )
        .dlong_quantiles <- vapply(seq_len(R),
            function(j) {
                quantile(
                    params@chains[[k]]@dlong[j, range],
                    c(0.025, 0.5, 0.975)
                )
            },
            FUN.VALUE = numeric(3)
        )
    }

    ## take means across chains
    blong <- colMeans(mean_blong[, , drop = FALSE])
    pilong <- colMeans(mean_pilong[, , drop = FALSE])
    qlong <- colMeans(mean_qlong[, , drop = FALSE])
    dlong <- colMeans(mean_dlong[, , drop = FALSE])

    .diagnoistics <- matrix(NA, 1, 1)
    .sigma_diagnostics <- vector("list", length = length(whichChains))
    for (i in seq_along(whichChains)) {
        j <- whichChains[i]
        .sigma_diagnostics[[i]] <- coda::as.mcmc(params@chains[[j]]@Sigma[range])
    }

    if (length(whichChains) > 1) {
        sigma_mcmc <- coda::as.mcmc.list(.sigma_diagnostics)
        gd <- coda::gelman.diag(sigma_mcmc)
        .diagnoistics <- gd$psrf
        rownames(.diagnoistics) <- c("sigma")
    }

    names(blong) <- names(pilong) <- names(qlong) <- names(dlong) <- seq.int(R)

    # make error predictions
    err <- error_prediction(
        res = HdxData,
        blong = blong,
        pilong = pilong,
        qlong = qlong,
        dlong = dlong,
        phi = phi
    )

    ## compute resolution metrics
    index <- prepareIndexes(res = HdxData)
    ARE <- signedARE <- TRE <- matrix(nrow = numTimepoints - 1, ncol = R)
    for (i in seq_len(numTimepoints)[-1]) {
        for (j in seq_len(ncol(C))) {
            C[index[[j]], j] <- err[j, i]
        }
        C2 <- (t((t(C) / sapply(index, length))))
        redundancy <- colSums(apply(C, 1, function(x) !is.na(x)))
        signedARE[i - 1, ] <- rowMeans(C2, na.rm = TRUE)
        ARE[i - 1, ] <- rowMeans(abs(C2), na.rm = TRUE)
        TRE[i - 1, ] <- rowSums(C2, na.rm = TRUE)
    }

    Rex.resolution <- DataFrame(
        Residues = Residues,
        signedARE = t(signedARE)[Residues,],
        ARE = t(ARE)[Residues,],
        TRE = t(TRE)[Residues,]
    )

    colnames(Rex.resolution) <- c(
        "Resdiues",
        paste0(rep(c("signedARE_", "ARE_", "TRE_"), each = numTimepoints - 1), timepoints[-1])
    )

    Rex.peptideError <- DataFrame(unique(HdxData$Sequence), err)

    colnames(Rex.peptideError) <- c("Peptide", timepoints)

    ## store the posterior estimates
    posteriorEstimates <- DataFrame(
        Residues = Residues,
        blong = blong[Residues],
        pilong = pilong[Residues],
        qlong = qlong[Residues],
        dlong = dlong[Residues]
    )

    colnames(posteriorEstimates) <- c(
        "Residues",
        "blong",
        "pilong",
        "qlong",
        "dlong"
    )


    Rex.quantiles <- DataFrame(
        Residues = Residues,
        blong_quantiles = t(.blong_quantiles)[Residues,],
        pilong_quantiles = t(.pilong_quantiles)[Residues,],
        qlong_quantiles = t(.qlong_quantiles)[Residues,],
        dlong_quantiles = t(.dlong_quantiles)[Residues,]
    )

    colnames(Rex.quantiles) <- c(
        "Residues",
        paste0(rep(c(
            "blong_quantiles",
            "pilong_quantiles",
            "qlong_quantiles",
            "dlong_quantiles"), each = 3
        ), rep(c("_2.5%", "_50%", "_97.5%"), times = 4))
    )


    Rex.globals <- DataFrame(
        sigma = sigma,
        likelihood = likelihood
    )
    colnames(Rex.globals) <- c("sigma", "likelihood")


    ## store the summary
    params@summary <- .RexSummary(
        posteriorEstimates = posteriorEstimates,
        diagnostics = .diagnoistics,
        Rex.quantiles = Rex.quantiles,
        Rex.globals = Rex.globals,
        Rex.resolution = Rex.resolution,
        Rex.peptideError = Rex.peptideError
    )

    return(params)
}
