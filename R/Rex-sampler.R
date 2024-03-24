##' These functions are used to implement the ReX model from Crook et al 2024.
##' using HDX-MS data using RJ-MCMC for inference.
##'
##' The 'ReX' function generates the sample from the joint posterior distribution
##' (object or class 'Rexparams') based on annotated HDX-MS data (object of class
##' [DataFrame]). Both are then passed to the `RexProcess` function to generate
##' summary data for further analysis and visualisation
##'
##' @title residue-resolved HDX-MS data analysis using the Rex yes
##' @param res An object of class DataFrame containing the HDX-MS data
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
##' @param seed An integer value indicating the random number generation seed.
##'  Default is NULL. This is used to set the seed for reproducibility.
##' @return `Rex` returns an object of class `RexParams` containing the results
##' @md
##' @rdname rex
resolver <- function(res,
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
    seed = NULL) {
    # Setting seed manually
    if (is.null(seed)) {
        message("You haven't provided a seed, you may wish to provide a seed")
    }

    # Storage
    Sigma <- vector("numeric", length = numIter)
    blongstore <- matrix(NA, nrow = R, ncol = numIter)
    pilongstore <- matrix(NA, nrow = R, ncol = numIter)
    qlongstore <- matrix(NA, nrow = R, ncol = numIter)
    dlongstore <- matrix(NA, nrow = R, ncol = numIter)
    numBreak <- vector("numeric", length = numIter)
    loglikelihoodstore <- vector("numeric", length = numIter)
    numPeptides <- length(unique(res$Sequence))
    numRep <- table(res$Sequence)[1] / numtimepoints

    # build priors from generalised inverse:
    numExch <- maxUptakes(res = res)
    index <- prepareIndexes(res = res)
    C <- coverageHeatmap(res = res, plot = FALSE)

    # make an initial guess of uptake at each residue
    uptake_guess <- uptakeGuess(
        res = res,
        numRep = numRep,
        numPeptides = numPeptides,
        numtimepoints = numtimepoints,
        R = R,
        C = C,
        phi = phi
    )


    # attempt to initalise functional paramters
    params <- initialiser(
        uptake_guess = uptake_guess,
        timepoints = timepoints,
        numRep = numRep,
        numtimepoints = numtimepoints,
        R = R,
        phi = phi
    )

    currentblong <- params[, 1]
    currentqlong <- params[, 2]
    currentdlong <- params[, 3]
    currentpilong <- params[, 4]


    # fit optimal fused lasso
    if (init_param == "d") {
        longstart <- log(currentdlong)
    } else if (init_param == "b") {
        longstart <- log(currentblong)
    }

    # clever initialisation using fused lasso based on initalisation using intial
    # parameter guesses
    naloc <- is.na(longstart) | is.infinite(longstart)
    fusedres <- genlasso::trendfilter(y = longstart[!naloc], pos = seq.int(R)[!naloc], ord = 0)
    cv <- genlasso::cv.trendfilter(fusedres)
    fusedpred <- predict(object = fusedres, pos = seq.int(R), lambda = cv$lambda.min)
    guess <- cbind(y = fusedpred$fit[, 1], seq.int(R)[!naloc])

    # current dlong
    stepfun <- approxfun(x = guess[, 2], y = guess[, 1], method = "constant")
    longstart <- stepfun(seq.int(R))
    longstart[is.na(longstart)] <- mean(longstart, na.rm = TRUE)
    longstart <- round(longstart, 6)

    currentbp <- c(0, which(diff(longstart) != 0), R)
    sub <- currentbp[-1]

    # force parameters to be consistent
    currentblong <- cleanLong(currentlong = currentblong, sub = sub, currentbp = currentbp)
    currentqlong <- cleanLong(currentlong = currentqlong, sub = sub, currentbp = currentbp)
    currentdlong <- cleanLong(currentlong = currentdlong, sub = sub, currentbp = currentbp)
    currentpilong <- cleanLong(currentlong = currentpilong, sub = sub, currentbp = currentbp)

    # short params, i.e dont include the break points
    currentb <- currentblong[sub]
    currentq <- currentqlong[sub]
    currentd <- currentdlong[sub]
    currentpi <- currentpilong[sub]

    # initilise, global uncertainty, b rate and d rate
    currentSigma <- 0.01
    rate <- 1
    drate <- 100

    # set priors
    lambda <- priors$lambda
    meanlog <- priors$meanlog
    logsd <- priors$logsd
    rho <- priors$rho
    shape1 <- priors$shape1
    shape2 <- priors$shape2
    shape <- priors$shape
    b_alpha <- priors$b_alpha
    b_beta <- priors$b_beta
    dshape <- priors$dshape
    d_alpha <- priors$d_alpha
    d_beta <- priors$d_beta
    sigma_sd <- priors$sigma_sd
    pishape1 <- priors$pishape1
    pishape2 <- priors$pishape2


    K <- length(currentbp)


    # sampler here

    for (j in seq.int(numIter)) {
        currentSigma <- metropolisSigma(
            res = res,
            currentSigma = currentSigma,
            pilong = currentpilong,
            blong = currentblong,
            qlong = currentqlong,
            dlong = currentdlong,
            mean = meanlog,
            sd = sigma_sd,
            meanprop = 0,
            sdprop = 1,
            timepoints = timepoints,
            tCoef = tCoef,
            index = index,
            numExch = numExch,
            density = density,
            phi = phi
        )
        currentb <- metropolisRate(
            res = res,
            pilong = currentpilong,
            currentb = currentb,
            qlong = currentqlong,
            dlong = currentdlong,
            currentSigma = currentSigma,
            bp = currentbp,
            mean = 0,
            sd = 0.03,
            shape = shape,
            rate = rate,
            timepoints = timepoints,
            tCoef = tCoef,
            index = index,
            numExch = numExch,
            density = density,
            phi = phi
        )
        currentblong <- rep(currentb, diff(floor(currentbp)))

        currentpi <- metropolisPi(
            res = res,
            currentpi = currentpi,
            blong = currentblong,
            qlong = currentqlong,
            dlong = currentdlong,
            currentSigma = currentSigma,
            bp = currentbp,
            mean = 0,
            sd = 0.03,
            pishape1 = pishape1,
            pishape2 = pishape2,
            timepoints = timepoints,
            tCoef = tCoef,
            index = index,
            numExch = numExch,
            density = density,
            phi = phi
        )
        currentpilong <- rep(currentpi, diff(floor(currentbp)))

        currentq <- metropolisStretch(
            res = res,
            pilong = currentpilong,
            blong = currentblong,
            currentq = currentq,
            dlong = currentdlong,
            currentSigma = currentSigma,
            bp = currentbp,
            mean = 0,
            sd = 0.03,
            shape1 = shape1,
            shape2 = shape2,
            timepoints = timepoints,
            tCoef = tCoef,
            index = index,
            numExch = numExch,
            density = density,
            phi = phi
        )

        currentqlong <- rep(currentq, diff(floor(currentbp)))

        # metropolis sampler for d.
        currentd <- metropolisdRate(
            res = res,
            pilong = currentpilong,
            blong = currentblong,
            qlong = currentqlong,
            currentd = currentd,
            currentSigma = currentSigma,
            bp = currentbp,
            mean = 0,
            sd = 0.001,
            dshape = dshape,
            drate = drate,
            timepoints = timepoints,
            tCoef = tCoef,
            index = index,
            numExch = numExch,
            density = density,
            phi = phi
        )
        currentdlong <- rep(currentd, diff(floor(currentbp)))



        for (k in seq.int(5)) {
            u <- runif(1)
            if (rho > u) {
                out <- reversibleJumpBreakPointsAdd(
                    res = res,
                    currentSigma = currentSigma,
                    currentpi = currentpi,
                    currentb = currentb,
                    currentq = currentq,
                    currentd = currentd,
                    currentbp = currentbp,
                    shape = shape,
                    rate = rate,
                    shape1 = shape1,
                    shape2, shape2,
                    K = K,
                    R = R,
                    rho = rho,
                    lambda = lambda,
                    logmean = logmean,
                    pishape1 = pishape1,
                    pishape2 = pishape2,
                    dshape = dshape,
                    drate = drate,
                    logsd = logsd,
                    multivariate = multivariate,
                    timepoints = timepoints,
                    tCoef = tCoef,
                    index = index,
                    numExch = numExch,
                    density = density,
                    R_lower = R_lower,
                    R_upper = R_upper,
                    phi = phi
                )
                currentpi <- out$currentpi
                currentb <- out$currentb
                currentq <- out$currentq
                currentd <- out$currentd
                currentbp <- out$currentbp
                K <- out$K
            } else if ((rho < u) & (K != 1)) {
                # upto here with edits
                out <- reversibleJumpBreakPointsRemove(
                    res = res,
                    currentSigma = currentSigma,
                    currentpi = currentpi,
                    currentb = currentb,
                    currentq = currentq,
                    currentd = currentd,
                    currentbp = currentbp,
                    shape = shape,
                    rate = rate,
                    shape1 = shape1,
                    shape2 = shape2,
                    K = K,
                    R = R,
                    rho = rho,
                    lambda = lambda,
                    logmean = logmean,
                    logsd = logsd,
                    pishape1 = pishape1,
                    pishape2 = pishape2,
                    dshape = dshape,
                    drate = drate,
                    multivariate = multivariate,
                    timepoints = timepoints,
                    tCoef = tCoef,
                    index = index,
                    numExch = numExch,
                    density = density,
                    phi = phi
                )
                currentb <- out$currentb
                currentpi <- out$currentpi
                currentq <- out$currentq
                currentd <- out$currentd
                currentbp <- out$currentbp
                currentd <- out$currentd
                K <- out$K
            }
        }
        # clean up break points
        sub <- which(diff(floor(currentbp)) != 0)
        currentb <- currentb[sub]
        currentpi <- currentpi[sub]
        currentq <- currentq[sub]
        currentd <- currentd[sub]
        currentbp <- unique(floor(currentbp))

        K <- length(currentbp)

        currentblong <- rep(currentb, diff(floor(currentbp)))
        currentpilong <- rep(currentpi, diff(floor(currentbp)))
        currentqlong <- rep(currentq, diff(floor(currentbp)))
        currentdlong <- rep(currentd, diff(floor(currentbp)))

        # update rate
        rate <- rgamma(
            n = 1, shape = b_alpha + length(currentb) * shape,
            rate = b_beta + sum(currentb)
        )
        drate <- rgamma(
            n = 1, shape = d_alpha + length(currentd) * dshape,
            rate = d_beta + sum(currentd)
        )

        Sigma[j] <- currentSigma
        blongstore[, j] <- currentblong
        pilongstore[, j] <- currentpilong
        qlongstore[, j] <- currentqlong
        dlongstore[, j] <- currentdlong
        numBreak[j] <- K
        loglikelihoodstore[j] <- loglik(
            res = res,
            pilong = currentpilong,
            blong = currentblong,
            qlong = currentqlong,
            dlong = currentdlong,
            sigmasq = currentSigma,
            timepoints = timepoints,
            tCoef = tCoef,
            index = index,
            numExch = numExch,
            density = density,
            phi = phi
        )
    }

    .out <- .RexChain(
        dataset = "RexExperiment",
        blong = blongstore,
        pilong = pilongstore,
        qlong = qlongstore,
        dlong = dlongstore,
        loglikelihood = loglikelihoodstore,
        phi = phi,
        UptakeGuess = uptake_guess,
        numBreak = numBreak,
        Sigma = Sigma,
        R = as.integer(R),
        timepoints = timepoints,
        numIter = as.integer(numIter),
        numPeptides = as.integer(numPeptides),
        numTimepoints = as.integer(numtimepoints)
    )


    return(.out)
}
