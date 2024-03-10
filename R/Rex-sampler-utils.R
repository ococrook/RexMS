logisticFunction <- function(t, b, a, q = 1){
  
  out <- a * (1 - exp(- b * t^q))
  return(out)
}

doublelogisticFunction <- function(t, b, a = 1, q = 1, pi, d){
  out <- a * ( (1 - pi) * (1 - exp(- b * t^q)) + pi * (1 - exp(- d * t)))
  return(out)
}


loglik <- function(res,
                   pilong,
                   blong,
                   qlong,
                   dlong,
                   sigmasq,
                   timepoints,
                   tCoef = 1,
                   index,
                   numExch,
                   density = "Gaussian",
                   phi){
  
  num_timepoints <- length(timepoints)
  loglik <- vector(mode = "numeric", length = length(unique(res$Sequence)))
  numRep <- table(res$Sequence)[1]/num_timepoints
  d <- matrix(NA, nrow = length(unique(res$Sequence)), ncol =  num_timepoints * numRep)
  .sd <- matrix(NA, nrow = length(unique(res$Sequence)), ncol = length(tCoef) * numRep)
  
  for (j in seq_along(unique(res$Sequence))){
    
    br <- blong[index[[j]]]
    pir <- pilong[index[[j]]]
    qr <- qlong[index[[j]]]
    dr <- dlong[index[[j]]]
    mu <- rowSums(sapply(seq.int(length(br)),
                         function(z) 
                           doublelogisticFunction(t = timepoints,
                                                  b = br[z],
                                                  a = phi,
                                                  q = qr[z], 
                                                  pi = pir[z],
                                                  d = dr[z])))
    d[j,] <- res$Uptake[res$Sequence == unique(res$Sequence)[j]] - rep(mu, each = numRep)
    
    .sd[j,] <- rep(tCoef * numExch[[j]] * sqrt(sigmasq), each = numRep)
  }
  
  if (density == "Gaussian"){
    if (ncol(.sd) > 1){
      
      rmv_zero <- seq.int(numRep)
      
      loglik <- vapply(seq.int(nrow(d)),
                       function(z) sum(dnorm(d[z, -rmv_zero], sd = .sd[z, -rmv_zero], log = TRUE)),
                       FUN.VALUE = numeric(1))
    } else{
      loglik <- vapply(seq.int(nrow(d)),
                       function(z) sum(dnorm(d[z,], sd = .sd[z,], log = TRUE)),
                       FUN.VALUE = numeric(1))
    }
  } else{
    if (ncol(.sd) > 1){
      
      rmv_zero <- seq.int(numRep)
      
      loglik <- vapply(seq.int(nrow(d)),
                       function(z) 
                         sum(LaplacesDemon::dlaplace(x = d[z, -rmv_zero], scale = .sd[z, -rmv_zero], log = TRUE)),
                       FUN.VALUE = numeric(1))
    } else{
      loglik <- vapply(seq.int(nrow(d)),
                       function(z) 
                         sum(LaplacesDemon::dlaplace(x = d[z,], scale = .sd[z,], log = TRUE)),
                       FUN.VALUE = numeric(1))
    }
  }
  
  
  return(sum(loglik))
}

cleanLong <- function(currentlong, sub, currentbp){
  
  currentlong[is.na(currentlong)] <- mean(currentlong, na.rm = TRUE)
  current_param <- currentlong[sub]
  
  # force consistency
  currentlong <- rep(current_param, diff(floor(currentbp)))
  
  return(currentlong)
}


loglikeMultivariate <- function(res,
                                along,
                                blong,
                                sigmasq,
                                cormat,
                                timepoints){
  
  loglik <- vector(mode = "numeric", length = length(unique(res$Sequence)))
  centredmu <- matrix(NA, nrow = length(unique(res$Sequence)), ncol = length(timepoints))
  maxUptake <- vector(mode = "numeric", length = length(unique(res$Sequence)))
  
  for (j in seq_along(unique(res$Sequence))){
    
    start <- res$Start[res$Sequence == unique(res$Sequence)[j]][1]
    end <- res$End[res$Sequence == unique(res$Sequence)[j]][1]
    maxUptake[j] <- res$MaxUptake[res$Sequence == unique(res$Sequence)[j]][1]
    seq <- strsplit(unique(res$Sequence)[j], "")[[1]]
    seq <- seq[-c(1:2)]
    index <- which(seq != "P") + start + 1
    ar <- along[index]
    br <- blong[index]
    mu <- rowSums(sapply(seq.int(length(br)),
                         function(z) 
                           logisticFunction(t = timepoints, b = br[z], a = ar[z])))
    centredmu[j, ] <- res$Uptake[res$Sequence == unique(res$Sequence)[j]][seq.int(length(timepoints))] - mu
    
    
  }
  # something like the following, we need to check this works
  # double check eigenvalue decomposition
  SigmaMat <- diag(sqrt(currentSigma * maxUptake)) %*%
    cormat %*%
    diag(sqrt(currentSigma * maxUptake))
  eDecomp <- eigen(SigmaMat)
  sub <- which(eDecomp$values < .Machine$double.eps)
  mu <- diag(1/eDecomp$values[-sub]) %*% muT[-sub,]
  
  loglik <- sum(dnorm(x = mu, mean = 0, sd = 1, log = TRUE))
  
  return(loglik)
}


metropolisSigma <- function(res,
                            currentSigma,
                            pilong,
                            blong,
                            qlong,
                            dlong,
                            mean,
                            sd,
                            meanprop = 0,
                            sdprop = 1,
                            multivariate = FALSE,
                            cormat = NULL,
                            timepoints,
                            tCoef,
                            index,
                            numExch,
                            density = "Gaussian",
                            phi){
  
  # propose new sigma
  proposedSigma <- exp(log(currentSigma) + rnorm(1, mean = meanprop, sd = sdprop))
  
  # compute metropolis ratio
  if (isFALSE(multivariate)){
    numerator <- loglik(res = res,
                        pilong = pilong,
                        blong = blong,
                        qlong = qlong,
                        dlong = dlong,
                        sigmasq = proposedSigma,
                        timepoints = timepoints,
                        tCoef = tCoef, 
                        index = index,
                        numExch = numExch,
                        density = density,
                        phi = phi) + 
      dlnorm(proposedSigma, meanlog = mean, sdlog = sd, log = TRUE)
    denominator <- loglik(res = res,
                          pilong = pilong,
                          blong = blong,
                          qlong = qlong, 
                          dlong = dlong,
                          sigmasq = currentSigma,
                          timepoints = timepoints,
                          tCoef = tCoef, 
                          index = index,
                          numExch = numExch,
                          density = density,
                          phi = phi) + 
      dlnorm(currentSigma, meanlog = mean, sdlog = sd, log = TRUE)
    logratio <- numerator - denominator
  } else{
    numerator <- loglikeMultivariate(res = res,
                                     along = along,
                                     blong = blong,
                                     qlong = qlong,
                                     sigmasq = proposedSigma,
                                     cormat = cormat,
                                     timepoints = timepoints) + 
      dlnorm(proposedSigma, meanlog = mean, sdlog = sd, log = TRUE)
    denominator <- loglikeMultivariate(res = res,
                                       along = along, 
                                       blong = blong,
                                       qlong = qlong,
                                       sigmasq = currentSigma,
                                       cormat = cormat,
                                       timepoints = timepoints) + 
      dlnorm(currentSigma, meanlog = mean, sdlog = sd, log = TRUE)
    logratio <- numerator - denominator
  }
  
  
  # acceptance probability
  u <- log(runif(1))
  if (logratio > u){
    sigma <- proposedSigma
  } else{
    sigma <- currentSigma    
  }
  return(sigma)
}

metropolisRate <- function(res,
                           pilong,
                           qlong, 
                           dlong,
                           currentb,
                           currentSigma,
                           bp,
                           mean = 0,
                           sd = 1,
                           shape = 1,
                           rate = 5,
                           multivariate = FALSE,
                           cormat = NULL,
                           timepoints = timepoints,
                           tCoef,
                           index,
                           numExch,
                           density,
                           phi){
  
  # propose new b
  update <- rnorm(length(currentb), mean = mean, sd = sd)
  randind <- sample(floor(length(currentb)), size = floor(length(currentb))/4)
  proposedb <- currentb
  proposedb[randind] <- currentb[randind] + update[randind] 
  
  for (j in randind){
    tempb <- currentb
    tempb[j] <- proposedb[j]
    proposedblong <- rep(tempb, diff(floor(bp)))
    currentblong <- rep(currentb, diff(floor(bp)))
    
    # compute metropolis ratio
    if (isFALSE(multivariate)){
      numerator <- loglik(res = res,
                          pilong = pilong,
                          blong = proposedblong,
                          qlong = qlong,
                          dlong = dlong,
                          sigmasq = currentSigma,
                          timepoints = timepoints,
                          tCoef = tCoef, 
                          index = index,
                          numExch = numExch,
                          density = density,
                          phi = phi) + 
        dgamma(proposedb[j], shape = shape, rate = rate, log = TRUE)
      denominator <- loglik(res = res,
                            pilong = pilong,
                            blong = currentblong,
                            qlong = qlong,
                            dlong = dlong,
                            sigmasq = currentSigma,
                            timepoints = timepoints,
                            tCoef = tCoef, 
                            index = index,
                            numExch = numExch,
                            density = density,
                            phi = phi) + 
        dgamma(currentb[j], shape = shape, rate = rate, log = TRUE)
      logratio <- numerator - denominator 
    } else {
      numerator <- loglikeMultivariate(res = res,
                                       along = along,
                                       blong = proposedblong,
                                       qlong = qlong,
                                       sigmasq = currentSigma,
                                       cormat = cormat,
                                       timepoints = timepoints) + 
        dgamma(proposedb[j], shape = shape, rate = rate, log = TRUE)
      denominator <- loglikeMultivariate(res = res,
                                         along = along,
                                         blong = currentblong,
                                         qlong = qlong,
                                         sigmasq = currentSigma,
                                         cormat = cormat,
                                         timepoints = timepoints) + 
        dgamma(currentb[j], shape = shape, rate = rate, log = TRUE)
      logratio <- numerator - denominator 
      
    }
    
    
    # acceptance probability
    u <- log(runif(1))
    if (logratio > u){
      currentb <- tempb
    } else{
      currentb <- currentb
    }
  }
  
  return(currentb)
}

metropolisdRate <- function(res,
                            pilong,
                            qlong, 
                            blong,
                            currentd,
                            currentSigma,
                            bp,
                            mean = 0,
                            sd = 1,
                            dshape = 1,
                            drate = 5,
                            multivariate = FALSE,
                            cormat = NULL,
                            timepoints = timepoints,
                            tCoef,
                            index,
                            numExch,
                            density,
                            phi){
  
  # propose new d
  update <- rnorm(length(currentd), mean = mean, sd = sd)
  randind <- sample(floor(length(currentd)), size = floor(length(currentd))/4)
  proposedd <- currentd
  proposedd[randind] <- currentd[randind] + update[randind] 
  
  for (j in randind){
    tempd <- currentd
    tempd[j] <- proposedd[j]
    proposeddlong <- rep(tempd, diff(floor(bp)))
    currentdlong <- rep(currentd, diff(floor(bp)))
    
    # compute metropolis ratio
    if (isFALSE(multivariate)){
      numerator <- loglik(res = res,
                          pilong = pilong,
                          blong = blong,
                          qlong = qlong,
                          dlong = proposeddlong,
                          sigmasq = currentSigma,
                          timepoints = timepoints,
                          tCoef = tCoef, 
                          index = index,
                          numExch = numExch,
                          density = density,
                          phi = phi) + 
        dgamma(proposedd[j], shape = dshape, rate = drate, log = TRUE)
      denominator <- loglik(res = res,
                            pilong = pilong,
                            blong = blong,
                            qlong = qlong,
                            dlong = currentdlong,
                            sigmasq = currentSigma,
                            timepoints = timepoints,
                            tCoef = tCoef, 
                            index = index,
                            numExch = numExch,
                            density = density,
                            phi = phi) + 
        dgamma(currentd[j], shape = dshape, rate = drate, log = TRUE)
      logratio <- numerator - denominator 
    } 
    
    
    # acceptance probability
    u <- log(runif(1))
    if (logratio > u){
      currentd <- tempd
    } else{
      currentd <- currentd
    }
  }
  
  return(currentd)
}

metropolisPlateau <- function(res,
                              currenta,
                              blong,
                              qlong,
                              currentSigma,
                              bp,
                              mean = 0,
                              sd = 1,
                              logmean = 0,
                              logsd = 1,
                              multivariate = FALSE,
                              cormat = NULL,
                              timepoints = timepoints,
                              tCoef,
                              index,
                              numExch,
                              density,
                              phi){
  
  # propose new a
  update <- rnorm(length(currenta), mean = mean, sd = sd)
  randind <- sample(floor(length(currenta)), size = floor(length(currenta))/4)
  proposeda <- currenta
  proposeda[randind] <- exp(log(currenta[randind]) + update[randind]) 
  
  for (j in randind){
    tempa <- currenta
    tempa[j] <- proposeda[j]
    proposedalong <- rep(tempa, diff(floor(bp)))
    currentalong <- rep(currenta, diff(floor(bp)))
    
    # compute metropolis ratio
    if (isFALSE(multivariate)){
      numerator <- loglik(res = res,
                          along = proposedalong,
                          blong = blong,
                          qlong = qlong,
                          sigmasq = proposedSigma,
                          timepoints = timepoints,
                          tCoef = tCoef, 
                          index = index,
                          numExch = numExch,
                          density = density,
                          phi = phi) + 
        dlnorm(proposeda[j], meanlog = logmean, sdlog = logsd, log = TRUE)
      denominator <- loglik(res = res,
                            along = currentalong,
                            blong = blong,
                            qlong = qlong,
                            sigmasq = currentSigma,
                            timepoints = timepoints,
                            tCoef = tCoef,
                            index = index,
                            numExch = numExch,
                            density = density,
                            phi = phi) + 
        dlnorm(currenta[j], meanlog = logmean, sdlog = logsd, log = TRUE)
      logratio <- numerator - denominator  
    } else{
      numerator <- loglikeMultivariate(res = res,
                                       along = proposedalong,
                                       blong = blong,
                                       qlong = qlong,
                                       sigmasq = currentSigma,
                                       cormat = cormat,
                                       timepoints = timepoints) + 
        dlnorm(proposeda[j], meanlog = logmean, sdlog = logsd, log = TRUE)
      denominator <- loglikeMultivariate(res = res,
                                         along = currentalong,
                                         blong = blong,
                                         qlong = qlong,
                                         sigmasq = currentSigma,
                                         cormat = cormat,
                                         timepoints = timepoints) + 
        dlnorm(currenta[j], meanlog = logmean, sdlog = logsd, log = TRUE)
      logratio <- numerator - denominator
    }
    
    
    # acceptance probability
    u <- log(runif(1))
    if (logratio > u){
      currenta <- tempa
    } else{
      currenta <- currenta
    }
  }
  
  return(currenta)
  
}

metropolisPlateau2 <- function(res,
                               currenta,
                               blong,
                               qlong,
                               currentSigma,
                               bp,
                               mean = 0,
                               sd = 1,
                               ashape1 = 1,
                               ashape2 = 1,
                               multivariate = FALSE,
                               cormat = NULL,
                               timepoints = timepoints,
                               tCoef,
                               index,
                               numExch,
                               density,
                               phi){
  
  # propose new a
  update <- rnorm(length(currenta), mean = mean, sd = sd)
  randind <- sample(floor(length(currenta)), size = floor(length(currenta))/4)
  proposeda <- currenta
  proposeda[randind] <- exp(log(currenta[randind]) + update[randind]) 
  
  for (j in randind){
    tempa <- currenta
    tempa[j] <- proposeda[j]
    proposedalong <- rep(tempa, diff(floor(bp)))
    currentalong <- rep(currenta, diff(floor(bp)))
    
    # compute metropolis ratio
    if (isFALSE(multivariate)){
      numerator <- loglik(res = res,
                          along = proposedalong,
                          blong = blong,
                          qlong = qlong,
                          sigmasq = currentSigma,
                          timepoints = timepoints,
                          tCoef = tCoef,
                          index = index,
                          numExch = numExch,
                          density = density,
                          phi = phi) + 
        dbeta(proposeda[j], shape1 = ashape1, shape2 = ashape2, log = TRUE)
      denominator <- loglik(res = res,
                            along = currentalong,
                            blong = blong,
                            qlong = qlong,
                            sigmasq = currentSigma,
                            timepoints = timepoints,
                            tCoef = tCoef,
                            index = index,
                            numExch = numExch,
                            density = density,
                            phi = phi) + 
        dbeta(currenta[j], shape1 = ashape1, shape2 = ashape2, log = TRUE)
      logratio <- numerator - denominator  
    } else{
      numerator <- loglikeMultivariate(res = res,
                                       along = proposedalong,
                                       blong = blong,
                                       qlong = qlong,
                                       sigmasq = currentSigma,
                                       cormat = cormat,
                                       timepoints = timepoints) + 
        dbeta(proposeda[j], shape1 = ashape1, shape2 = ashape2, log = TRUE)
      denominator <- loglikeMultivariate(res = res,
                                         along = currentalong,
                                         blong = blong,
                                         qlong = qlong,
                                         sigmasq = currentSigma,
                                         cormat = cormat,
                                         timepoints = timepoints) + 
        dbeta(currenta[j], shape1 = ashape1, shape2 = ashape2, log = TRUE)
      logratio <- numerator - denominator
    }
    
    
    # acceptance probability
    u <- log(runif(1))
    if (logratio > u){
      currenta <- tempa
    } else{
      currenta <- currenta
    }
  }
  
  return(currenta)
  
}   

metropolisPi <- function(res,
                         currentpi,
                         blong,
                         qlong,
                         dlong,
                         currentSigma,
                         bp,
                         mean = 0,
                         sd = 1,
                         pishape1 = 1,
                         pishape2 = 1,
                         multivariate = FALSE,
                         cormat = NULL,
                         timepoints = timepoints,
                         tCoef,
                         index,
                         numExch,
                         density,
                         phi){
  
  # propose new a
  update <- rnorm(length(currentpi), mean = mean, sd = sd)
  randind <- sample(floor(length(currentpi)), size = floor(length(currentpi))/4)
  proposedpi <- currentpi
  proposedpi[randind] <- exp(log(currentpi[randind]) + update[randind]) 
  
  for (j in randind){
    temppi <- currentpi
    temppi[j] <- proposedpi[j]
    proposedpilong <- rep(temppi, diff(floor(bp)))
    currentpilong <- rep(currentpi, diff(floor(bp)))
    
    # compute metropolis ratio
    if (isFALSE(multivariate)){
      numerator <- loglik(res = res,
                          pilong = proposedpilong,
                          blong = blong,
                          qlong = qlong,
                          dlong = dlong,
                          sigmasq = currentSigma,
                          timepoints = timepoints,
                          tCoef = tCoef,
                          index = index,
                          numExch = numExch,
                          density = density,
                          phi = phi) + 
        dbeta(proposedpi[j], shape1 = pishape1, shape2 = pishape2, log = TRUE)
      denominator <- loglik(res = res,
                            pilong = currentpilong,
                            blong = blong,
                            qlong = qlong,
                            dlong = dlong,
                            sigmasq = currentSigma,
                            timepoints = timepoints,
                            tCoef = tCoef,
                            index = index,
                            numExch = numExch,
                            density = density,
                            phi = phi) + 
        dbeta(currentpi[j], shape1 = pishape1, shape2 = pishape2, log = TRUE)
      logratio <- numerator - denominator  
    }
    
    # acceptance probability
    u <- log(runif(1))
    if (logratio > u){
      currentpi <- temppi
    } else{
      currentpi <- currentpi
    }
  }
  
  return(currentpi)
  
}   


metropolisStretch <- function(res,
                              pilong,
                              blong,
                              dlong,
                              currentq,
                              currentSigma,
                              bp,
                              mean = 0,
                              sd = 1,
                              shape1 = 1,
                              shape2 = 5,
                              multivariate = FALSE,
                              cormat = NULL,
                              timepoints = timepoints,
                              tCoef,
                              index,
                              numExch,
                              density,
                              phi){
  
  # propose new q
  update <- rnorm(length(currentq), mean = mean, sd = sd)
  randind <- sample(floor(length(currentq)), size = floor(length(currentq))/4)
  proposedq <- currentq
  proposedq[randind] <- currentq[randind] + update[randind] 
  
  for (j in randind){
    tempq <- currentq
    tempq[j] <- proposedq[j]
    proposedqlong <- rep(tempq, diff(floor(bp)))
    currentqlong <- rep(currentq, diff(floor(bp)))
    
    # compute metropolis ratio
    if (isFALSE(multivariate)){
      numerator <- loglik(res = res,
                          pilong = pilong,
                          blong = blong,
                          qlong = proposedqlong,
                          dlong = dlong,
                          sigmasq = currentSigma,
                          timepoints = timepoints,
                          tCoef = tCoef,
                          index = index,
                          numExch = numExch,
                          density = density,
                          phi = phi) + 
        dbeta(proposedq[j], shape1 = shape1, shape2 = shape2, log = TRUE)
      denominator <- loglik(res = res,
                            pilong = pilong,
                            blong = blong,
                            qlong = currentqlong,
                            dlong = dlong,
                            sigmasq = currentSigma,
                            timepoints = timepoints,
                            tCoef = tCoef,
                            index = index,
                            numExch = numExch,
                            density = density,
                            phi = phi) + 
        dbeta(currentq[j], shape1 = shape1, shape2 = shape2, log = TRUE)
      logratio <- numerator - denominator 
    } else {
      numerator <- loglikeMultivariate(res = res,
                                       along = along,
                                       blong = blong,
                                       qlong = proposedqlong,
                                       sigmasq = currentSigma,
                                       cormat = cormat,
                                       timepoints = timepoints) + 
        dbeta(proposedq[j], shape1 = shape1, shape2 = shape2, log = TRUE)
      denominator <- loglikeMultivariate(res = res,
                                         along = along,
                                         blong = blong,
                                         qlong = currentqlong,
                                         sigmasq = currentSigma,
                                         cormat = cormat,
                                         timepoints = timepoints) + 
        dbeta(currentq[j], shape1 = shape1, shape2 = shape2, log = TRUE)
      logratio <- numerator - denominator 
      
    }
    
    # acceptance probability
    u <- log(runif(1))
    if (logratio > u){
      currentq <- tempq
    } else{
      currentq <- currentq
    }
  }
  
  return(currentq)
  
}


reversibleJumpBreakPointsAdd <- function(res,
                                         currentSigma,
                                         currentb,
                                         currentpi,
                                         currentbp,
                                         currentq,
                                         currentd,
                                         shape = 1,
                                         rate = 5,
                                         logmean = 0,
                                         logsd = 1,
                                         shape1 = 1,
                                         shape2 = 5,
                                         pishape1 = 1,
                                         pishape2 = 1,
                                         dshape = 1,
                                         drate = 1,
                                         K,
                                         R,
                                         rho,
                                         lambda,
                                         multivariate = FALSE,
                                         cormat = NULL,
                                         timepoints = timepoints,
                                         tCoef,
                                         index,
                                         numExch,
                                         density, 
                                         R_lower,
                                         R_upper,
                                         phi){
  
  j <- sample.int(n = K - 1, size = 1)
  tau <- runif(n = 1, min = currentbp[j], max = currentbp[j + 1])
  
  if ((tau < R_lower) | (tau > R_upper)){
    return(list(currentpi = currentpi,
                currentb = currentb,
                currentq = currentq,
                currentd = currentd,
                currentbp = currentbp, 
                K = K))
    
  }
  
  # breakpoints
  proposedbp <- sort(c(currentbp, tau))
  # new rates and plateau
  additionalb <- rgamma(n = 2, shape = shape, rate = rate)
  #additionala <- rlnorm(n = 2, meanlog = logmean, sdlog = logsd)
  additionalpi <- rbeta(n = 2, shape1 = pishape1, shape2 = pishape2)
  additionalq <- rbeta(n = 2, shape1 = shape1, shape2 = shape2)
  additionald <- rgamma(n = 2, shape = dshape, rate = drate)
  
  
  # new b vector and a vector
  proposedb <- c(currentb[seq_len(j - 1)], additionalb, currentb[-seq_len(j)])
  proposedpi <- c(currentpi[seq_len(j - 1)], additionalpi, currentpi[-seq_len(j)])
  proposedq <- c(currentq[seq_len(j - 1)], additionalq, currentq[-seq_len(j)])
  proposedd <- c(currentd[seq_len(j - 1)], additionald, currentd[-seq_len(j)])
  
  # generated new blong and along
  proposedblong <- rep(proposedb, diff(floor(proposedbp)))
  currentblong <- rep(currentb, diff(floor(currentbp)))
  proposedpilong <- rep(proposedpi, diff(floor(proposedbp)))
  currentpilong <- rep(currentpi, diff(floor(currentbp)))
  proposedqlong <- rep(proposedq, diff(floor(proposedbp)))
  currentqlong <- rep(currentq, diff(floor(currentbp)))
  proposeddlong <- rep(proposedd, diff(floor(proposedbp)))
  currentdlong <- rep(currentd, diff(floor(currentbp)))
  
  # transition probabilities
  logqAdd <- log(rho) + log((1/(K + 1))) + log(R/(currentbp[j + 1] - currentbp[j])) + 
    sum(dgamma(additionalb, shape = shape , rate = rate, log = TRUE)) +
    sum(dbeta(additionalpi, shape1 = pishape1, shape2 = pishape2, log = TRUE)) + 
    sum(dbeta(additionalq, shape1 = shape1, shape2 = shape2, log = TRUE)) + 
    sum(dgamma(additionald, shape = dshape, rate = drate, log = TRUE))
  
  
  #sum(dlnorm(additionala, meanlog = logmean, sdlog = logsd, log = TRUE)) + 
  
  
  logqAddrev <- log(1 - rho) + log((1/(K + 1))) + 
    dgamma(currentb[j], shape = shape, rate = rate, log = TRUE) + 
    dbeta(currentpi[j], shape1 = pishape1, shape2 = pishape2, log = TRUE) + 
    dbeta(currentq[j], shape1 = shape1, shape2 = shape2, log = TRUE) + 
    dgamma(currentd[j], shape = dshape, rate = drate, log = TRUE)
  
  # dlnorm(currenta[j], meanlog = logmean, sdlog = logsd, log = TRUE) + 
  
  
  # factorial term cancels
  if (isFALSE(multivariate)){
    numerator <- loglik(res = res,
                        pilong = proposedpilong,
                        blong = proposedblong,
                        qlong = proposedqlong,
                        sigmasq = currentSigma,
                        dlong = proposeddlong,
                        timepoints = timepoints,
                        tCoef = tCoef,
                        index = index,
                        numExch = numExch,
                        density = density,
                        phi = phi) + 
      sum(dgamma(proposedb, shape = shape, rate = rate, log = TRUE)) + 
      sum(dbeta(proposedpi, shape1 = pishape1, shape2 = pishape2, log = TRUE)) + 
      sum(dbeta(proposedq, shape1 = shape1, shape2 = shape2, log = TRUE)) + 
      sum(dgamma(proposedd, shape = dshape, rate = drate, log = TRUE)) +   
      log((K + 1)/(R - 1)) + 
      dpois(K + 1, lambda = lambda, log = TRUE) + logqAdd
    
    denominator <- loglik(res = res,
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
                          phi = phi) + 
      sum(dgamma(currentb, shape = shape, rate = rate, log = TRUE)) + 
      sum(dbeta(currentpi, shape1 = pishape1, shape2 = pishape2, log = TRUE)) + 
      sum(dbeta(currentq, shape1 = shape1, shape2 = shape2, log = TRUE)) +
      sum(dgamma(currentd, shape = dshape, rate = drate, log = TRUE)) + 
      dpois(K, lambda = lambda, log = TRUE) + logqAddrev
    
    logratio <- numerator - denominator 
  } else {
    print("here")
    numerator <- loglikeMultivariate(res = res,
                                     along = proposedalong,
                                     blong = proposedblong,
                                     qlong = proposedqlong,
                                     sigmasq = currentSigma,
                                     cormat = cormat,
                                     timepoints = timepoints) + 
      sum(dgamma(proposedb, shape = shape, rate = rate, log = TRUE)) + 
      sum(dlnorm(proposeda, meanlog = logmean, sdlog = logsd, log = TRUE)) + 
      sum(dbeta(proposedq, shape1 = shape1, shape2 = shape2, log = TRUE)) + 
      log((K + 1)/(R - 1)) + 
      dpois(K + 1, lambda = lambda, log = TRUE) + logqAdd
    denominator <- loglikeMultivariate(res = res,
                                       along = currentalong,
                                       blong = currentblong,
                                       qlong = currentqlong,
                                       sigmasq = currentSigma,
                                       cormat = cormat,
                                       timepoints = timepoints) + 
      sum(dgamma(currentb, shape = shape, rate = rate, log = TRUE)) + 
      sum(dlnorm(currenta, meanlog = logmean, sdlog = logsd, log = TRUE)) + 
      sum(dbeta(currentq, shape1 = shape1, shape2 = shape2, log = TRUE)) + 
      dpois(K, lambda = lambda, log = TRUE) + logqAddrev
    
    logratio <- numerator - denominator     
  }
  
  
  # acceptance probability
  u <- log(runif(1))
  if (logratio > u){
    currentb <- proposedb
    currentpi <- proposedpi
    currentq <- proposedq
    currentd <- proposedd
    currentbp <- proposedbp
    K <- K + 1
  } else{
    currentb <- currentb
    currentpi <- currentpi
    currentq <- currentq
    currentd <- currentd
    currentbp <- currentbp
    K <- K
  }
  
  return(list(currentpi = currentpi,
              currentb = currentb,
              currentq = currentq,
              currentbp = currentbp,
              currentd = currentd,
              K = K))
}

reversibleJumpBreakPointsRemove <- function(res, 
                                            currentSigma,
                                            currentb,
                                            currentpi,
                                            currentq,
                                            currentd,
                                            currentbp,
                                            logmean = 0,
                                            logsd = 1,
                                            shape = 1,
                                            rate = 5,
                                            shape1 = 1,
                                            shape2 = 5,
                                            pishape1 = 1,
                                            pishape2 = 1,
                                            dshape = 1,
                                            drate = 1,
                                            K,
                                            R,
                                            rho,
                                            lambda,
                                            multivariate = FALSE,
                                            cormat = cormat,
                                            timepoints = timepoints,
                                            tCoef,
                                            index,
                                            numExch,
                                            density,
                                            phi){
  
  # interval to delete, last break point that can be removed is K - 1
  # Can't remove first break point which is 0.
  j <- sample.int(n = K - 2, size = 1) + 1
  
  # breakpoints
  proposedbp <- currentbp[-j]
  # new rates for merged interval and new plateau
  newb <- rgamma(n = 1, shape = shape, rate = rate)
  newpi <- rbeta(n = 1, shape1 = pishape1, shape2 = pishape2)
  newq <- rbeta(n = 1, shape1 = shape1, shape2 = shape2)
  newd <- rgamma(n = 1, shape = dshape, rate = drate)
  
  # new b vector, a and q vector
  proposedb <- c(currentb[seq_len(j - 2)], newb, currentb[-seq_len(j)])
  proposedpi <- c(currentpi[seq_len(j - 2)], newpi, currentpi[-seq_len(j)])
  proposedq <- c(currentq[seq_len(j - 2)], newq, currentq[-seq_len(j)])
  proposedd <- c(currentd[seq_len(j - 2)], newd, currentd[-seq_len(j)])
  
  # generated new blong, along and qlong
  proposedblong <- rep(proposedb, diff(floor(proposedbp)))
  currentblong <- rep(currentb, diff(floor(currentbp)))
  
  proposedpilong <- rep(proposedpi, diff(floor(proposedbp)))
  currentpilong <- rep(currentpi, diff(floor(currentbp)))
  
  proposedqlong <- rep(proposedq, diff(floor(proposedbp)))
  currentqlong <- rep(currentq, diff(floor(currentbp)))
  
  proposeddlong <- rep(proposedd, diff(floor(proposedbp)))
  currentdlong <- rep(currentd, diff(floor(currentbp)))
  
  # transition probabilities removing dimension
  logqremove <- log( 1- rho) + log((1/(K + 1))) + 
    dgamma(newb, shape = shape , rate = rate, log = TRUE) + 
    dbeta(newpi, shape1 = pishape1, shape2 = pishape2, log = TRUE) + 
    dbeta(newq, shape1 = shape1, shape2 = shape2, log = TRUE) + 
    dgamma(newd, shape = dshape, rate = drate, log = TRUE)
  
  logqRemoverev <- log(1 - rho) + log((1/(K))) + log(R/(currentbp[j + 1] - currentbp[j - 1])) + 
    sum(dgamma(currentb[c(j - 1, j)], shape = shape, rate = rate, log = TRUE)) + 
    sum(dbeta(currentpi[c(j - 1, j)], shape1 = pishape1, shape2 = pishape2, log = TRUE)) + 
    sum(dbeta(currentq[c(j - 1, j)], shape1 = shape1, shape2 = shape2, log = TRUE)) + 
    sum(dgamma(currentd[c(j - 1, j)], shape = dshape, rate = drate, log = TRUE))
  
  # factorial term cancels
  if (isFALSE(multivariate)){
    numerator <- loglik(res = res,
                        pilong = proposedpilong,
                        blong = proposedblong,
                        qlong = proposedqlong, 
                        dlong = proposeddlong,
                        sigmasq = currentSigma,
                        timepoints = timepoints,
                        tCoef = tCoef,
                        index = index,
                        numExch = numExch,
                        density = density,
                        phi = phi) + 
      sum(dgamma(proposedb, shape = shape, rate = rate, log = TRUE))  + 
      sum(dbeta(proposedpi, shape1 = pishape1, shape2 = pishape2, log = TRUE)) + 
      sum(dbeta(proposedq, shape1 = shape1, shape2 = shape2, log = TRUE)) +
      sum(dgamma(proposedd, shape = dshape, rate = drate, log = TRUE)) + 
      dpois(K - 1, lambda = lambda, log = TRUE) + logqremove
    
    denominator <- loglik(res = res,
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
                          phi = phi) + 
      sum(dgamma(currentb, shape = shape, rate = rate, log = TRUE)) + log((R - 1)/K) +
      sum(dbeta(currentpi, shape1 = pishape1, shape2 = pishape2, log = TRUE)) +
      sum(dbeta(currentq, shape1 = shape1, shape2 = shape2, log = TRUE)) + 
      sum(dgamma(currentd, shape = dshape, rate = drate, log = TRUE)) + 
      dpois(K, lambda = lambda, log = TRUE) + logqRemoverev
    logratio <- numerator - denominator
  }else{
    numerator <- loglikeMultivariate(res = res,
                                     along = proposedalong,
                                     blong = proposedblong,
                                     qlong = proposedqlong,
                                     sigmasq = currentSigma,
                                     cormat = cormat,
                                     timepoints = timepoints) + 
      sum(dgamma(proposedb, shape = shape, rate = rate, log = TRUE))  + 
      sum(dlnorm(proposeda, meanlog = logmean, sdlog = logsd, log = TRUE)) + 
      sum(dbeta(proposedq, shape1 = shape1, shape2 = shape2, log = TRUE)) +  
      dpois(K - 1, lambda = lambda, log = TRUE) + logqremove
    denominator <- loglikeMultivariate(res = res,
                                       along = currentalong,
                                       blong = currentblong,
                                       qlong = currentqlong,
                                       sigmasq = currentSigma,
                                       cormat = cormat,
                                       timepoints = timepoints) + 
      sum(dgamma(currentb, shape = shape, rate = rate, log = TRUE)) + log((R - 1)/K) +
      sum(dlnorm(currenta, meanlog = logmean, sdlog = logsd, log = TRUE)) +
      sum(dbeta(currentq, shape1 = shape1, shape2 = shape2, log = TRUE)) +
      dpois(K, lambda = lambda, log = TRUE) + logqRemoverev
    logratio <- numerator - denominator
  }
  
  
  # acceptance probability
  u <- log(runif(1))
  if (logratio > u){
    currentb <- proposedb
    currentpi <- proposedpi
    currentq <- proposedq
    currentd <- proposedd
    currentbp <- proposedbp
    K <- K - 1
  } else{
    currentb <- currentb
    currentpi <- currentpi
    currentq <- currentq
    currentd <- currentd
    currentbp <- currentbp
    K <- K
  }
  
  return(list(currentpi = currentpi,
              currentb = currentb,
              currentq = currentq,
              currentd = currentd,
              currentbp = currentbp, 
              K = K))
}