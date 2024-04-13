##' Function that prepares exchange indexes for a given sequence
##' @param res An object of class DataFrame containing the results of the
##'  hdx-ms  experiment
##' @return Returns exchange indexes as a list
##' @md
##'
##' @examples
##' data("BRD4_apo")
##' BRD4_apo <- DataFrame(BRD4_apo)
##'
##' prepareIndexes(res = BRD4_apo)
##'
##' @export
prepareIndexes <- function(res) {
    stopifnot("res must be a DataFrame" = is(res, "DFrame"))


    index <- vector(mode = "list", length = length(unique(res$Sequence)))
    for (j in seq_along(unique(res$Sequence))) {
        start <- res$Start[res$Sequence == unique(res$Sequence)[j]][1]
        end <- res$End[res$Sequence == unique(res$Sequence)[j]][1]
        numExch <- res$MaxUptake[res$Sequence == unique(res$Sequence)[j]][1]
        seq <- strsplit(unique(res$Sequence)[j], "")[[1]]
        seq <- seq[-c(1:2)]
        index[[j]] <- which(seq != "P") + start + 1
    }
    return(index)
}

##' Function that compuate the maximum uptakes for all sequences in the dataset
##' @param res An object of class DataFrame containing the results of the
##' hdx-ms  experiment
##' @return Returns a numeric vector of maximum uptakes
##' @md
##'
##' @examples
##'
##' data("BRD4_apo")
##' BRD4_apo <- DataFrame(BRD4_apo)
##'
##' maxUptakes(res = BRD4_apo)
##'
##' @export
maxUptakes <- function(res) {
  stopifnot("res must be a DataFrame" = is(res, "DFrame"))

    .out <- vapply(seq_along(unique(res$Sequence)),
        function(z) res$MaxUptake[res$Sequence == unique(res$Sequence)[z]][1],
        FUN.VALUE = numeric(1)
    )
    return(.out)
}

##' Plot a coverage heatmap
##'
##' @param res An object of class DataFrame containing the results of the
##' hdx-ms  experiment
##' @param plot A logical value indicating whether to plot the heatmap
##' @return Returns a matrix of coverage
##' @md
##'
##' @examples
##' data("BRD4_apo")
##' BRD4_apo <- DataFrame(BRD4_apo)
##'
##' coverageHeatmap(res = BRD4_apo)
##'
##' @export
coverageHeatmap <- function(res, plot = TRUE) {
  stopifnot("res must be a DataFrame" = is(res, "DFrame"))

    # Get the maximum residue number
    R <- max(res$End)
    # Get the number of unique peptides
    numPeptides <- length(unique(res$Sequence))
    # Prepare the indexes
    index <- prepareIndexes(res = res)

    # Create a matrix to store the coverage
    C <- matrix(0, nrow = R, ncol = numPeptides)

    # Fill the matrix
    for (j in seq_len(ncol(C))) {
        C[index[[j]], j] <- 1
    }

    # Plot the coverage
    if (isTRUE(plot)) {
        image(C,
            col = c("white", "darkblue"),
            xlab = "Residue",
            main = "Peptide Coverage"
        )
    }

    return(C)
}

##' Function that checks a hdx-ms dataset for correct formatting and will
##' clean up the dataset if specified
##'
##'
##'
##'
##' @param res An object of class DataFrame containing the results of the
##'   hdx-ms  experiment
##' @param clean If TRUE will filter multiple charge states and remove missing
##'  values
##' @return Returns a cleaned DataFrame
##'
##' @md
##'
##' @examples
##' data("BRD4_apo")
##' BRD4_apo <- DataFrame(BRD4_apo)
##'
##' cleanHDX(res = BRD4_apo, clean = TRUE)
##'
##'
##'
cleanHDX <- function(res, clean = TRUE) {
  stopifnot("res must be a DataFrame" = is(res, "DFrame"))

    # define the needed columns
    neccessaryColumns <- c("State", "Sequence", "Start", "End", "MaxUptake", "Exposure", "replicate", "Uptake")

    # check if the needed columns are present

    if (all(neccessaryColumns %in% colnames(res)) == FALSE) {
        stop("The dataset is not formatted correctly ensure
         you have the following columns: State, Sequence, Start, End, MaxUptake, Exposure, replicate, Uptake.
         The case is sensitive")
    }

    if (isTRUE(clean)) {
        # check duplications
        duplicated_rows <- duplicated(res[, c(
            "State",
            "Sequence",
            "replicate",
            "Exposure"
        )])

        # remove duplicated rows
        res <- res[-which(duplicated_rows), ]

        # remove missing values
        res_which_complete <- data.frame(res) %>% complete.cases() %>% which()
        res_omitted <- res[-res_which_complete, ]
        missing_seqs <- res_omitted[, "Sequence"] %>% unique()

        res <- res[!res$Sequence %in% missing_seqs, ]
    }

    return(res)
}

##' Function that predict the uptake from the parameters of a Rex model
##'
##' @param params An object of class RexParams containing the parameters of the Rex model
##' @return Returns a DataFrame of the predicted uptake
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
##' rex_example <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                           params = rex_example,
##'                           thin = 1,
##'                           range = 1:4,
##'                           whichChains = c(1))
##'                           
##' uptakePredict(rex_example)                                      
##' @export
uptakePredict <- function(params){
  
  stopifnot("params must be a RexParams object" = is(params, "RexParams"))
  
  # Get the parameters
  pe <- posteriorEstimates(params)
  blong <- pe$blong
  pilong <- pe$pilong
  qlong <- pe$qlong
  dlong <- pe$dlong
  
  # get timepoints
  timepoints <- params@chains@chains[[1]]@timepoints
  phi <- params@chains@chains[[1]]@phi
  
  out <- sapply(seq.int(length(blong)), function(z) {
    doublelogisticFunction(timepoints,
                           b = blong[z],
                           a = phi, q = qlong[z],
                           pi = pilong[z], d = dlong[z]
    )
  })
  
  out_long <- DataFrame(Residue = rep(pe$Residues, each = length(timepoints)),
                        timepoints = rep(timepoints, times = length(pe$Residues)),
                        Uptake = as.vector(out))
  
  return(out_long)
  
}

##' Plot the uptake from the output of a Rex model from the uptakePredict function
##' @param Uptake An object of class DataFrame containing the results of the
##' uptakePredict function
##' @param facet A logical value indicating whether to facet the plot 
##'   (Facets by timepoints)
##' @param values A vector of colours to use in the plot of timepoints
##' @param nrow An integer value indicating the number of rows in the facet.
##' @return Returns a ggplot object
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
##' rex_example <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                           params = rex_example,
##'                           thin = 1,
##'                           range = 1:4,
##'                           whichChains = c(1))
##'                           
##' Uptake <- uptakePredict(rex_example)
##' plotUptake(Uptake)                                      
##' @export
plotUptake <- function(Uptake,
                       facet = FALSE,
                       values = brewer.pal(9, "Dark2"),
                       nrow = 2){
  
  stopifnot("Uptake must be a DataFrame" = is(Uptake, "DFrame"))
  
  gg <- ggplot(Uptake, aes(x = Residue,
                           y = Uptake,
                           group = timepoints,
                           col = factor(timepoints))) +
    geom_line(lwd = 1.1, alpha = 0.5) +
    geom_point() + 
    theme_minimal() +
    labs(x = "Residues", y = "Uptake") + 
    labs(color = "Timepoints") + 
    scale_x_continuous(n.breaks = 10) + 
    scale_y_continuous(n.breaks = 10, limits = c(0, 1)) + 
    scale_color_manual(values = values )
    
  if (isTRUE(facet)) {
    
    gg <- gg + facet_wrap(~timepoints, nrow = nrow)
    
  }
   
  return(gg)
}

##' Function to sample uncertainty of the Uptake
##' 
##' @param params An object of class RexParams containing the parameters of
##'  the Rex model
##' @param method A character value indicating the method to use.
##'  Either "fitted" or "predict". Default is "fitted" which only include
##'  uncertainty in the model parameters and "predict" includes uncertainty
##'  in the model parameters and the error term (e.g. observation level error)
##' @param whichChains An integer value indicating which chain to sample from
##' @param tCoef A numeric vector of coefficients to use in the prediction
##'   
##' @return Returns a list of samples 
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
##' samples <- marginalEffect(rex_example)                           
##' @export
marginalEffect <- function(params,
                           method = "fitted",
                           whichChains = 1,
                           tCoef = NULL){
  
  stopifnot("params must be a RexParams object" = is(params, "RexParams"))
  
  # get timepoints
  timepoints <- params@chains@chains[[1]]@timepoints
  phi <- params@chains@chains[[1]]@phi
  numIter <- params@chains@chains[[1]]@numIter
  Residues <- params@summary@Rex.resolution$Resdiues
  
  out_long <- vector(mode = "list", length = numIter)
  
  if (is.null(tCoef)) {
    tCoef <- c(0, rep(1, length(timepoints) -1))
  }
  
  # Get the parameters
  for (i in seq.int(numIter)){
    
    blong <- params@chains@chains[[whichChains]]@blong[, i]
    pilong <- params@chains@chains[[whichChains]]@pilong[, i]
    qlong <- params@chains@chains[[whichChains]]@qlong[, i]
    dlong <- params@chains@chains[[whichChains]]@dlong[, i]
    sigma <- params@chains@chains[[whichChains]]@Sigma[i]
    .sd <- sqrt(sigma) * tCoef
    
  if (method == "fitted"){
    
    out <- sapply(seq.int(length(blong)), function(z) {
      doublelogisticFunction(timepoints,
                             b = blong[z],
                             a = phi, q = qlong[z],
                             pi = pilong[z], d = dlong[z])
    })
    
  } else if (method == "predict"){
    
    if (.sd[1] == 0){
      
      err <- c(0, rlaplace(n = length(timepoints) - 1,
                           location = 0,
                           scale = .sd[-1]))
      
    } else {
      
      err <- rlaplace(n = length(timepoints),
                      location = 0,
                      scale = .sd)
    }
      
    
    out <- sapply(seq.int(length(blong)), function(z) {
      doublelogisticFunction(timepoints,
                             b = blong[z],
                             a = phi, q = qlong[z],
                             pi = pilong[z], d = dlong[z]) + err
        
    })
    
  }
  
  out_long[[i]]  <- DataFrame(Residue = rep(Residues, each = length(timepoints)),
                              timepoints = rep(timepoints, times = length(Residues)),
                              Uptake = as.vector(out[,Residues]),
                              mcmcIter = rep(i, each = length(Residues) * length(timepoints)))
    
    
  }  
  
  samples <- do.call(rbind, out_long)
  
  return(samples)
  
}

##' Function to plot the uncertainty of the Uptake
##' 
##' @param samples An object of class DataFrame containing the samples of the Uptake
##' @param quantiles A numeric vector of quantiles to calculate. 
##'   Default is c(0.025, 0.05, 0.5, 0.95, 0.975)
##' @param values A vector of colours to use in the plot of timepoints
##' @param nrow An integer value indicating the number of rows in the facet.
##' @return Returns a ggplot object
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
##' rex_example <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                           params = rex_example,
##'                           thin = 1,
##'                           range = 1:4,
##'                           whichChains = c(1))
##' samples <- marginalEffect(rex_example)
##' plotUptakeUncertainty(samples)                           
##' @export
plotUptakeUncertainty <- function(samples,
                                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                  values = brewer.pal(9, "Dark2"),
                                  nrow = 2){
  
  stopifnot("samples must be a DataFrame" = is(samples, "DFrame"))
  
  df <- data.frame(samples) %>% group_by(timepoints, Residue) %>% 
    reframe(Uptake = quantile(Uptake, quantiles),
            quantile = quantiles)
  
  df <- df %>% pivot_wider(names_from = quantile, values_from = Uptake)
  
  gg <- df %>%
    ggplot(aes(x = Residue,
               y = `0.5`,
               group = timepoints,
               col = factor(timepoints),
               fill = factor(timepoints))) +
    geom_line(alpha = 1, lwd = 1.2) +
    theme_minimal() +
    geom_ribbon(aes(x = Residue, ymin = `0.025`, ymax = `0.975`), alpha = 0.2) +
    geom_ribbon(aes(x = Residue, ymin = `0.05`, ymax = `0.95`), alpha = 0.5) +
    labs(x = "Residues", y = "Uptake") +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks = 10) +
    facet_wrap(~timepoints, nrow = nrow) + 
    scale_color_manual(values = values ) + 
    scale_fill_manual(values = values ) + 
    labs(color = "Timepoints", fill = "Timepoints")
  
  return(gg)
  
}























