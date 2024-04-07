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




























