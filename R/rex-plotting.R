##' Rex volcano plots for differential analysis
##' 
##' @param diff_params An object of class RexDifferential
##' @param nrow The number of rows in the facet (to seperate timepoints)
##' @param quantity The quantity to plot either "TRE" or "ARE"
##' @return Returns a ggplot object
##' @md
##' 
##' @examples
##' require(RexMS)
##' require(dplyr)
##' require(ggplot2)
##' 
##' data("BRD4_apo")
##' data("BRD4_ibet")
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' BRD4_ibet <- BRD4_ibet %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##'
##' rex_test <- rex(HdxData = DataFrame(BRD4_ibet),
##'                 numIter = 10,
##'                 R = max(BRD4_ibet$End),
##'                 numtimepoints = numTimepoints,
##'                 timepoints = Timepoints,
##'                 seed = 1L,
##'                 tCoef = c(0, rep(1, 5)),
##'                 BPPARAM = SerialParam())
##' 
##' rex_test <- RexProcess(HdxData = DataFrame(BRD4_ibet),
##'                        params = rex_test,
##'                        range = 5:10,
##'                        thin = 1, whichChains = c(1,2))
##' 
##' rex_diff <- processDifferential(params = rex_test,
##'                                 HdxData = DataFrame(BRD4_apo),
##'                                  whichChain = c(1))
##'
##' gg1 <- plotVolcano(diff_params = rex_diff,
##'                    nrow = 5,
##'                    quantity = "TRE")
##' print(gg1)
##'
##'
##' @export
plotVolcano <- function(diff_params,
                        nrow = 5,
                        quantity = "TRE") {
  
    if (quantity == "TRE") {
        TRE <- diff_params@Rex.estimates[, grep("TRE", colnames(diff_params@Rex.estimates))]
        probs <- diff_params@Rex.probs[, -c(1, ncol(diff_params@Rex.probs))]
        timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(TRE)))



        df_volcanoe <- data.frame(
            TRE = unlist(TRE),
            probs = unlist(probs),
            residues = rep(diff_params@Rex.probs$Residues,
                times = length(timepoints)
            ),
            timepoints = rep(timepoints,
                each = length(diff_params@Rex.probs$Residues)
            )
        )


        gg1 <- df_volcanoe %>%
            ggplot(aes(x = TRE, y = probs, col = TRE < 0)) +
            geom_point(size = 2) +
            theme_bw() +
            facet_wrap(. ~ timepoints, nrow = nrow) +
            xlim(c(-max(abs(df_volcanoe$TRE)), max(abs(df_volcanoe$TRE)))) +
            xlab("Effect Size (TRE)") +
            ylab("Probability") +
            scale_color_manual(
                values = alpha(c("darkred", "darkblue"), 0.5),
                labels = c("Deprotection", "Protection")
            ) +
            labs(color = "") +
            guides(
                alpha = guide_legend(override.aes = list(size = 5)),
                color = guide_legend(override.aes = list(size = 5))
            ) +
            theme(
                text = element_text(size = 20),
                panel.spacing = unit(1, "lines"),
                strip.background = element_rect(fill = "steelblue")
            )
    } else if (quantity == "ARE") {
        ARE <- diff_params@Rex.estimates[, grep("ARE", colnames(diff_params@Rex.estimates))]
        probs <- diff_params@Rex.probs[, -c(1, ncol(diff_params@Rex.probs))]
        timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(ARE)))



        df_volcanoe <- data.frame(
            ARE = unlist(ARE),
            probs = unlist(probs),
            residues = rep(diff_params@Rex.probs$Residues,
                times = length(timepoints)
            ),
            timepoints = rep(timepoints,
                each = length(diff_params@Rex.probs$Residues)
            )
        )

        gg1 <- df_volcanoe %>%
            ggplot(aes(x = ARE, y = probs, col = ARE < 0)) +
            geom_point(size = 2) +
            theme_bw() +
            facet_wrap(. ~ timepoints, nrow = nrow) +
            xlim(c(-max(abs(df_volcanoe$ARE)), max(abs(df_volcanoe$ARE)))) +
            xlab("Effect Size (ARE)") +
            ylab("Probability") +
            scale_color_manual(
                values = alpha(c("darkred", "darkblue"), 0.5),
                labels = c("Deprotection", "Protection")
            ) +
            labs(color = "") +
            guides(
                alpha = guide_legend(override.aes = list(size = 5)),
                color = guide_legend(override.aes = list(size = 5))
            ) +
            theme(
                text = element_text(size = 20),
                panel.spacing = unit(1, "lines"),
                strip.background = element_rect(fill = "steelblue")
            )
    }



    return(gg1)
}

##' Rex butterfly plots for differential analysis
##' 
##' @param diff_params An object of class RexDifferential
##' @param nrow The number of rows in the facet (to seperate timepoints)
##' @param quantity The quantity to plot either "TRE" or "signedARE"
##' @param interval The interval to plot (Residues)
##' @return Returns a ggplot object
##' @md
##' 
##' @examples
##' require(RexMS)
##' require(dplyr)
##' require(ggplot2)
##' 
##' data("BRD4_apo")
##' data("BRD4_ibet")
##' 
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' BRD4_ibet <- BRD4_ibet %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##'
##' rex_test <- rex(HdxData = DataFrame(BRD4_ibet),
##'                 numIter = 10,
##'                 R = max(BRD4_ibet$End),
##'                 numtimepoints = numTimepoints,
##'                 timepoints = Timepoints,
##'                 seed = 1L,
##'                 tCoef = c(0, rep(1, 5)),
##'                 BPPARAM = SerialParam())
##' 
##' rex_test <- RexProcess(HdxData = DataFrame(BRD4_ibet),
##'                        params = rex_test,
##'                        range = 5:10,
##'                        thin = 1, whichChains = c(1,2))
##' 
##' rex_diff <- processDifferential(params = rex_test,
##'                                 HdxData = DataFrame(BRD4_apo),
##'                                 whichChain = c(1))
##'
##' gg1 <- plotButterfly(diff_params = rex_diff,
##'                    nrow = 5,
##'                    quantity = "TRE")
##' print(gg1)
##'
##'
##' @export
plotButterfly <- function(diff_params,
                          nrow = 5,
                          quantity = "TRE",
                          interval = NULL) {
  
  if (is.null(interval)) {
    interval <- seq.int(1, nrow(diff_params@Rex.estimates))
  } else if(max(interval) > max(diff_params@Rex.probs$Residues)){
    interval <- seq.int(min(interval), max(diff_params@Rex.probs$Residues))
  }   
  
  
  
  if (quantity == "TRE") {
    TRE <- diff_params@Rex.estimates[interval, grep("TRE", colnames(diff_params@Rex.estimates))]
    probs <- diff_params@Rex.probs[interval, -c(1, ncol(diff_params@Rex.probs))]
    timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(TRE)))
    
    
    
    df_butterfly <- data.frame(
      TRE = unlist(TRE),
      probs = unlist(probs),
      residues = rep(diff_params@Rex.probs$Residues[interval],
                     times = length(timepoints)
      ),
      timepoints = rep(timepoints,
                       each = length(diff_params@Rex.probs$Residues[interval]))
      )

  
    gg2 <- df_butterfly %>% 
      ggplot(aes(x = residues, y = TRE, col = TRE < 0, group = timepoints)) +
      geom_point(aes(alpha = probs), size = 2) + 
      theme_bw() +
      geom_line() + 
      facet_wrap(.~timepoints, nrow = nrow) + 
      ylim(c(-max(abs(df_butterfly$TRE)), max(abs(df_butterfly$TRE)))) +
      ylab("Effect Size (TRE)") +
      scale_alpha_continuous(range = c(0,1)) + 
      scale_color_manual(values = c("darkred", "darkblue"), 
                         labels = c("Deprotection", "Protection")) +
      labs(color = "", alpha = "Probability") + 
      guides(alpha = guide_legend(override.aes = list(size = 5)),
             color = guide_legend(override.aes = list(size = 5))) + 
      theme(text = element_text(size = 20),
            panel.spacing= unit(1, "lines"),
            strip.background =element_rect(fill="steelblue")) + 
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
    
    } else if(quantity == "signedARE"){
      
      ARE <- diff_params@Rex.estimates[interval, grep("signedARE", colnames(diff_params@Rex.estimates))]
      probs <- diff_params@Rex.probs[interval, -c(1, ncol(diff_params@Rex.probs))]
      timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(ARE)))
      
      
      df_butterfly <- data.frame(
        ARE = unlist(ARE),
        probs = unlist(probs),
        residues = rep(diff_params@Rex.probs$Residues[interval],
                       times = length(timepoints)
        ),
        timepoints = rep(timepoints,
                         each = length(diff_params@Rex.probs$Residues[interval]))
        )
      
      df_butterfly$ARE[is.nan(df_butterfly$ARE)] <- 0
      
      
      gg2 <- df_butterfly %>% 
        ggplot(aes(x = residues, y = ARE, col = ARE < 0, group = timepoints)) +
        geom_point(aes(alpha = probs), size = 2) + 
        theme_bw() +
        geom_line() + 
        facet_wrap(.~timepoints, nrow = nrow) + 
        ylim(c(-max(abs(df_butterfly$ARE)), max(abs(df_butterfly$ARE)))) +
        ylab("Effect Size (ARE)") +
        scale_alpha_continuous(range = c(0,1)) + 
        scale_color_manual(values = c("darkred", "darkblue"), 
                           labels = c("Deprotection", "Protection")) +
        labs(color = "", alpha = "Probability") + 
        guides(alpha = guide_legend(override.aes = list(size = 5)),
               color = guide_legend(override.aes = list(size = 5))) + 
        theme(text = element_text(size = 20),
              panel.spacing= unit(1, "lines"),
              strip.background =element_rect(fill="steelblue")) + 
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
      
      
    }
    
    return(gg2)
    
}

##' Rex time averaged butterfly plots for differential analysis
##' 
##' @param diff_params An object of class RexDifferential
##' @param interval The interval to plot (Residues)
##' @return Returns a ggplot object
##' @md
##' 
##' @examples
##' require(RexMS)
##' require(dplyr)
##' require(ggplot2)
##' 
##' data("BRD4_apo")
##' data("BRD4_ibet")
##' 
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' BRD4_ibet <- BRD4_ibet %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##'
##' rex_test <- rex(HdxData = DataFrame(BRD4_ibet),
##'                 numIter = 10,
##'                 R = max(BRD4_ibet$End),
##'                 numtimepoints = numTimepoints,
##'                 timepoints = Timepoints,
##'                 seed = 1L,
##'                 tCoef = c(0, rep(1, 5)),
##'                 BPPARAM = SerialParam())
##' 
##' rex_test <- RexProcess(HdxData = DataFrame(BRD4_ibet),
##'                        params = rex_test,
##'                        range = 5:10,
##'                        thin = 1, whichChains = c(1,2))
##' 
##' rex_diff <- processDifferential(params = rex_test,
##'                                 HdxData = DataFrame(BRD4_apo),
##'                                 whichChain = c(1))
##'
##' gg1 <- plotTimeAveragedButterfly(diff_params = rex_diff)
##' print(gg1)
##'
##'
##' @export
plotTimeAveragedButterfly <- function(diff_params,
                                      interval = NULL){
  
  if (is.null(interval)) {
    interval <- seq.int(1, nrow(diff_params@Rex.estimates))
  } else if(max(interval) > max(diff_params@Rex.probs$Residues)){
    interval <- seq.int(min(interval), max(diff_params@Rex.probs$Residues))
  }        

  # get important quantities
  TRE <- diff_params@Rex.estimates[interval, grep("TRE", colnames(diff_params@Rex.estimates))]
  totalprobs <- diff_params@Rex.probs[interval, ncol(diff_params@Rex.probs)]
  timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(TRE)))


  df_butterfly2 <- data.frame(
                      TRE = rowMeans(as.matrix(TRE)),
                      totalprobs = totalprobs,
                      residues = diff_params@Rex.probs$Residues[interval])

  gg3 <- df_butterfly2 %>%
    ggplot(aes(x = residues, y = TRE, col = TRE < 0, group = 1)) +
    geom_point(aes(alpha = totalprobs), size = 2) + 
    theme_bw() + 
    geom_line() + 
    ylim(c(-max(abs(df_butterfly2$TRE)), max(abs(df_butterfly2$TRE)))) +
    ylab("Time Averaged Effect Size (TRE)") + 
    scale_alpha_continuous(range = c(0,1)) + 
    scale_color_manual(values = c("darkred", "darkblue"),
                       labels = c("Deprotection", "Protection")) +
    labs(color = "", alpha = "Probability") + 
    guides(alpha = guide_legend(override.aes = list(size = 5)),
           color = guide_legend(override.aes = list(size = 5))) + 
    theme(text = element_text(size = 20),
          panel.spacing = unit(1, "lines"),
          strip.background = element_rect(fill = "steelblue")) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  
  return(gg3)

}
##' Plot Residue level error rates to identify points where modelling is poor
##' 
##' @param rex_params An object of class RexParams
##' @param nrow The number of rows in the facet (to seperate timepoints)
##' @param interval The interval to plot (Residues)
##' @return Returns a ggplot object
##' 
##' @examples
##' require(RexMS)
##' require(dplyr)
##' require(ggplot2)
##' 
##' data("BRD4_apo")
##' 
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##' 
##' rex_test <- rex(HdxData = DataFrame(BRD4_apo),
##'                numIter = 10,
##'                R = max(BRD4_apo$End),
##'                numtimepoints = numTimepoints,
##'                timepoints = Timepoints,
##'                seed = 1L,
##'                tCoef = c(0, rep(1, 5)),
##'                BPPARAM = SerialParam())
##'                
##' rex_test <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                       params = rex_test,
##'                       range = 5:10,
##'                       thin = 1, whichChains = c(1,2))
##' 
##' gg1 <- plotResidueResolution(rex_params = rex_test, nrow = 5)
##'                                                                    
##' print(gg1)
##' 
##' @export  
plotResidueResolution <- function(rex_params, 
                                  nrow = 1,
                                  interval = NULL) {
  
  if (is.null(interval)) {
    interval <- seq.int(1, nrow(rex_params@summary@posteriorEstimates))
  } else if(max(interval) > max(nrow(rex_params@summary@posteriorEstimates))){
    interval <- seq.int(min(interval), max(rex_params@summary@Rex.resolution$Resdiues))
  }     
  
  numTimepoints <- rex_params@chains@chains[[1]]@numTimepoints
  timepoints <- rex_params@chains@chains[[1]]@timepoints
  
  ARE <- Rex.resolution(rex_params)[interval, grep("ARE", colnames(Rex.resolution(rex_params)))]
  ARE <- ARE[,-seq.int(numTimepoints - 1)]
  timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(ARE)))
  
  
  df_butterfly <- data.frame(
    ARE = unlist(ARE),
    residues = rep(Rex.resolution(rex_params)$Resdiues[interval],
                   times = length(timepoints)
    ),
    timepoints = rep(timepoints,
                     each = length(Rex.resolution(rex_params)$Resdiues[interval]))
  )
  
  gg3 <- df_butterfly %>%
    ggplot(aes(x = residues, y = ARE, group = 1)) +
    geom_point(alpha = 0.9, size = 2) + 
    theme_bw() + 
    geom_line() + 
    ylim(0, max(abs(df_butterfly$ARE))*1.1) +
    ylab("ARE") + 
    facet_wrap(.~timepoints, nrow = nrow) +
    scale_alpha_continuous(range = c(0,1)) + 
    scale_color_manual(values = c("darkgreen")) +
    labs(color = "") + 
    guides(alpha = guide_legend(override.aes = list(size = 5)),
           color = guide_legend(override.aes = list(size = 5))) + 
    theme(text = element_text(size = 20),
          panel.spacing = unit(1, "lines"),
          strip.background = element_rect(fill = "steelblue")) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

  return(gg3)

}
##' Plot Peptide level error rates
##' 
##' @param rex_params An object of class RexParams
##' @param HdxData An object of class DataFrame, containing the HDX-MS data used.
##' @param relative Logical, if TRUE the error rates are normalised
##'  by the maximum uptake
##' @return Returns a pheatmap object 
##' 
##' @examples
##' require(RexMS)
##' require(dplyr)
##'
##' data("BRD4_apo")
##' 
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##' 
##' rex_test <- rex(HdxData = DataFrame(BRD4_apo),
##'               numIter = 10,
##'               R = max(BRD4_apo$End),
##'               numtimepoints = numTimepoints,
##'               timepoints = Timepoints,
##'               seed = 1L,
##'               tCoef = c(0, rep(1, 5)),
##'               BPPARAM = SerialParam())
##'               
##' rex_test <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                     params = rex_test,
##'                     range = 5:10,
##'                     thin = 1, whichChains = c(1,2))
##'                     
##' plotPeptideError(rex_params = rex_test, HdxData = DataFrame(BRD4_apo))
##' 
##' @export
plotPeptideError <- function(rex_params,
                             HdxData,
                             relative = FALSE){
  
  
  df <- as.matrix(Rex.peptideError(rex_params)[,-1])
  rownames(df) <- Rex.peptideError(rex_params)[,1] 
  
  if (relative) {
    
    Sequence <- Rex.peptideError(rex_params)[, "Peptide"]
    norm <- str_count(str_sub(Sequence, 3, -1)) - str_count(str_sub(Sequence, 3, -1), "P")
    df <- df / norm
  }
  
  rg <- max(abs(df))
  
  pp <- pheatmap(df, 
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 show_rownames = TRUE,
                 breaks = seq(-rg, rg, length.out = 100))
  
  
  return(pp)
}

##' plot loglikelihoods across chains
##' 
##' @param rex_params An object of class RexParams
##' @return Returns a ggplot object
##' 
##' 
##' @examples
##' require(RexMS)
##' require(dplyr)
##' require(ggplot2)
##' 
##' data("BRD4_apo")
##' 
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##' 
##' rex_test <- rex(HdxData = DataFrame(BRD4_apo),
##'              numIter = 10,
##'              R = max(BRD4_apo$End),
##'              numtimepoints = numTimepoints,
##'              timepoints = Timepoints,
##'              seed = 1L,
##'              tCoef = c(0, rep(1, 5)),
##'              BPPARAM = SerialParam())
##'
##' rex_test <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                    params = rex_test,
##'                    range = 5:10,
##'                    thin = 1, whichChains = c(1,2))
##'                    
##' plotLogLikelihoods(rex_params = rex_test)
##' 
##' @export
plotLogLikelihoods <- function(rex_params){
  
  
  logLikelihoods <- lapply(rex_params@chains@chains, function(x)
                                    x@loglikelihood)
  
  df <- data.frame(
    logLikelihoods = unlist(logLikelihoods),
    chain = as.factor(rep(seq.int(length(rex_params@chains@chains)),
                each = rex_params@chains@chains[[1]]@numIter)),
    iter = rep(seq.int(rex_params@chains@chains[[1]]@numIter),
               times = length(rex_params@chains@chains)))
  
  gg <- ggplot(df, aes(x = iter, y = logLikelihoods, group = chain, col = chain)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line() +
    theme_bw() +
    ylab("log Likelihood") +
    xlab("iteration") +
    theme(text = element_text(size = 20)) + 
    scale_color_manual(values = brewer.pal(n = max(3, length(unique(df$chain))), 
                                           name = "Set2"))
  
  return(gg)
}

##' plot sigmas across chains
##' 
##' @param rex_params An object of class RexParams
##' @return Returns a ggplot object
##' 
##' @examples
##' 
##' require(RexMS)
##' require(dplyr)
##' require(ggplot2)
##' 
##' data("BRD4_apo")
##' 
##' BRD4_apo <- BRD4_apo %>% filter(End < 100)
##' numTimepoints <- length(unique(BRD4_apo$Exposure))
##' Timepoints <- unique(BRD4_apo$Exposure)
##' numPeptides <- length(unique(BRD4_apo$Sequence))
##' 
##' rex_test <- rex(HdxData = DataFrame(BRD4_apo),
##'              numIter = 10,
##'              R = max(BRD4_apo$End),
##'              numtimepoints = numTimepoints,
##'              timepoints = Timepoints,
##'              seed = 1L,
##'              tCoef = c(0, rep(1, 5)),
##'              BPPARAM = SerialParam())
##'
##' rex_test <- RexProcess(HdxData = DataFrame(BRD4_apo),
##'                    params = rex_test,
##'                    range = 5:10,
##'                    thin = 1, whichChains = c(1,2))
##'                    
##' plotSigma(rex_params = rex_test)
plotSigma <- function(rex_params){
  
  sigmas <- lapply(rex_params@chains@chains, function(x)
                     x@Sigma)
  
  df <- data.frame(
    sigmas = unlist(sigmas),
    chain = as.factor(rep(seq.int(length(rex_params@chains@chains)),
                          each = rex_params@chains@chains[[1]]@numIter)),
    iter = rep(seq.int(rex_params@chains@chains[[1]]@numIter),
               times = length(rex_params@chains@chains)))
  
  gg <- ggplot(df, aes(x = iter, y = log(sigmas), group = chain, col = chain)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line() +
    theme_bw() +
    ylab("log sigma") +
    xlab("iteration") +
    theme(text = element_text(size = 20)) + 
    scale_color_manual(values = brewer.pal(n = max(3, length(unique(df$chain))), 
                                           name = "Set2"))
  
  return(gg)
}

