##' Funcition for performing unsupervised conformational signature analysis. 
##' This performs PCA on the TRE values of the RexDifferentialList and 
##' returns the PCA object and the TRE values in wide format.
##' 
##' @param RexDifferentialList A list of RexDifferential objects
##' @param quantity The quantity to use for the analysis. Default is "TRE"
##' @param states The state name to use for the analysis/ e.g. ligand used in 
##'  differential analysis
##' @param whichTimepoint The timepoint to use for the analysis. Default is 600
##' @param pca_params The parameters to use for the PCA. 
##' Default is list(scale = FALSE, center = TRUE)
##' 
##' 
##' 
##' @return A list containing the PCA object and the TRE values in wide format
##' 
##' @examples
##' library("RexMS")
##' data("out_lxr_compound_proccessed")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds)
##' 
##' 
##' ucsa <- UnsupervisedCSA(out_lxr_compound_proccessed,
##'                        quantity = "TRE",
##'                        states = states,
##'                        whichTimepoint = 600,
##'                        pca_params = list(scale = FALSE,
##'                                          center = TRUE))
##'                        
##' 
##' 
##' @export
UnsupervisedCSA <- function(RexDifferentialList,
                            quantity = "TRE",
                            states,
                            whichTimepoint = 600,
                            pca_params = list(scale = FALSE,
                                              center = TRUE)){
  
  TRE <- lapply(RexDifferentialList, function(diff_params) 
    diff_params@Rex.estimates[,
                              grep(quantity,
                                   colnames(diff_params@Rex.estimates))])
  
  probs <- lapply(RexDifferentialList, function(diff_params)
    diff_params@Rex.probs[, -c(1, ncol(diff_params@Rex.probs))])
  
  timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(TRE[[1]])))
  
  residue_list <- lapply(RexDifferentialList, function(x) x@Rex.estimates$Resdiues)
  
  df_TRE <- DataFrame(TRE = do.call(rbind, TRE),
                      probs = do.call(rbind, probs),
                      Residues = unlist(residue_list),
                      states = rep(states, times = sapply(residue_list, length)))
  
  values_from <- grep(paste0(quantity, "_", whichTimepoint), colnames(df_TRE))[1]
  
  TRE_states_wide <- data.frame(df_TRE) %>% 
    pivot_wider(names_from = states,
                values_from = all_of(values_from),
                id_cols = "Residues")
  
  quant_df <- t(TRE_states_wide[, -1])
  quant_df_full <- quant_df[,colSums(is.na(quant_df)) == 0 ]
  
  pca_states <- prcomp(quant_df_full,
                       scale = pca_params$scale,
                       center = pca_params$center)
  
  
  return(list(pca_states = pca_states, states_wide = TRE_states_wide))
}

##' Function to plot the PCA results from the UnsupervisedCSA function
##' 
##' @param pca_states The PCA object from the UnsupervisedCSA function
##' @param states_wide The values in wide format from the UnsupervisedCSA function
##' @param labels The labels to use for the plot. Construct labels carefully
##' using example below. The labels should be a data frame with the rownames
##' as the states and the columns as the labels to use for the colouring of the
##' labels.
##' @param pca_params The parameters to use for the PCA plot. 
##'   Default is list(x = 1, y = 2, whichlabel = NULL). x indicates principal
##'   component to plot on x-axis, y indicates principal component to plot on
##'    y-axis. whichlabel indicates the column in the labels to use for the
##'    colouring of the labels.
##' @param values The values to use for the colouring of the labels. Default is
##' `brewer.pal(n = 8, name = "OrRd")[c(8,3)]` and "grey" (for unknown labels)
##'    
##' @return A ggplot object
##' 
##' @examples
##' library("RexMS")
##' library(ggfortify)
##' # Construct labels carefully using known properties of the states (Ligands)
##' 
##' data("out_lxr_compound_proccessed")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds) 
##'
##' labels <- data.frame(ABCA1 = rep("Unknown", length(states)),
##' lipogenic = rep("Unknown", length(states)))
##' rownames(labels) <- states
##' 
##' labels$ABCA1[rownames(labels) %in% c("LXR.623", "AZ9", "AZ8", "AZ5")] <- "low"
##' labels$ABCA1[rownames(labels) %in% c("Az1", "AZ2", "AZ3", "AZ4", "AZ6",
##'                                     "AZ7", "AZ876", "T0.901317", "WAY.254011",
##'                                     "F1", "GW3965", "BMS.852927")] <- "high"
##' 
##' 
##' labels$lipogenic[rownames(labels) %in% c("AZ6", "AZ7", "AZ9",
##'                                         "AZ8", "GW3965", "BMS.852927",
##'                                         "LXR.623")] <- "Non-Lipogenic"
##'
##' labels$lipogenic[rownames(labels) %in% c("AZ876", "AZ1",
##'                                         "T0.901317", "F1", "WAY.254011")] <- "Lipogenic"
##'
##' labels$ABCA1 <- factor(labels$ABCA1,
##'                        levels = c("low", "high", "Unknown"))
##' labels$lipogenic <- factor(labels$lipogenic,
##'                           levels = c("Non-Lipogenic", "Lipogenic", "Unknown"))
##' 
##' 
##' 
##' ucsa <- UnsupervisedCSA(out_lxr_compound_proccessed,
##'                        quantity = "TRE",
##'                        states = states,
##'                        whichTimepoint = 600,
##'                        pca_params = list(scale = FALSE,
##'                                          center = TRUE))
##'                        
##'                        
##' plotUCSA(pca_states = ucsa$pca_states,
##'          states_wide = ucsa$states_wide,
##'          labels = labels)
##'
##'
##' @export
plotUCSA <- function(pca_states,
                     states_wide,
                     labels,
                     pca_params = list(x = 1,
                                       y = 2,
                                       whichlabel = NULL),
                     values = c(brewer.pal(n = 8,
                                           name = "OrRd")[c(8,3)],
                                "grey")){
  
  if (is.null(pca_params$whichlabel)) {
    colours_labels <- factor("Unknown", levels = c("Unknown"))
  } else {
    colours_labels <- labels[,pca_params$whichlabel]
  }
  
  
  hh1 <- autoplot(pca_states,
                  loadings = TRUE,
                  shape = TRUE,
                  data = t(states_wide[,-1]),
                  loadings.colour = 'blue',
                  loadings.label.size = 3, x = pca_params$x, y = pca_params$y) +
    geom_point(pch = 19, size = 0)  + 
    theme_bw() + 
    geom_text(aes(label = rownames(labels),
                  col = colours_labels,
                  size = 4,
                  fontface = "bold")) + 
    scale_color_manual(values = values, "grey") + 
    theme(legend.position = "none") + 
    theme(text = element_text(size = 14))
  return(hh1)
  
}

##' Function to perform supervised conformational signature analysis. In this
##' case the construction of the signatures using labels as part of the 
##' dimensionality reduction which are defined in labels. The type of the labels
##' can be catagorical or continuous. The function returns the OPLS-DA object.
##' The labels are allowed to contain "Unknown" values which are ignored in the
##' analysis if the type is catagorical. If the type is continuous, the "Unknown"
##' values should be recorded as NAs.
##' 
##' @param RexDifferentialList A list of RexDifferential objects
##' @param quantity The quantity to use for the analysis. Default is "TRE"
##' @param states The state name to use for the analysis. e.g. ligand used in
##' differential analysis
##' @param labels The labels to use for the analysis. Construct labels carefully
##' using example below. The labels should be a data frame with the rownames
##' as the states and the columns as the labels to use for the analysis.
##' @param whichlabel The column in the labels to use for the analysis. Default
##' is 1
##' @param whichTimepoint The timepoint to use for the analysis. Default is 600
##' @param type The type of the labels. Default is "catagorical" but could be
##' "continuous". If "continuous" the "Unknown" values should be recorded as NAs.
##' @param orthoI The number of orthogonal components to use in the OPLS-DA
##' analysis. Default is 1
##' 
##' @return An OPLS-DA object
##' 
##' @examples
##' 
##' # Construct labels carefully using known properties of the states (Ligands)
##' # first a catagorical example
##' library("RexMS")
##' data("out_lxr_compound_proccessed")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds) 
##'
##' labels <- data.frame(ABCA1 = rep("Unknown", length(states)),
##' lipogenic = rep("Unknown", length(states)))
##' rownames(labels) <- states
##' 
##' labels$ABCA1[rownames(labels) %in% c("LXR.623", "AZ9", "AZ8", "AZ5")] <- "low"
##' labels$ABCA1[rownames(labels) %in% c("Az1", "AZ2", "AZ3", "AZ4", "AZ6",
##'                                     "AZ7", "AZ876", "T0.901317", "WAY.254011",
##'                                     "F1", "GW3965", "BMS.852927")] <- "high"
##' 
##' 
##' labels$lipogenic[rownames(labels) %in% c("AZ6", "AZ7", "AZ9",
##'                                         "AZ8", "GW3965", "BMS.852927",
##'                                         "LXR.623")] <- "Non-Lipogenic"
##'
##' labels$lipogenic[rownames(labels) %in% c("AZ876", "AZ1",
##'                                         "T0.901317", "F1", "WAY.254011")] <- "Lipogenic"
##'
##' labels$ABCA1 <- factor(labels$ABCA1,
##'                        levels = c("low", "high", "Unknown"))
##' labels$lipogenic <- factor(labels$lipogenic,
##'                           levels = c("Non-Lipogenic", "Lipogenic", "Unknown"))
##'
##' # First using ABCA1 as an example
##' 
##' scsa <- supervisedCSA(RexDifferentialList = out_lxr_compound_proccessed,
##'                         quantity = "TRE",
##'                         states = states,
##'                         labels = labels,
##'                         whichlabel = "ABCA1",
##'                         whichTimepoint = 600,
##'                         orthoI = 1)
##'                         
##' # Now using lipogenic as an example
##' 
##' scsa2 <- supervisedCSA(RexDifferentialList = out_lxr_compound_proccessed,
##'                        quantity = "TRE",
##'                        states = states,
##'                        labels = labels,
##'                        whichlabel = "lipogenic",
##'                        whichTimepoint = 600,
##'                        orthoI = 1)
##'
##' # Now using a continuous example, add additional annotation to the labels
##' # Here we use the ED50 values of the ligands, using the log values to 
##' # because of the large range
##' 
##' labels$ED50 <- NA
##' 
##' labels[, "ED50"] <- log(c(4.11, NA, NA, 0.956, NA, 9.64, 1.49,
##'                          5.65, 0.969, 2.10, 11.3, NA, 31.5, 341, 32.2, 17.2,
##'                           NA))
##'                           
##' scsa3 <- supervisedCSA(RexDifferentialList = out_lxr_compound_proccessed,
##'                        quantity = "TRE",
##'                        states = states,
##'                        labels = labels,
##'                        whichlabel = "ED50",
##'                        whichTimepoint = 600,
##'                        type = "continuous",
##'                        orthoI = 1)             
##'
##' @export
supervisedCSA <- function(RexDifferentialList,
                          quantity = "TRE",
                          states,
                          labels,
                          whichlabel = 1,
                          whichTimepoint = 600,
                          type = "catagorical",
                          orthoI = 1){
  
  TRE <- lapply(RexDifferentialList, function(diff_params) 
    diff_params@Rex.estimates[,
                              grep(quantity,
                                   colnames(diff_params@Rex.estimates))])
  
  probs <- lapply(RexDifferentialList, function(diff_params)
    diff_params@Rex.probs[, -c(1, ncol(diff_params@Rex.probs))])
  
  timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(TRE[[1]])))
  
  residue_list <- lapply(RexDifferentialList, function(x) x@Rex.estimates$Resdiues)
  
  df_TRE <- DataFrame(TRE = do.call(rbind, TRE),
                      probs = do.call(rbind, probs),
                      Residues = unlist(residue_list),
                      states = rep(states, times = sapply(residue_list, length)))
  
  df_opls <- cbind(df_TRE, labels[df_TRE$states, ])
  
  values_from <- grep(paste0(quantity, "_", whichTimepoint), colnames(df_opls))[1]
  
  opls_states_wide <- data.frame(df_opls) %>% 
    pivot_wider(names_from = states,
                values_from = all_of(c(values_from)),
                id_cols = "Residues")
  
  if (type == "catagorical") {
    
    col_subset <-  seq.int(nrow(labels))[labels[, whichlabel] != "Unknown"] + 1
    df_reduced <- opls_states_wide[,
                                   col_subset]
    annotations <- labels[labels[, whichlabel] != "Unknown", whichlabel]
  } else {
    
    col_subset <-  seq.int(nrow(labels))[!is.na(as.numeric(labels[, whichlabel]))] + 1
    df_reduced <- opls_states_wide[,
                                   col_subset]
    annotations <- as.numeric(labels[!is.na(as.numeric(labels[, whichlabel])),
                                     whichlabel])
    
  }
  
  cross_val <- min(7, ncol(df_reduced)/2)
  
  if (type == "catagorical") {
    
    df_reduced <- data.frame(df_reduced)
    rownames(df_reduced) <- opls_states_wide$Residues
    df <- t(df_reduced)
    df_na_remove <- df[, colSums(is.na(df)) == 0]
    
    states.plsda <- opls(x = df_na_remove,
                         y =  factor(as.numeric(annotations)),
                         crossvalI = cross_val,
                         orthoI = orthoI)
  } else {
    
    df_reduced <- data.frame(df_reduced)
    rownames(df_reduced) <- opls_states_wide$Residues
    df <- t(df_reduced)
    df_na_remove <- df[, colSums(is.na(df)) == 0]
    
    states.plsda <- opls(x = df_na_remove,
                         y =  as.numeric(annotations),
                         crossvalI = cross_val,
                         orthoI = orthoI)
  }
  
  return(states.plsda = states.plsda)
}

##' Function to plot the supervised conformational signature analysis results
##' 
##' @param states.plsda The OPLS-DA object from the supervisedCSA function
##' @param labels The labels to use for the plot. Construct labels carefully
##' using example below. The labels should be a data frame with the rownames
##' as the states and the columns as the labels to use for the colouring of the
##' labels.
##' @param whichlabel The column in the labels to use for the analysis. Default
##' is 1
##' @param values The values to use for the colouring of the labels. Default is
##' `brewer.pal(n = 8, name = "PRGn")[c(1,8)]` and "grey" (for unknown labels)
##' @param type The type of the labels. Default is "catagorical" but could be
##' "continuous". If "continuous" the "Unknown" values should be recorded as NAs.
##' 
##' @return A ggplot object
##' 
##' @examples
##' 
##' # Construct labels carefully using known properties of the states (Ligands)
##' # first a catagorical example
##' library("RexMS")
##' data("out_lxr_compound_proccessed")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds) 
##'
##' labels <- data.frame(ABCA1 = rep("Unknown", length(states)),
##' lipogenic = rep("Unknown", length(states)))
##' rownames(labels) <- states
##' 
##' labels$ABCA1[rownames(labels) %in% c("LXR.623", "AZ9", "AZ8", "AZ5")] <- "low"
##' labels$ABCA1[rownames(labels) %in% c("Az1", "AZ2", "AZ3", "AZ4", "AZ6",
##'                                     "AZ7", "AZ876", "T0.901317", "WAY.254011",
##'                                     "F1", "GW3965", "BMS.852927")] <- "high"
##' 
##' 
##' labels$lipogenic[rownames(labels) %in% c("AZ6", "AZ7", "AZ9",
##'                                         "AZ8", "GW3965", "BMS.852927",
##'                                         "LXR.623")] <- "Non-Lipogenic"
##'
##' labels$lipogenic[rownames(labels) %in% c("AZ876", "AZ1",
##'                                         "T0.901317", "F1", "WAY.254011")] <- "Lipogenic"
##'
##' labels$ABCA1 <- factor(labels$ABCA1,
##'                        levels = c("low", "high", "Unknown"))
##' labels$lipogenic <- factor(labels$lipogenic,
##'                           levels = c("Non-Lipogenic", "Lipogenic", "Unknown"))
##'
##' # First using ABCA1 as an example
##' 
##' scsa <- supervisedCSA(RexDifferentialList = out_lxr_compound_proccessed,
##'                         quantity = "TRE",
##'                         states = states,
##'                         labels = labels,
##'                         whichlabel = "ABCA1",
##'                         whichTimepoint = 600,
##'                         orthoI = 1)
##' 
##' plotSCSA(states.plsda = scsa,
##'         labels = labels,
##'         whichlabel = "ABCA1")
##'
##' @export                                   
plotSCSA <- function(states.plsda,
                     labels,
                     whichlabel = 1,
                     values = c(brewer.pal(n = 8,
                                           name = "PRGn")[c(1,8)],
                                "grey"),
                     type = "catagorical"){
  
  if (is.null(whichlabel)) {
    colours_labels <- factor("Unknown", levels = c("Unknown"))
  } else {
    colours_labels <- labels[,whichlabel]
  }
  
  if (type == "catagorical"){
    df_pls2 <- as.data.frame(cbind(states.plsda@scoreMN,
                                   states.plsda@orthoScoreMN))
    df_pls2$annotation <- labels[rownames(df_pls2), whichlabel]
    df_pls2$names <- rownames(df_pls2)
    
    
    gg1 <- df_pls2 %>% ggplot(aes(x = p1, y = o1, col = annotation)) +
      theme_bw() + geom_point(size = 0) + 
      scale_color_manual(values = values,
                         labels = levels(df_pls2$annotation)) +
      labs(col = "annotation") +
      theme(text = element_text(size = 20)) + 
      geom_text(aes(label = names), size = 4, fontface = "bold") +  
      labs(x = paste0("Predictive Dimension 1, Variance explained:",
                      states.plsda@modelDF[1,1] * 100)) + 
      labs(y = paste0("Orthogonal Dimension 1, Variance explained:",
                      states.plsda@modelDF[2,1] * 100)) 
    
  } else{
    df_pls2 <- as.data.frame(cbind(states.plsda@scoreMN,
                                   states.plsda@orthoScoreMN))
    df_pls2$annotation <- labels[rownames(df_pls2), whichlabel]
    df_pls2$names <- rownames(df_pls2)
    
    gg1 <- df_pls2 %>% ggplot(aes(x = p1, y = o1, col = annotation)) +
      theme_bw() + geom_point(size = 0) + 
      scale_color_continuous(type = "viridis") +
      labs(col = "annotation") +
      theme(text = element_text(size = 20)) + 
      geom_text(aes(label = names), size = 4, fontface = "bold") +  
      labs(x = paste0("Predictive Dimension 1, Variance explained:",
                      states.plsda@modelDF[1,1] * 100)) + 
      labs(y = paste0("Orthogonal Dimension 1, Variance explained:",
                      states.plsda@modelDF[2,1] * 100))
    
    
    
  }

  return(gg1)
}

##' Function to plot the loading of the supervised conformational signature
##' analysis results
##' 
##' @param states.plsda The OPLS-DA object from the supervisedCSA function
##' @param labels The labels to use for the plot. Construct labels carefully
##' using example below. The labels should be a data frame with the rownames
##' as the states and the columns as the labels to use for the colouring of the
##' labels.
##' @param whichlabel The column in the labels to use for the analysis. Default
##' is 1
##' @param whichLoading The type of loading to plot. Default is "predictive" but
##' could be "orthogonal"
##' @param threshold The threshold to use for the critical residues. Default is
##' 0.125.
##' 
##' @return A ggplot object
##' 
##' @examples
##' # Construct labels carefully using known properties of the states (Ligands)
##' # first a catagorical example
##' library("RexMS")
##' data("out_lxr_compound_proccessed")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds) 
##'
##' labels <- data.frame(ABCA1 = rep("Unknown", length(states)),
##' lipogenic = rep("Unknown", length(states)))
##' rownames(labels) <- states
##' 
##' labels$ABCA1[rownames(labels) %in% c("LXR.623", "AZ9", "AZ8", "AZ5")] <- "low"
##' labels$ABCA1[rownames(labels) %in% c("Az1", "AZ2", "AZ3", "AZ4", "AZ6",
##'                                     "AZ7", "AZ876", "T0.901317", "WAY.254011",
##'                                     "F1", "GW3965", "BMS.852927")] <- "high"
##' 
##' 
##' labels$lipogenic[rownames(labels) %in% c("AZ6", "AZ7", "AZ9",
##'                                         "AZ8", "GW3965", "BMS.852927",
##'                                         "LXR.623")] <- "Non-Lipogenic"
##'
##' labels$lipogenic[rownames(labels) %in% c("AZ876", "AZ1",
##'                                         "T0.901317", "F1", "WAY.254011")] <- "Lipogenic"
##'
##' labels$ABCA1 <- factor(labels$ABCA1,
##'                        levels = c("low", "high", "Unknown"))
##' labels$lipogenic <- factor(labels$lipogenic,
##'                           levels = c("Non-Lipogenic", "Lipogenic", "Unknown"))
##'
##' # First using ABCA1 as an example
##' 
##' scsa <- supervisedCSA(RexDifferentialList = out_lxr_compound_proccessed,
##'                         quantity = "TRE",
##'                         states = states,
##'                         labels = labels,
##'                         whichlabel = "ABCA1",
##'                         whichTimepoint = 600,
##'                         orthoI = 1)
##' 
##' plotLoadingSCSA(states.plsda = scsa,
##'                labels = labels,
##'                whichlabel = "ABCA1",
##'                whichLoading = "predictive",
##'                threshold = 0.125)
##'                
##' @export                                         
plotLoadingSCSA <- function(states.plsda,
                            labels,
                            whichlabel = 1,
                            whichLoading = "predictive",
                            threshold = 0.125){
  
  df_loadings <- data.frame(loadings1 = states.plsda@loadingMN[,1], 
                            loadings2 = states.plsda@orthoLoadingMN[,1],
                            residues = as.numeric(rownames(states.plsda@loadingMN)))
  
  #df_loadings$residues <- as.numeric(sub(".", "", df_loadings$residues))
  
  if (whichLoading == "predictive") {
    df_loadings$loadings1 <- states.plsda@loadingMN[,1]
  } else {
    df_loadings$loadings1 <- states.plsda@orthoLoadingMN[,1]
  }
  
  gg2 <- df_loadings %>% ggplot(aes(x = residues,
                                    xend = residues, y = 0, 
                                    yend = loadings1,
                                    col = abs(loadings1) > threshold)) + 
    geom_point() + 
    theme_bw() + 
    theme(text = element_text(size = 20)) + 
    geom_segment() + 
    ylab("loading") + 
    labs(col = "Critical Residues") + 
    ggtitle(paste0("Loadings of ", whichLoading ," dimension 1"))
  
  return(gg2)
}

##' Function to calculate the uncertainty in the TRE values by sampling the
##' posterior distribution of the Rex parameters. Typically used to obtain
##' uncertainty in the conformational signatures.
##' 
##' @param HdxData The HdxData object
##' @param RexParamsList A list of RexParams objects
##' @param quantity The quantity to use for the analysis. Default is "TRE"
##' @param states The state name to use for the analysis. e.g. ligand used in
##'  differential analysis
##' @param whichChain The chain to use for the analysis. Default is 1
##' @param whichSamples The samples to use for the analysis. 
##'  Default is seq.int(1, 50)
##' @param whichTimepoint The timepoint to use for the analysis. Default is 600. 
##' @param num_montecarlo The number of montecarlo samples to use for the
##'  analysis
##'  
##' @return A data frame containing the TRE values and the residue number over
##' monte carlo interations
##'  
##'   
##' @examples
##' 
##'   
##' \dontrun{
##' library("RexMS")
##' data("LXRalpha_processed")
##' data("rex_lxr")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds)
##' 
##' TRE_dist <- sampleTREuncertainty(HdxData = LXRalpha_processed,
##'                                RexParamsList = rex_lxr,
##'                                quantity = "TRE",
##'                                states = states,
##'                                whichChain = 1,
##'                                whichSamples = seq.int(1, 50),
##'                                whichTimepoint = 600,
##'                                num_montecarlo = 5000)
##'                                
##' 
##' 
##' }
##'
##' @export
sampleTREuncertainty <- function(HdxData,
                                 RexParamsList,
                                 quantity = "TRE",
                                 states,
                                 whichChain = 1,
                                 whichSamples = seq.int(1, 50),
                                 whichTimepoint = 600,
                                 num_montecarlo = 5000){
  
  # process uncertainty
  out_sample_uncertainty <- lapply(RexParamsList, function(x)
    processTREuncertainty(HdxData = DataFrame(HdxData),
                          params = x,
                          whichChain = whichChain,
                          whichSamples = whichSamples,
                          num_montecarlo = num_montecarlo))
  R <- max(HdxData$End)
  
  df <- vector(mode = "list", length = length(RexParamsList))
  for (i in seq_along(out_sample_uncertainty)){
    
    df[[i]] <- data.frame(TRE = c(t(out_sample_uncertainty[[i]][,2,])),
                          residue = rep(seq.int(R), times = length(whichSamples)),
                          mc_num = rep(whichSamples, each = R))
    
  }
  
  df_all <- do.call(rbind, df)
  df_all$state <- rep(states, each = R * length(whichSamples))
  
  return(df_all)
  
}  

##' Function to plot the uncertainty in the TRE values for different states 
##'   e.g. ligands. Typically used to obtain uncertainty in the conformational
##'   landscape.
##'
##' @param df_all The data frame containing the TRE values and the residue number
##' over monte carlo interations. Results from the sampleTREuncertainty function.
##' @param pca_states The PCA object from the UnsupervisedCSA function
##' @param states The state name to use for the analysis. e.g. ligand used in
##' differential analysis
##' @param whichSamples The samples to use for the analysis.
##' Default is seq.int(1, 50)
##' 
##' @return A ggplot object
##' 
##' @examples
##' library("RexMS")
##' data(TRE_dist)
##' data("out_lxr_compound_proccessed")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds)
##' 
##' 
##' ucsa <- UnsupervisedCSA(out_lxr_compound_proccessed,
##'                        quantity = "TRE",
##'                        states = states,
##'                        whichTimepoint = 600,
##'                        pca_params = list(scale = FALSE,
##'                                          center = TRUE))
##'                        
##'                        
##' plotTREuncertainty(df_all = TRE_dist,
##'                   pca_states = ucsa$pca_states,
##'                   states = states,
##'                   whichSamples = seq.int(1, 50))
##' 
##' @export                                            
plotTREuncertainty <- function(df_all,
                               pca_states,
                               states,
                               whichSamples = seq.int(1, 50)){
  ## Create coordinates
  d <- 2
  dims <- c(1,2)
  coords <- matrix(NA, nrow = length(states), ncol = 2)
  res0 <- pca_states
  R <- max(df_all$residue)
  #eigs <- colnames(.pca)
  
  ## Create inital dataset, computing average location in PCA plot
  Y0 <- matrix(NA, nrow = length(states), ncol = 2)
  for (j in seq.int(length(states))) {
    Y0[j, ] <- res0$x[j, seq_len(d)]
  }
  
  ## Repeat for different monte-carlo samples
  Y.lst <- list()
  res.rot <- list()
  for ( i in seq.int(max(df_all$mc_num))) {
    
    data <- matrix(df_all$TRE[df_all$mc_num == unique(df_all$mc_num)[i]], ncol = R, byrow = TRUE)
    res <- prcomp(data, center = TRUE)
    res.proc <- vegan::procrustes(Y0, res$x[, dims])
    Y.df <- data.frame(res.proc$Yrot, states = states)
    Y.lst[[i]] <- Y.df
    res.rot[[i]] <- res$rotation[, 1]
  }
  
  ## create long data formats
  names(Y.lst) <- seq_len(length(whichSamples))
  Y.lst.df <- plyr::ldply(Y.lst, .fun = function(x) x, .id = "mcmcIter")
  cols <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(length(states))
  eigs <- round(res0$sdev[c(1,2)]^2/sum(res0$sdev^2) * 100, 2)
  
  ## ggplot
  gg <- ggplot(
    data = dplyr::mutate(Y.lst.df, states = factor(states)),
    aes(x = X1, y = X2, color = states)) +
    coord_fixed() +
    geom_density2d(contour = TRUE) +
    geom_point(alpha = 0.7) +
    xlab(paste0(eigs[1], "%")) +
    ylab(paste0(eigs[2], "%")) +
    theme(legend.position = "right",
          text = element_text(size = 12)) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 0.5,
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 14),
          text = element_text(size = 20)) +
    ggtitle(label = "Uncertainty in conformational landscape")
  
  return(gg)
}  
##' Function to plot the uncertainty in the loading in the conformational 
##' landscape. Typically used to obtain uncertainty in the conformational 
##' signatures.
##' 
##' @param df_all The data frame containing the TRE values and the residue number
##' over monte carlo interations. Results from the sampleTREuncertainty function.
##' @param pca_states The PCA object from the UnsupervisedCSA function
##' @param states The state name to use for the analysis. e.g. ligand used in
##'  differential analysis
##' @param whichSamples The samples to use for the analysis.
##'  Default is seq.int(1, 50)
##' 
##' @return A ggplot object
##' 
##' @examples
##' 
##' library("RexMS")
##' data(TRE_dist)
##' data("out_lxr_compound_proccessed")
##' data("LXRalpha_compounds")
##' 
##' states <- names(LXRalpha_compounds)
##' 
##' ucsa <- UnsupervisedCSA(out_lxr_compound_proccessed,
##'                        quantity = "TRE",
##'                        states = states,
##'                        whichTimepoint = 600,
##'                        pca_params = list(scale = FALSE,
##'                                          center = TRUE))
##'                        
##' 
##' plotTREuncertaintyLoadings(df_all = TRE_dist,
##'                           pca_states = ucsa$pca_states,
##'                           states = states,
##'                           whichSamples = seq.int(1, 50))
##'   
##' @export                                                     
plotTREuncertaintyLoadings <- function(df_all,
                                       pca_states,
                                       states,
                                       whichSamples = seq.int(1, 50)){
  
  ## Create coordinates
  d <- 2
  dims <- c(1,2)
  coords <- matrix(NA, nrow = length(states), ncol = 2)
  res0 <- pca_states
  #eigs <- colnames(.pca)
  R <- max(df_all$residue)
  
  ## Create inital dataset, computing average location in PCA plot
  Y0 <- matrix(NA, nrow = length(states), ncol = 2)
  for (j in seq.int(length(states))) {
    Y0[j, ] <- res0$x[j, seq_len(d)]
  }
  
  ## Repeat for different monte-carlo samples
  Y.lst <- list()
  res.rot <- list()
  for ( i in seq.int(max(df_all$mc_num))) {
    
    data <- matrix(df_all$TRE[df_all$mc_num == unique(df_all$mc_num)[i]], 
                   ncol = R, byrow = TRUE)
    res <- prcomp(data, center = TRUE)
    res.proc <- vegan::procrustes(Y0, res$x[, dims])
    Y.df <- data.frame(res.proc$Yrot, states = states)
    Y.lst[[i]] <- Y.df
    res.rot[[i]] <- res$rotation[, 1]
  }
  
  names(res.rot) <- seq_len(length(whichSamples))
  resrot.lst.df <- plyr::ldply(res.rot, .fun = function(x) x, .id = "mcmcIter")
  resrot.lst.df_long <- resrot.lst.df %>% pivot_longer(cols = 2:(R + 1))
  resrot.lst.df_long$name <- as.numeric(sub(".", "", resrot.lst.df_long$name))
  loadingquants <- resrot.lst.df_long %>% group_by(name) %>% 
    reframe(level = c(0.025, 0.5, 0.975), 
            quantile = quantile(value, c(0.025, 0.5, 0.975)))
  
  
  
  gg <- ggplot(
    data = loadingquants,
    aes(x = name,
        y = quantile,
        group = factor(level),
        col = factor(level))) +
    geom_line(lwd = 1.5) + theme_bw() + 
    scale_color_manual(values = c("black",
                                  brewer.pal(n = 3, name = "PuOr"))[c(4,1,2)]) +
    theme(text = element_text(size = 20)) + 
    labs(col = "quantile") + xlab("Residue") + ylab("Loading")
  
  return(gg)
}
