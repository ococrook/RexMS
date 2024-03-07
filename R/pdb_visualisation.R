##' Function to set colour ranges for HDX data visualisation onto a pdb 
##' structure
##' 
##' 
##' Colour function
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##'
##' @param dataset The dataset for which numerical values will be colour enconded
##' @param scale_limits You can force a range of numerical values to be mapped only
##' @param cmap_name Specifies the name of the color map. Current options:
##'  "ProtDeprot", "viridis", "CSA", "ProbProtDeProt"
##' @return Returns a color function
##' @md
##' 
##' @rdname pdb-visualisation
##' @export
define_color_function <- function(dataset, 
                                  scale_limits = NULL,
                                  cmap_name = NULL){
  
  if (is.null(scale_limits)) {
    # By default work out the numerical domain from the data
    scale_limits <- c(min(dataset[!is.na(dataset)]), 
                      max(dataset[!is.na(dataset)]))
  }
  
  rlog::log_info(paste(" Your scale limits are ", round(scale_limits, 2)))
  
  if (!is.null(cmap_name) && cmap_name == "ProtDeprot") {
    rlog::log_info(" Negative values will be coloured in Blue, and positive ones on Red")
    
    colormap <- colorRamp(brewer.pal(8, "Blues"))
    domain <- c(0.0, abs(scale_limits[1]))
    color_function_deprotected <- col_bin(colormap, domain, na.color="#808080")
    
    colormap <- colorRamp(brewer.pal(8, "Reds"))
    domain <- c(0.0, abs(scale_limits[2]))
    color_function_protected <- col_bin(colormap, domain, na.color="#808080")
    
    output_function <- function(x){
      if (!is.na(x) && x >= 0.0){
        return(color_function_protected(abs(x)))
      }else{
        return(color_function_deprotected(abs(x)))
      }
    }
    
    rlog::log_warn("NA values will be coloured in grey")
    
  }else if (is.null(cmap_name) || cmap_name == "viridis") {
    rlog::log_info(" Your values will be coloured using Viridis")
    
    n_values <- length(unique(sort(dataset)))
    col_pal <- c("white", viridis(n_values))
    output_function <- col_bin(col_pal, scale_limits, na.color="#808080")
    
    rlog::log_warn("NA values will be coloured in grey")
  }else if(is.null(cmap_name) || cmap_name == "CSA") {
    
    rlog::log_info(" Your values will be coloured using custom format")
    
    colormap <- colorRamp(brewer.pal(8, "Greens"))
    domain <- c(0.0, abs(scale_limits[1]))
    color_function_deprotected <- col_bin(colormap, domain, na.color="#808080")
    
    colormap <- colorRamp(brewer.pal(8, "Purples"))
    domain <- c(0.0, abs(scale_limits[2]))
    color_function_protected <- col_bin(colormap, domain, na.color="#808080")
    
    output_function <- function(x){
      if (!is.na(x) && x >= 0.0){
        return(color_function_protected(abs(x)))
      }else{
        return(color_function_deprotected(abs(x)))
      }
    }
    
    rlog::log_warn("NA values will be coloured in grey")
    
  }else if(is.null(cmap_name) || cmap_name == "ProbProtDeProt") {
    
    rlog::log_info(" Your values will be coloured using custom format")
    
    colormap <- colorRamp(brewer.pal(8, "Blues"))
    domain <- c(0.0, abs(scale_limits[1]))
    color_function_deprotected <- col_bin(colormap, domain, na.color="#808080")
    
    colormap <- colorRamp(brewer.pal(8, "Reds"))
    domain <- c(0.0, abs(scale_limits[2]))
    color_function_protected <- col_bin(colormap, domain, na.color="#808080")
    
    output_function <- function(x){
      if (!is.na(x) && x >= 0.0){
        return(color_function_protected(abs(x)))
      }else{
        return(color_function_deprotected(abs(x)))
      }
    }
    
    rlog::log_warn("NA values will be coloured in grey")
    
  }
  
  return(output_function)
}

##' A function to map estimated (de)protection values to Blue-White-Red colormap
##'  values
##' @author Broncio Aguilar-Sanjuan and Oliver Crook
##' 
##' @param dataset A list of estimated protection (positive) and
##'  deprotection (negative) values per residue
##' @param pdb_filepath The path to the PDB file for mapping of
##'  (de)protection values
##' @param scale_limits You can force a range of numerical
##'  values to be mapped only
##' @param cmap_name Specifies the name of the color map.
##'  Current options: "ProtDeprot, "viridis", "CSA", "ProbProtDeProt"
##' @return Returns a list of colours and residue numbers to be inputted into NGLVieweR, 
##' with protection values in Red and deprotection values in Blue, and NA values in Grey.
##' @md
##' 
##' @examples
##' 
##' v <- matrix(rnorm(n = 477), nrow = 1)
##' colnames(v) <- seq.int(ncol(v))
##' 
##' v2 <- v[,seq.int(344, 477), drop = FALSE]
##' colnames(v2) <- seq.int(ncol(v2)) 
##' 
##' pdb_filepath <- system.file("extdata", "test_BRD4.pdb", mustWork = TRUE,
##'  package = "ReX")
##'
##' mycolor_parameters <- hdx_to_pdb_colours(v2, pdb = pdb_filepath,
##'  cmap_name = "ProtDeprot")
##'
##' 
##' 
##' @rdname pdb-visualisation
##' @export
hdx_to_pdb_colours <- function(dataset,
                               pdb_filepath,
                               scale_limits = NULL,
                               cmap_name = NULL){
  
  if (!file.exists(pdb_filepath)){
    stop(paste(" PDB filepath does not exist:", pdb_filepath))
  }
  
  # Extract residue numbers from available residues in PDB coordinates ---
  
  pdb_content <- read.pdb(pdb_filepath)
  sequence_from_pdb <- pdbseq(pdb_content)
  sequence_residue_numbers_pdb <- strtoi(row.names(data.frame(sequence_from_pdb)))
  
  # Report available data ---
  
  n_residues_from_pdb <- nchar(paste(sequence_from_pdb, collapse=""))
  
  # logging
  rlog::log_info(paste(" Your HDX input dataset has",
                       length(colnames(dataset)), "entries"))
  rlog::log_info(paste(" And excluding NA data you only have",
                       length(dataset[!is.na(dataset)]), "entries"))
  rlog::log_info(paste(" However, your input PDB has only",
                       n_residues_from_pdb, "residues in total"))
  
  # Work out residues and (de)protection values that can be mapped onto PDB ---
  # Note: some residues are likely to be missing in the PDB
  
  sequence_residue_numbers_hdx <- as.numeric(colnames(dataset))
  residue_numbers_hdx_pdb <- dplyr::intersect(sequence_residue_numbers_hdx,
                                              sequence_residue_numbers_pdb)
  
  dataset_for_pdb_mapping <- dataset[colnames(dataset) %in% residue_numbers_hdx_pdb]
  
  aa_sequence_from_pdb <- sequence_from_pdb[sequence_residue_numbers_pdb %in% residue_numbers_hdx_pdb]
  aa_sequence_from_pdb <- as.vector(aa_sequence_from_pdb)
  
  pdb_viewer_data <- list("values"= dataset_for_pdb_mapping,
                          "residues"= residue_numbers_hdx_pdb,
                          "aa"= aa_sequence_from_pdb)
  
  # Define colormap according to scale_limits--
  
  color_function <- define_color_function(dataset, scale_limits, cmap_name)
  
  residue_selections <- pdb_viewer_data$residues
  df <- data.frame(x = unlist(lapply(pdb_viewer_data$values, color_function)),
                   y = residue_selections)
  color_parameters <- to_list(for(i in seq_len(length(residue_selections)))
    c(df$x[i], df$y[i]))
  
  return(color_parameters)
}

##' Function to view HDX data onto a pdb structure
##' 
##' 
##' @param pdb_filepath The path to the PDB file that you wish to view
##' @param representation A character representing the representation.
##'  Default is "cartoon".
##' @param quality The quality of the plot (default is "high").
##' @param color_parameters A list of colours and residue numbers to be inputted
##'  into NGLViewer. This is the output of the `hdx_to_pdb_colours` function.
##' @return Returns a structure in the viewer panel
##' @md
##' 
##' @examples
##' 
##' library(NGLVieweR)
##' v <- matrix(rnorm(n = 477), nrow = 1)
##' colnames(v) <- seq.int(ncol(v))
##' 
##' v2 <- v[,seq.int(344, 477), drop = FALSE]
##' colnames(v2) <- seq.int(ncol(v2)) 
##' 
##' pdb_filepath <- system.file("extdata", "test_BRD4.pdb", mustWork = TRUE,
##'  package = "ReX")
##'
##' mycolor_parameters <- hdx_to_pdb_colours(v2, pdb = pdb_filepath,
##'  cmap_name = "ProtDeprot")
##'
##' view_structure(pdb_filepath = pdb_filepath,
##'  color_parameters = mycolor_parameters)
##' 
##' @export
##' @rdname pdb-visualisation
##' 
view_structure <- function(pdb_filepath,
                           representation = "cartoon",
                           quality = "high",
                           color_parameters){
  
  if (!file.exists(pdb_filepath)){
    stop(paste(" PDB filepath does not exist:", pdb_filepath))
  }
  
  # Set up the viewer ---
  
  graphics <- NGLVieweR(pdb_filepath) %>%
    stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
    addRepresentation(representation) %>%   
    setQuality(quality) %>% 
    addRepresentation(representation,
                      param = list(color=color_parameters,
                                   backgroundColor="white"))
  
  return(graphics)
}