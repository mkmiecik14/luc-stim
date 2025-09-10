# extract_spectral_data.R
# Matt Kmiecik
# Purpose: extract spectral data from MATLAB and prepare for R

# libraries ----
library(R.matlab)
library(tidyverse)

extract_spectral_data <- function(fpath){
  
  # Input validation ----
  if (!file.exists(fpath)) {
    stop("File not found: ", fpath)
  }
  if (!grepl("\\.mat$", fpath, ignore.case = TRUE)) {
    stop("File must have .mat extension: ", fpath)
  }
  
  # trigger definitions
  triggers <-
    tibble(
      trigger = c(103, 105, 107, 109, 203, 205, 207, 209),
      eyes = c(
        "open", "closed", "open", "closed", "open", 
        "closed", "open", "closed"
      ),
      block = seq_along(trigger),
      task = c(rep("pre", 4), rep("post", 4))
    )
  
  
  # Reading in and unpacking spectral results ----
  spec_data <- readMat(fpath) 
  d <- spec_data$spec.results[,,1] # extracts further
  
  # SUBJECT ID ----
  ssid <- as.vector(d$subject.id)
  cat(sprintf("Processing participant: %s ...", ssid))
  
  # DATETIME ----
  proc_time <- as.POSIXct(d$meta.processing.complete)
  
  # PARAMS ----
  param.wsize <- as.vector(d$param.wsize) # WINDOW SIZE
  param.overlap <- as.vector(d$param.overlap) # OVERLAP
  param.frange <- as.vector(d$param.frange) # FREQUENCY-RANGE
  param.alpha.window <- as.vector(d$param.alpha.window) # ALPHA WINDOW
  param.cmin <- as.vector(d$param.cmin) # Min channels for cross-channel avg
  param.fw <- as.vector(d$param.fw) # SGF frame width
  param.poly <- as.vector(d$param.poly) # SGF polynomial order
  
  # Helper function for array reshaping ----
  reshape_array_data <- function(array, elec, freqs, block_names, value_name) {
    dimnames(array) <- list(electrode = elec, frequency = freqs, block = block_names)
    reshape2::melt(array, value.name = value_name) %>% 
      as_tibble() %>%
      rename(trigger = block) %>%
      left_join(triggers, by = "trigger")
  }
  
  # PREALLOCATIONS ----
  block_names <- as.vector(unlist(d$block.names))
  elec <- as.vector(unlist(d$meta.channels.after.interp))
  freqs <- as.vector(d$freqs[,1])
  
  # METADATA ----
  # Extract all metadata with consistent pattern
  metadata <- list(
    original.channels = as.vector(d$meta.original.channels),
    final.channels = as.vector(d$meta.final.channels),
    channels.before.interp = as.vector(unlist(d$meta.channels.before.interp)),
    channels.after.interp = as.vector(unlist(d$meta.channels.after.interp)),
    rejected.ics = as.vector(d$meta.rejected.ics),
    rejected.ic.labels = as.vector(unlist(d$meta.rejected.ic.labels)),
    n.blocks.requested = as.vector(d$meta.n.blocks.requested),
    n.blocks.processed = as.vector(d$meta.n.blocks.processed),
    n.frequency.bins = as.vector(d$meta.n.frequency.bins),
    frequency.resolution = as.vector(d$meta.frequency.resolution),
    mean.psd.power = as.vector(d$meta.mean.psd.power)
  )
  
  # SPECTRAL DATA ----
  psd_db <- reshape_array_data(d$spectra.db, elec, freqs, block_names, "psd_db")
  psd_uv <- reshape_array_data(d$spectra.psd, elec, freqs, block_names, "psd_uv")
  
  # FREQUENCY METRICS ----
  # Helper function for electrode-wise metrics
  reshape_elec_data <- function(matrix, elec, block_names, value_name) {
    rownames(matrix) <- elec 
    colnames(matrix) <- block_names
    as_tibble(matrix, rownames = "elec") %>%
      pivot_longer(cols = -elec, names_to = "trigger", values_to = value_name) %>%
      mutate(trigger = as.numeric(trigger)) %>%
      left_join(triggers, by = "trigger")
  }
  
  paf <- reshape_elec_data(d$paf, elec, block_names, "paf")
  cog <- reshape_elec_data(d$cog, elec, block_names, "cog")
  
  # IAF (individual alpha frequency) ----
  colnames(d$iaf) <- c("paf", "cog")
  rownames(d$iaf) <- block_names
  iaf <- 
    as_tibble(d$iaf, rownames = "trigger") %>%
    mutate(trigger = as.numeric(trigger)) %>%
    left_join(triggers, by = "trigger")
  
  # COMBINES INTO RESULTS ----
  res <- list(
    # Core data
    subject = ssid,
    date = proc_time,
    psd_db = psd_db,
    psd_uv = psd_uv,
    paf = paf,
    cog = cog,
    iaf = iaf,
    
    # Processing parameters
    param.wsize = param.wsize,
    param.overlap = param.overlap,
    param.frange = param.frange,
    param.alpha.window = param.alpha.window,
    param.cmin = param.cmin,
    param.fw = param.fw, 
    param.poly = param.poly
  )
  
  # Add metadata using list concatenation
  res <- c(res, setNames(metadata, paste0("meta.", names(metadata))))
  
  cat("Finished!\n")
  
  return(res)
  
}
