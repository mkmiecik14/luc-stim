# extract_spectral_data.R
# Matt Kmiecik
# Purpose: extract spectral data from MATLAB and prepare for R

# libraries ----
library(R.matlab)
library(tidyverse)

extract_spectral_data <- function(fpath){
  
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
  
  # PREALLOCATIONS ----
  block_names <- as.vector(unlist(d$block.names)) # BLOCK NAMES
  block_desc <- as.vector(unlist(d$block.descriptions)) # BLOCK DESCRIPTIONS
  elec <- as.vector(unlist(d$meta.channels.after.interp)) # ELECTRODES
  freqs <- as.vector(d$freqs[,1]) # FREQUENCIES
  
  # METADATA ----
  meta.original.channels <- as.vector(d$meta.original.channels)
  meta.final.channels <- as.vector(d$meta.final.channels)
  meta.channels.before.interp <- as.vector(unlist(d$meta.channels.before.interp))
  meta.channels.after.interp <- as.vector(unlist(d$meta.channels.after.interp))
  meta.rejected.ics <- as.vector(d$meta.rejected.ics)
  meta.rejected.ic.labels <- as.vector(unlist(d$meta.rejected.ic.labels))
  meta.n.blocks.requested <- as.vector(d$meta.n.blocks.requested)
  meta.n.blocks.processed <- as.vector(d$meta.n.blocks.processed)
  meta.n.frequency.bins <- as.vector(d$meta.n.frequency.bins)
  meta.frequency.resolution <- as.vector(d$meta.frequency.resolution)
  meta.mean.psd.power <- as.vector(d$meta.mean.psd.power)
  
  # BLOCK-WISE PSD IN dB ----
  
  # Add dimension names first
  dimnames(d$spectra.db) <- 
    list(electrode = elec, frequency = freqs, block = block_names)
  
  # Melt to long format
  psd_db <- 
    reshape2::melt(d$spectra.db, value.name = "psd_db") %>% 
    as_tibble() %>%
    rename(trigger = block) %>%
    left_join(., triggers, by = "trigger")
  
  # BLOCK-WISE PSD in uV^2/Hz ----
  
  # Add dimension names first
  dimnames(d$spectra.psd) <- 
    list(electrode = elec, frequency = freqs, block = block_names)
  
  # Melt to long format
  psd_uv <- 
    reshape2::melt(d$spectra.psd, value.name = "psd_uv") %>% 
    as_tibble() %>%
    rename(trigger = block) %>%
    left_join(., triggers, by = "trigger")
  
  # PAF ----
  rownames(d$paf) <- elec 
  colnames(d$paf) <- block_names 
  paf <- 
    as_tibble(d$paf, rownames = "elec") %>%
    pivot_longer(cols = -elec, names_to = "trigger", values_to = "paf") %>%
    mutate(trigger = as.numeric(trigger)) %>%
    left_join(., triggers, by = "trigger")
    
  # COG ----
  rownames(d$cog) <- elec 
  colnames(d$cog) <- block_names
  cog <- 
    as_tibble(d$cog, rownames = "elec") %>%
    pivot_longer(cols = -elec, names_to = "trigger", values_to = "cog") %>%
    mutate(trigger = as.numeric(trigger)) %>%
    left_join(., triggers, by = "trigger")
  
  # IAF ----
  colnames(d$iaf) <- c("paf", "cog")
  rownames(d$iaf) <- block_names
  iaf <- 
    as_tibble(d$iaf, rownames = "trigger") %>%
    mutate(trigger = as.numeric(trigger)) %>%
    left_join(., triggers, by = "trigger")
  
  # COMBINES INTO RESULTS ----
  res <- list(
    subject = ssid,
    date = proc_time,
    psd_db = psd_db,
    psd_uv = psd_uv,
    paf = paf,
    cog = cog,
    iaf = iaf,
    param.wsize = param.wsize,
    param.overlap = param.overlap,
    param.frange = param.frange,
    param.alpha.window = param.alpha.window,
    param.cmin = param.cmin,
    param.fw = param.fw, 
    param.poly = param.poly,
    meta.original.channels = meta.original.channels, 
    meta.final.channels  = meta.final.channels,
    meta.channels.before.interp = meta.channels.before.interp, 
    meta.channels.after.interp = meta.channels.after.interp,
    meta.rejected.ics = meta.rejected.ics,
    meta.rejected.ic.labels = meta.rejected.ic.labels,
    meta.n.blocks.requested = meta.n.blocks.requested,
    meta.n.blocks.processed  = meta.n.blocks.processed,
    meta.n.frequency.bins = meta.n.frequency.bins,
    meta.frequency.resolution = meta.frequency.resolution,
    meta.mean.psd.power = meta.mean.psd.power
  )
  
  cat("Finished!\n")
  
  return(res)
  
}
