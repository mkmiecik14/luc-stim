# bandwise-analyses.R
# Matt Kmiecik
# Purpose: conduct bandwise analyses on the LUC stim data

# libraries----
library(tidyverse)

# functions ----
source("fns/bandwise_analysis.R")

# CONFIG ----
CONFIG <- list(
  bin_width = 0.25
)

# Bands ----
bands <- list(
  "delta" = seq(1, 4, CONFIG$bin_width),
  "theta" = seq(4, 8, CONFIG$bin_width),
  "alpha" = seq(8, 13, CONFIG$bin_width),
  "beta_low" = seq(13, 20, CONFIG$bin_width),
  "beta_high" = seq(20, 30, CONFIG$bin_width)
)

# data ----

# helper function to load data
load_data <- function(data_file = NULL, refresh = FALSE){
  if (refresh) {
    cat("Refreshing data...\n")
    source("scripts/prepro-to-r.r")
  } else{
    cat("Using previously saved data.\n")
  }
  res <- read_rds(file = data_file)
  return(res)
}

# loads data
cat("Loading data...\n")
f <- file.path("output", "r-prepro", "psd_db.rds")
dd <- load_data(data_file = f, refresh = FALSE)

# proc ----
res <- 
  bands %>%
  imap(
    ~bandwise_analysis(
      data_file = dd, this_unit = "dB", band_name = .y, band = .x, p_cor = FALSE,
      thresh = 0.99
      )
    )

# saves out ----
f_out <- file.path("output", "r-analysis", "bandwise-analyses-db.rds")
write_rds(res, f_out)
