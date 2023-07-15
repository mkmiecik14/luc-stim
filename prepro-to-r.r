# Preprocessing Spectral Data
# Matt Kmiecik
# Started 15 July 2023

# Purpose: to convert the data prepared in MATLAB into R

# Packages ----
library(R.matlab)
library(tidyverse)
library(readxl)

# Loads in matlab files ----
# Loading in spectral results from matlab (.mat)
spec_files <- list.files("../output/", pattern = "*spec-res.mat")

# Loading in channel locations
chan_locs <- 
  read_xlsx("../doc/ss-info.xlsx", sheet = "elec") %>%
  mutate(labels = gsub("'", "", labels)) # removes apostrophe

# saving out channel locations
save(chan_locs, file = "../output/chan-locs.rda")
write_csv(chan_locs, file = "../output/chan-locs.csv")

# Gathering subject numbers
subjs <- 
  tibble(ss = gsub("-spec-res.mat", "", spec_files)) %>%
  mutate(
    ss_i = 1:n(), # to help join below / for specificity
    )

# trigger definitions
triggers <-
  tibble(
    trigger = c(101, 103, 105, 107, 109, 203, 205, 207, 209),
    eyes = c(
      "open", "open", "closed", "open", "closed", "open", 
      "closed", "open", "closed"
    ),
    block = 1:length(trigger),
    task = c("iaf", rep("pre", 4), rep("post", 4))
  )

# Reading in and unpacking spectral results ----
spec_res <- 
  spec_files %>%
  map(~readMat(file.path("../output/", .x))) %>%
  map("spec.res") %>%
  map(~as.matrix(.x))

# Pulling out various spectral results - - - -
# resting-state blocks along z dim in numerical order
stim_spectra  <- spec_res %>% map(1) # broadband PSD (in dB)
stim_freqs <- spec_res %>% map(2) # frequencies
freqs <- stim_freqs[[1]][,,1] # vector of frequencies
# WORKING ON THIS: FIND A WAY TO EXTRACT FREQ VECTOR FROM FIRST SS
dim(stim_freqs[[1]])[3]
1


stim_paf <- spec_res %>% map(3) # peak alpha frequency
stim_cog <- spec_res %>% map(4) # center of gravity
stim_iaf <- spec_res %>% map(5) # individual alpha freq
# broadband PSD (converted to uV^2/Hz + surface Laplacian)
stim_psd <- spec_res %>% map(6) 
# broadband PSD (converted to uV^2/Hz + surface Laplacian) corrected for 
# pink and white noise
stim_psd_cor <- spec_res %>% map(7)
stim_freqvec <- spec_res %>% map(8) # frequencies for noise correction
pwn_freqs <-  stim_freqvec[[1]][,1] # frequency vector for noise correction

# Cleanup
rm(spec_res) # removes from memory 
gc() # garbage collection