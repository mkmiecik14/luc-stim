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
    ss_i = seq_len(n()), # to help join below / for specificity
    )

# trigger definitions
triggers <-
  tibble(
    trigger = c(101, 103, 105, 107, 109, 203, 205, 207, 209),
    eyes = c(
      "open", "open", "closed", "open", "closed", "open", 
      "closed", "open", "closed"
    ),
    block = seq_along(trigger),
    task = c("iaf", rep("pre", 4), rep("post", 4))
  )

# stimulation table
stim_table <- 
  read_excel("../doc/ss-info.xlsx", sheet = "main") %>%
  select(ss, session, stim_type, stim_version)

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

# this loop finds the first set of non-NULL frequency vector
for (i in 1:length(stim_freqs)){
  for (j in 1:dim(stim_freqs[[i]])[3]){
    this_len <- length(stim_freqs[[i]][, , j])
    if(sum(is.na(stim_freqs[[i]][, , j])) == this_len) {
      # skips here
    } else {
      freqs <- stim_freqs[[i]][, , j] # saves the freqs here
      break # breaks the loop once a non-NULL set is found
    }
  }
}

# IAF measures
stim_paf <- spec_res %>% map(3) # peak alpha frequency
stim_cog <- spec_res %>% map(4) # center of gravity
stim_iaf <- spec_res %>% map(5) # individual alpha freq

# broadband PSD (converted to uV^2/Hz)
stim_psd <- spec_res %>% map(6)

# broadband PSD (converted to uV^2/Hz) corrected for
# pink and white noise
stim_psd_cor <- spec_res %>% map(7)
stim_freqvec <- spec_res %>% map(8) # frequencies for noise correction
pwn_freqs <-  stim_freqvec[[1]][, 1] # frequency vector for noise correction

# Cleanup
rm(spec_res) # removes from memory
gc() # garbage collection

# PEAK ALPHA FREQUENCY ----

# Cleaning up peak alpha frequency
paf_res <-
  stim_paf %>%
  map_df(~as_tibble(.x) %>% mutate(elec = chan_locs$labels), .id = "ss_i") %>%
  mutate(ss_i = as.numeric(ss_i)) %>%
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>%
  select(-ss_i) %>%
  pivot_longer(cols = c(-ss, -elec), names_to = "block", values_to = "paf") %>%
  mutate(block = as.numeric(regmatches(block, regexpr("\\d", block)))) %>%
  left_join(., triggers %>% select(block, eyes, task), by = "block") %>%
  left_join(., stim_table, by = "ss") %>%
  select(ss, session, stim_type, stim_version, elec, block, eyes, task, paf)

# Saving out paf results
save(paf_res, file = "../output/paf-res.rda")
write_csv(paf_res, file = "../output/paf-res.csv")
rm(paf_res) # removes from memory
gc() # garbage collection

# CENTER OF GRAVITY ----

# Cleaning up cog results
cog_res <-
  stim_cog %>%
  map_df(~as_tibble(.x) %>% mutate(elec = chan_locs$labels), .id = "ss_i") %>%
  mutate(ss_i = as.numeric(ss_i)) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>%
  select(-ss_i) %>%
  pivot_longer(cols = c(-ss, -elec), names_to = "block", values_to = "cog") %>%
  mutate(block = as.numeric(regmatches(block, regexpr("\\d", block)))) %>%
  left_join(., triggers %>% select(block, eyes, task), by = "block") %>%
  left_join(., stim_table, by = "ss") %>%
  select(ss, session, stim_type, stim_version, elec, block, eyes, task, cog)

# Saving out cog results
save(cog_res, file = "../output/cog-res.rda")
write_csv(cog_res, file = "../output/cog-res.csv")
rm(cog_res) # removes from memory
gc() # garbage collection

# INDIVIDUAL ALPHA FREQUENCY - GRANDAVERAGES ----

# Cleaning up iaf results
iaf_res <-
  stim_iaf %>%
  map_df(
    ~as_tibble(.x) %>%
      rename(paf = V1, cog = V2) %>%
      mutate(block = 1:n()),
    .id = "ss_i"
    ) %>%
  mutate(ss_i = as.numeric(ss_i)) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>%
  left_join(., triggers %>% select(block, eyes, task), by = "block") %>%
  left_join(., stim_table, by = "ss") %>%
  select(ss, session, stim_type, stim_version, block, eyes, task, paf, cog)

# Saving out iaf results
save(iaf_res, file = "../output/iaf-res.rda")
write_csv(iaf_res, file = "../output/iaf-res.csv")
rm(iaf_res) # removes from memory
gc() # garbage collection

# PSD data ----
# Extracting and tidying spectral results

# frequency resolution is determined in MATLAB spectral decomposition scripts
# so see those scripts for these numbers:
eeg_srate <- 256 # data sampled at 256 Hz
fft_window <- 4 # in seconds
n_points <- eeg_srate*fft_window # number of data points in FFT window
freq_bins <- eeg_srate/n_points # frequency bins

# defining frequency bands (for later)
delta <-  seq(.5, 4, freq_bins)
theta <-  seq(4, 7.5, freq_bins)
alpha <-  seq(7.5, 13, freq_bins)
beta <-   seq(13, 30, freq_bins)

# Each participant has:
# 64 rows (electrodes) * n (depends on frequency bins) * block (rs blocks)
psd_res <-
  stim_psd %>%
  map_dfr(
    ~bind_rows(apply(.x, 3, as_tibble), .id = "block") %>% 
      mutate(elec = rep(chan_locs$labels, nrow(triggers))) %>%
      relocate(elec, .after = block) %>%
      pivot_longer(c(-block, -elec)) %>%
      mutate(name = rep(freqs, nrow(triggers)*length(chan_locs$labels))),
    .id = "ss_i"
  ) %>%
  # converts to numeric
  mutate(
    ss_i = as.numeric(ss_i),
    block = as.numeric(block)
  ) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>% # joins with subject numbers
  left_join(., triggers %>% select(block, eyes, task), by = "block") %>%
  left_join(., stim_table, by = "ss") %>%
  rename(freq = name, psd = value) %>% # renaming
  # reordering + deselecting (ss_i)
  select(ss, session, stim_type, stim_version, block, eyes, elec, freq, psd) %>% 
  # injects frequency band labels here
  mutate(
    band = case_when(
      freq %in% delta ~ "delta",
      freq %in% theta ~ "theta",
      freq %in% alpha ~ "alpha",
      freq %in% beta ~ "beta",
      TRUE ~ "outside"
    )
  )

# Saving out
save(psd_res, file = "../output/psd-res.rda")
write_csv(psd_res, file = "../output/iaf-res.csv")
rm(psd_res) # removes from memory
gc() # garbage collection

# Cleans workspace objects ----
rm(
  chan_locs, freqs, i, j, pwn_freqs, spec_files, stim_cog, stim_freqs, 
  stim_freqvec, stim_iaf, stim_paf, stim_psd, stim_psd_cor, stim_spectra, subjs,
  this_len, triggers, stim_table, alpha, beta, delta, eeg_srate, fft_window, 
  freq_bins, n_points, theta
)
gc() # garbage collection
