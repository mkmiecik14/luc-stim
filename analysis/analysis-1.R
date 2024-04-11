# analysis-1.R
# Matt Kmiecik
# Started 06 April 2024

# Purpose: preliminary viz and analysis of LUC stim data

source("prepro-to-r.r")

# libraries ----
library(tidyverse)

# data ----
files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
walk(files, ~load(.x, .GlobalEnv))


# Total Peak Alpha Frequency (PAF) and Center of Gravity (COG) ----

## subject-wise data
iaf_ss <- 
  iaf_res %>% 
  filter(!task %in% c("iaf")) %>% # removes IAF data
  mutate(
    block = block - 1, # adjusts as IAF was not collected
    ss = sub("^(\\d+)_\\d+$", "\\1", ss) # removes session info
    ) 
# perhaps set contrasts here?

## subject-wise summary
iaf_ss_sum <- 
  iaf_ss %>%
  pivot_longer(cols = c(paf, cog)) %>%
  group_by(ss, stim_type, eyes, task, name) %>%
  summarise(m = mean(value), n = n()) %>%
  ungroup()

# study summary
iaf_sum <- 
  iaf_ss_sum %>% 
  group_by(stim_type, eyes, task, name) %>%
  summarise(
    M = mean(m), 
    SD = sd(m), 
    N = n(), 
    SEM = SD/sqrt(N), 
    MOE = qt(.975, df = N-1)*SEM
    ) %>%
  ungroup()

# plot distributions and mean plots here  