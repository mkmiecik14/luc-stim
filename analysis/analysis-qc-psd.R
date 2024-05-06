# Quality control of PSD Data
# Matt Kmiecik
# 05 May 2024

# Purpose: visualizes quality control of EEG data

# libraries ----
library(tidyverse)

# functions ----
source("fns/topo_interp.R"); source("fns/topo_plot.R")

# data ----
files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
walk(files, ~load(.x, .GlobalEnv)) # loads all files

# bad data
library(readxl)
bd <-
  read_excel("../doc/ss-info.xlsx", sheet = "bad_data") %>% 
  mutate(ss = as.character(ss), session = as.character(session))

# preps data
alpha <- seq(8, 12, .25) # alpha range
ss <- 
  psd_res %>%
  filter(freq %in% alpha, !is.na(psd)) %>%
  mutate(
    block = block - 1, # adjusts as IAF was not collected,
    task = if_else(block > 4, "post", "pre") # inserts task
  ) %>%
  # turns certain cols to factors
  mutate(
    across(.cols = c(eyes, task), .fns = ~factor(.x)),
    eyes = relevel(eyes, ref = "open"), 
    task = relevel(task, ref = "pre")
  ) %>%
  group_by(ss, session, stim, block, eyes, elec, task) %>%
  summarise(m = mean(psd), n = n()) %>% # takes mean PSD across alpha band
  ungroup() %>%
  left_join(., bd, by = c("ss", "session", "stim", "block", "eyes", "task")) %>%
  filter(is.na(drop)) # drops bad data here

# plot
this_data <- 
  ss %>% filter(block == 8, eyes == "closed", task == "post")

ggplot(this_data, aes(m)) + geom_histogram(binwidth = 1)

ggplot(this_data, aes(stim, m)) + geom_boxplot() + facet_wrap(~elec)

this_data %>% filter(m > 50) %>% count(ss, session, stim)
this_data %>% filter(m > 100)
