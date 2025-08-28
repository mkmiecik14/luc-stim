# Quality control of PSD Data
# Matt Kmiecik
# 05 May 2024

# Purpose: visualizes quality control of EEG data

# libraries ----
library(tidyverse)

# functions ----
source("fns/topo_interp.R"); source("fns/topo_plot.R")

# data ----
files <- 
  as.list(
    dir(
      path = "../output", 
      pattern = "2025\\-03\\-22*.\\.rda", # modify the date here as necessary
      full.names = TRUE
      )
    )
walk(files, ~load(.x, .GlobalEnv)) # loads all files
load("../output/chan-locs.rda")

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
  # fixes naming issue:
  mutate(ss = case_when(
    ss == "1461832050" ~ "1461831842050",
    ss == "1461842050" ~ "1461831842050",
    .default = ss
    )
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

this_data %>% filter(m > 50) %>% count(ss, session, stim) %>% View()
this_data %>% filter(m > 100)

this_ss <- "1461831842050"
this_data <- ss %>% filter(ss %in% this_ss) %>% mutate(m = if_else(m > 100, NA, m))
ggplot(this_data, aes(m)) + geom_histogram(binwidth = 1) + facet_wrap(~block, scales = "free")
ggplot(this_data, aes(factor(block), m)) + 
  geom_boxplot() +
  facet_wrap(~stim)

# taking a look at raw data
this_data %>% filter(is.na(m)) %>% View()
ss$ss %>% unique()

#* Problematic subjects
#* throw out 2050 for now
#* and fix the naming of 2050 in the prepro script