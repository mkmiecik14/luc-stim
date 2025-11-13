# analysis-mediation.R
# Matt Kmiecik
# Purpose: perform mediation analysis

# notes:
# AUT?
# AUT SemDis from OpenScoring
# GPT Originality
# BERT SemDisMean
# BERT SemDisMAD (maximum distance)
# BERT Volume
# DSI

# FF?
# Org SemDis: sem dis of the new word vs all the previous ones
# Word to Word Sem Dis: the mean Sem dis of each word to the next one
# DSI: the sphere that englobes 
# Originality is scored in openscoring 
# SemDis: Glove or Bert (word to word MAD and Mean)



# libraries ----
library(tidyverse); library(lme4); library(mediation)

# custom functions ----
source("fns/load_eeg_data_into_R.R")

# CONFIG ----
CONFIG <- list(
  bin_width = 0.25
)

# frequency bands definitions
bands <- list(
  "delta" = seq(1, 4, CONFIG$bin_width),
  "theta" = seq(4, 8, CONFIG$bin_width),
  "alpha" = seq(8, 13, CONFIG$bin_width),
  "beta_low" = seq(13, 20, CONFIG$bin_width),
  "beta_high" = seq(20, 30, CONFIG$bin_width)
)

# function toggles
thresh <- .99
# you can also make a toggle for the freq band
# you can also make a toggle for the data set (dB or uV)

# data ----

## EEG data
f <- file.path("output", "r-prepro", "psd_db.rds")
psd_data <- load_eeg_data_into_R(data_file = f, refresh = FALSE)

# Determine which PSD column is available ----
psd_cols <- c("psd_db", "psd_uv")
available_cols <- psd_cols[psd_cols %in% names(psd_data)]

if (length(available_cols) == 0) {
  stop("No PSD column found! Expected either 'psd_db' or 'psd_uv'")
} else if (length(available_cols) > 1) {
  stop("Multiple PSD columns found: ", paste(available_cols, collapse = ", "),
       ". Please ensure only one PSD column is present.")
}

psd_col <- available_cols[1]
cat("Using PSD column:", psd_col, "\n")

psd_data2 <-
  psd_data %>%
  filter(frequency %in% bands$alpha, !is.na(.data[[psd_col]])) %>%
  summarise(
    m = mean(.data[[psd_col]]), n = n(),
    .by = c(ss, session, stim, block, eyes, electrode, task)
  ) %>%
  # turns certain cols to factors
  mutate(
    across(.cols = c(eyes, task, stim), .fns = ~factor(.x)),
    eyes = relevel(eyes, ref = "open"), 
    task = relevel(task, ref = "pre"),
    stim = relevel(stim, ref = "sham"),
    ss = factor(ss),
    session = factor(session, levels = c(1:4))
  )

# Removing outliers from EEG
if (thresh == 1) {
  cat("Outlier exclusion not performed...\n")
} else{
  cat(sprintf("Outliers outside of %.0f%% threshold will be excluded...\n", thresh*100))
}
tail <- ( 1- thresh ) / 2
lower <- quantile(psd_data2$m, tail, na.rm = TRUE)
upper <- quantile(psd_data2$m, 1-tail, na.rm = TRUE)

# Removes outliers according to specified threshold ----
# replaces outliers with NA
ss_r <- psd_data2 %>% mutate(m = if_else(between(m, lower, upper), m, NA)) 
ss_dropped <- ss_r %>% filter(is.na(m))

# ss_r is the EEG data!



## AUT data
aut_data <- read_csv(file = "data/AUTTotalTable_BobEEG_AUTMM_MultiD_word_20250218.csv")
ff_data <- read_csv(file = "data/FFTotalTable_BobEEG_FFMM_MultiD_word_20250218.csv")

# helper function to help conform data to eeg structure
warp_to_eeg <- function(data, dv_cols){
  cat("Conforming dataset to naming conventions in EEG data...")
  res <- 
    data %>%
    mutate(
      ss = gsub("S", "146183184", Subject),
      session = str_extract(Visit, "\\d"),
      rt = endTime - startTime,
      stim = tolower(gsub("DMN", "", Condition))
    ) %>%
    select(
      ss, behav_task = Task, session, order = ResponseOrder, rt, stim, all_of(dv_cols)
    ) %>%
    # setting factors / levels / contrasts
    mutate(
      session = factor(session, levels = c(1:4)),
      stim = factor(stim, levels = c("sham", "tacs", "tdcs", "trns")),
      ss = factor(ss)
    )
  cat("Done!\n")
  return(res)
}

# Dependent variables in each dataset
aut_dvs <- c(
  "OS_SEMDIS_originality", "GPT4_originality", "AUT_bert_HD_word_SemDismean_word",
  "AUT_bert_HD_word_SemDisMAD_word", "AUT_bert_HD_word_volume_word", 
  "AUT_bert_HD_word_dsi_word"
)
ff_dvs <- c(
  "FF_bert_HD_word_volume_word", "FF_bert_HD_word_SemDis_word",
  "FF_bert_HD_word_SemDismean_word", "FF_bert_HD_word_dsi_word"
)

# conforms data
aut_data2 <- warp_to_eeg(data = aut_data, dv_cols = aut_dvs)
ff_data2 <- warp_to_eeg(data = ff_data, dv_cols = ff_dvs)

# Calculating change in PSD from pre- to post-stim
delta_psd <-
  ss_r %>%
  summarise(
    M = mean(m), N = n(), .by = c(ss, session, stim, eyes, electrode, task)
    ) %>%
  pivot_wider(
    id_cols = c(ss, session, stim, eyes, electrode),
    names_from = task,
    values_from = M
    ) %>%
  mutate(delta = post - pre)

# aggregates aut data
aut_data2_sum <- 
  aut_data2 %>% 
  summarise(
    across(
      .cols = all_of(c("rt", aut_dvs)), 
      .fns = list(m = ~mean(.x, na.rm = TRUE))
      ), 
    resps = n(), # calculates number of responses
    .by = c(ss, session, stim)
    )

# links all data together
all_data <- left_join(delta_psd, aut_data2_sum, join_by(ss, session, stim)) 

mediation_analysis <- function(data, this_elec, this_eyes, this_y){
  
  cat(sprintf("Requested %s %s %s \n", this_elec, this_eyes, this_y))
  
  # constructs data
  cat("Slimming dataset ... \n")
  this_mod_data <- 
    data %>% 
    filter(electrode %in% this_elec, eyes %in% this_eyes) %>%
    dplyr::select(ss, session, stim, eyes, electrode, delta, Y = all_of(this_y))
  
  # linear mixed models
  cat("Linear mixed modeling ... \n")
  med_mod <- lmer(delta ~ 1 + session + stim + (1 | ss), data = this_mod_data)
  out_mod <- lmer(Y ~ 1 + session + stim + delta + (1 | ss), data = this_mod_data)
  
  # preallocates
  stims <- c("tacs", "trns", "tdcs") 
  med_res_list <- as.list(set_names(stims))
  
  # calculates mediation model for each stim condition vs. sham
  cat("Mediation modeling ... \n")
  med_res <- 
    med_res_list %>% 
    map(
      ~mediate(
        med_mod, out_mod, 
        treat = "stim", mediator = "delta",
        control.value = "sham",    # Explicitly specify control
        treat.value = .x,          # Explicitly specify treatment
        sims = 1000
      )
    )
  
  # gathering results into list
  res_list <- list(
    elec = this_elec,
    eyes = this_eyes,
    Y = this_y,
    mediation_results = med_res
  )
  return(res_list)
}

# determine a way to loop through
test <- mediation_analysis(data = all_data, "Fp1", "open", "rt_m")











# OS_SEMDIS_originality; GPT4_originality; AUT_bert_HD_word_SemDismean_word
# AUT_bert_HD_word_SemDisMAD_word; AUT_bert_HD_word_volume_word; AUT_bert_HD_word_dsi_word
mod <- lmer(
  AUT_bert_HD_word_dsi_word ~ 1 + Condition + Visit + ResponseOrder + (1 | Subject), 
  data = mod_data
  )
summary(mod)


# FF

mod_data <- 
  ff_data %>%
  mutate(
    Visit = factor(Visit, levels = paste0("Visit", 1:4)),
    Condition = factor(Condition, levels = c("DMNSham", "DMNtDCS", "DMNtACS", "DMNtRNS")),
    Subject = factor(Subject),
    rt = endTime - startTime
  ) %>%
  relocate(rt, .after = endTime)

# word order (number of words); rt; FF_bert_HD_word_volume_word; FF_bert_HD_word_SemDis_word;
# FF_bert_HD_word_SemDismean_word; FF_bert_HD_word_dsi_word
mod <- lmer(
  FF_bert_HD_word_dsi_word ~ 1 + Condition + Visit + ResponseOrder + (1 | Subject),
  data = mod_data
)
summary(mod)
  
