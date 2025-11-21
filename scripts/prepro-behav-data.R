# prepro-behav-data.R
# Matt Kmiecik
# Purpose: preprocess the behavioral data

# libraries ----
library(tidyverse)

# data ----

## AUT data
f <- file.path("data", "AUTTotalTable_BobEEG_AUTMM_MultiD_word_20250218.csv")
aut_data <- read_csv(f)

## FF data
f <- file.path("data", "FFTotalTable_BobEEG_FFMM_MultiD_word_20250218.csv")
ff_data <- read_csv(f)

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

# Dependent variables in each dataset (these can be modified to add more)
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

# combines into list
res <- list(aut = aut_data2, ff = ff_data2)

# saves out
write_rds(res, file = file.path("output", "r-prepro", "behav-data.rds"))