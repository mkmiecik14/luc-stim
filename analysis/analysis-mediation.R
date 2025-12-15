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
library(tidyverse); library(lme4); library(glue)
# library(mediation) # commented out because select() conflicts with dplyr::

# custom functions ----
source("fns/load_eeg_data_into_R.R")
source("fns/mediation_analysis_pipeline.R")

# getting behavioral data (should be an input)
behav_data <- read_rds(file = file.path("output", "r-prepro", "behav-data.rds"))

# for later
aut_dvs <- c(
  "OS_SEMDIS_originality", "GPT4_originality", "AUT_bert_HD_word_SemDismean_word",
  "AUT_bert_HD_word_SemDisMAD_word", "AUT_bert_HD_word_volume_word", 
  "AUT_bert_HD_word_dsi_word", "rt"
)
ff_dvs <- c(
  "FF_bert_HD_word_volume_word", "FF_bert_HD_word_SemDis_word",
  "FF_bert_HD_word_SemDismean_word", "FF_bert_HD_word_dsi_word", "rt"
)

# CONFIG ----
CONFIG <- list(
  bin_width = 0.25,
  sims = 1000,
  thresh = 0.99
)

# frequency bands definitions
bands <- list(
  "delta" = seq(1, 4, CONFIG$bin_width),
  "theta" = seq(4, 8, CONFIG$bin_width),
  "alpha" = seq(8, 13, CONFIG$bin_width),
  "beta_low" = seq(13, 20, CONFIG$bin_width),
  "beta_high" = seq(20, 30, CONFIG$bin_width)
)


# data ----

## EEG data
f <- file.path("output", "r-prepro", "psd_db.rds")
psd_data <- load_eeg_data_into_R(data_file = f, refresh = FALSE)

# analysis ----

for (i in aut_dvs[2]) {
  
  # AUT MEDIATION ANALYSIS
  aut_med_res <- 
    mediation_analysis_pipeline(
      eeg_data = psd_data, 
      behav_data = behav_data$aut, 
      band_name = "alpha", 
      band = bands$alpha, 
      thresh = CONFIG$thresh, 
      med_sims = CONFIG$sims, 
      behav_dvs = i
    )
  
  # saves out AUT mediation results
  f_i <- gsub("\\_", "\\-", i)
  f <- file.path("output", "r-analysis", glue("aut-psd-{f_i}-med-res.rds"))
  write_rds(aut_med_res, file =  f)

}


# FF MEDIATION ANALYSIS
ff_med_res <- 
  mediation_analysis_pipeline(
    eeg_data = psd_data, 
    behav_data = behav_data$aut, 
    band_name = "alpha", 
    band = bands$alpha, 
    thresh = CONFIG$thresh, 
    med_sims = CONFIG$sims, 
    behav_dvs = ff_dvs
  )

# saves out AUT mediation results
f <- file.path("output", "analysis", "ff-psd-med-res.rds")
write_rds(aut_med_res, file = f)


# # OS_SEMDIS_originality; GPT4_originality; AUT_bert_HD_word_SemDismean_word
# # AUT_bert_HD_word_SemDisMAD_word; AUT_bert_HD_word_volume_word; AUT_bert_HD_word_dsi_word
# mod <- lmer(
#   AUT_bert_HD_word_dsi_word ~ 1 + Condition + Visit + ResponseOrder + (1 | Subject), 
#   data = mod_data
#   )
# summary(mod)
# 
# 
# # FF
# 
# mod_data <- 
#   ff_data %>%
#   mutate(
#     Visit = factor(Visit, levels = paste0("Visit", 1:4)),
#     Condition = factor(Condition, levels = c("DMNSham", "DMNtDCS", "DMNtACS", "DMNtRNS")),
#     Subject = factor(Subject),
#     rt = endTime - startTime
#   ) %>%
#   relocate(rt, .after = endTime)
# 
# # word order (number of words); rt; FF_bert_HD_word_volume_word; FF_bert_HD_word_SemDis_word;
# # FF_bert_HD_word_SemDismean_word; FF_bert_HD_word_dsi_word
# mod <- lmer(
#   FF_bert_HD_word_dsi_word ~ 1 + Condition + Visit + ResponseOrder + (1 | Subject),
#   data = mod_data
# )
# summary(mod)
  
