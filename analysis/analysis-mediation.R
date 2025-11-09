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
library(tidyverse); library(lme4); library(lmerTest)

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
psd_col <- "psd_db"
thresh <- .99

# data ----

## EEG data
f <- file.path("output", "r-prepro", "psd_db.rds")
psd_data <- load_eeg_data_into_R(data_file = f, refresh = FALSE)

psd_data2 <-
  psd_data %>%
  filter(frequency %in% bands$alpha, !is.na(.data[[psd_col]])) %>%
  summarise(
    m = mean(.data[[psd_col]]), n = n(),
    .by = c(ss, session, stim, block, eyes, electrode, task)
  ) %>%
  mutate(ss = str_extract(ss, "^\\d+(?=_)")) %>% # cleans ss number
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
  cat(sprintf("Outliers outside of %.0f%% threshold are excluded...\n", thresh*100))
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

# MIGHT NOT NEED THIS BLOCK AS I AM TRYING TO JOIN ALL THE DATA
# # Calculating change in PSD from pre- to post-stim
# delta_psd <- 
#   ss_r %>% 
#   summarise(
#     M = mean(m), N = n(), .by = c(ss, session, stim, eyes, electrode, task)
#     ) %>%
#   pivot_wider(
#     id_cols = c(ss, session, stim, eyes, electrode), 
#     names_from = task, 
#     values_from = M
#     ) %>%
#   mutate(delta = post - pre)
# 
# # joins 
# delta_psd %>% 
#   left_join(
#     ., aut_data2, join_by(ss, session, stim), relationship = "many-to-many"
#     )

# looks like we have an issue with condition key between behavior and EEG
eeg_cond <- ss_r %>% select(ss, session, stim) %>% distinct() %>% rename(EEG = stim)
aut_cond <- aut_data2 %>% select(ss, session, stim) %>% distinct() %>% rename(AUT = stim)
ff_cond <- ff_data2 %>% select(ss, session, stim) %>% distinct() %>% rename(FF = stim)

cond_check <- 
  reduce(list(eeg_cond, aut_cond, ff_cond), full_join, join_by(ss, session)) %>%
  mutate(all_equal = (EEG == AUT & AUT == FF))
# writing this out for Bob to double check:
write_csv(cond_check, file = file.path("output", "r-analysis", "stim-cond-check.csv"))
  

# trying to join all the data!
mod_data <- 
  full_join(
    aut_data2, ss_r, join_by(ss, session, stim), relationship = "many-to-many"
    )

# Step 1: Mixed model for mediator (alpha change ~ stim)
mod_data_nest <- mod_data %>% nest_by(eyes, electrode) %>% filter(complete.cases(.))
mod_data %>% filter(is.na(electrode)) %>% pull(ss) -> ss_miss

mods <- 
  mod_data_nest %>% 
  mutate(
    med_model = list(lmer(m ~ 1 + session + block + stim*task + (1 | ss), data = data))
    )

mod_data %>% filter(!complete.cases(.))






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
  


# Proc ---- 
library(mediation)
library(lme4)

# Step 1: Mixed model for mediator
med_model <- lmer(M ~ X + covariates + (1|subject_id), data = df)

# Step 2: Mixed model for outcome
out_model <- lmer(Y ~ X + M + covariates + (1|subject_id), data = df)

# Step 3: Mediation analysis
med_results <- mediate(med_model, out_model,
                       treat = "X", mediator = "M",
                       boot = TRUE, sims = 1000,
                       cluster = "subject_id")  # Cluster-robust bootstrap
