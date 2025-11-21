# mediation_analysis_pipeline.R
# Matt Kmiecik
# function to perform mediation analysis

# libraries ----
library(tidyverse); library(lme4)

mediation_analysis_pipeline <- function(
    eeg_data = NULL, behav_data = NULL, band_name = "alpha", 
    band = seq(8, 12, .25), thresh = 0.99, med_sims = 1000, behav_dvs = NULL){
  
  # FUNCTION START!
  # you could add validation here
  
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
  
  # preps EEG data
  psd_data2 <-
    psd_data %>%
    filter(frequency %in% band, !is.na(.data[[psd_col]])) %>%
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
  
  # aggregates BEHAVIORAL data
  behav_data_sum <- 
    behav_data %>% 
    summarise(
      across(
        .cols = all_of(c(behav_dvs)), # change this to the input
        .fns = list(m = ~mean(.x, na.rm = TRUE))
      ), 
      resps = n(), # calculates number of responses
      .by = c(ss, session, stim)
    )
  
  # links all data together
  all_data <- left_join(delta_psd, behav_data_sum, join_by(ss, session, stim)) 
  
  # helper function to preform mediation analysis
  mediation_analysis <- function(data, this_elec, this_eyes, this_y, n_sims = 5){
    
    cat(sprintf("Requested %s %s %s \n", this_elec, this_eyes, this_y))
    
    # constructs data
    cat("Slimming dataset ... \n")
    this_mod_data <- 
      data %>% 
      filter(electrode %in% this_elec, eyes %in% this_eyes) %>%
      select(ss, session, stim, eyes, electrode, delta, Y = all_of(this_y))
    
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
        ~mediation::mediate(
          med_mod, out_mod, 
          treat = "stim", mediator = "delta",
          control.value = "sham",    # Explicitly specify control
          treat.value = .x,          # Explicitly specify treatment
          sims = n_sims # change to 50-100 for development
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
  
  # determine a way to loop through ----
  
  # Define outcome variables (can be one or more)
  outcome_vars <- paste0(behav_dvs, "_m") # dynamic fron function input
  eyes_conds <- c("open", "closed")
  elecs <- as.character(unique(all_data$electrode))
  
  # Create combinations including outcome variables
  analysis_grid <- expand_grid(
    outcome = outcome_vars,
    eyes = eyes_conds,
    electrode = elecs
  )
  
  # Run mediation analysis for all combinations
  res <- 
    analysis_grid %>%
    pmap(
      function(electrode, eyes, outcome) {
        mediation_analysis(
          data = all_data, 
          this_elec = electrode, 
          this_eyes = eyes, 
          this_y = outcome,
          n_sims = med_sims # function input
        )
      },
      .progress = TRUE
    ) %>%
    set_names(
      paste(analysis_grid$electrode, analysis_grid$eyes, analysis_grid$outcome, sep = ".")
    )
  
  # Simple function to extract mediation results to a tibble
  tidy_mediate <- function(med_obj) {
    s <- summary(med_obj)
    
    tibble::tibble(
      effect = c("ACME", "ADE", "Total Effect", "Prop. Mediated"),
      estimate = c(s$d.avg, s$z.avg, s$tau.coef, s$n.avg),
      ci_lower = c(s$d.avg.ci[1], s$z.avg.ci[1], s$tau.ci[1], s$n.avg.ci[1]),
      ci_upper = c(s$d.avg.ci[2], s$z.avg.ci[2], s$tau.ci[2], s$n.avg.ci[2]),
      p_value = c(s$d.avg.p, s$z.avg.p, s$tau.p, s$n.avg.p)
    )
  }
  
  # Extract results with all identifiers preserved
  results_tidy <- 
    res %>%
    imap(~ {
      # Split the name to get electrode, eyes, and outcome
      name_parts <- str_split(.y, "\\.", simplify = TRUE)
      
      # Extract mediation results
      .x$mediation_results %>%
        map(tidy_mediate) %>%
        list_rbind(names_to = "stim_condition") %>%
        mutate(
          electrode = name_parts[1],
          eyes = name_parts[2],
          outcome = name_parts[3]
        )
    }) %>%
    list_rbind() %>%
    relocate(electrode, eyes, outcome, stim_condition)
  
  # returns tidy table of mediation results with function settings
  return_obj <- list(
    behav_data_sum = behav_data_sum, 
    band_name = band_name, 
    band = band, 
    thresh = thresh, 
    med_sims = med_sims, 
    behav_dvs = behav_dvs,
    mediation_results = results_tidy
  )
  
  return(return_obj)
}