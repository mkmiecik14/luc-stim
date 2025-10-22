# bandwise_analysis.R
# Matt Kmiecik
# Purpose: analyze the PSD data using linear mixed models after QC

bandwise_analysis <- function(data_file = NULL,
    this_unit = "dB", band_name = "alpha", band = seq(8, 12, .25), 
    thresh = 0.99){
  
  # libraries ----
  library(tidyverse); library(readxl); library(patchwork)
  library(lme4); library(lmerTest); library(broom.mixed); library(rlang)
  
  # functions ----
  source("fns/topo_interp.R")
  source("fns/topo_plot.R")
  
  # CONSTANTS ----
  CONFIG <- list(
    rdgy = RColorBrewer::brewer.pal(11, "RdGy"),
    dark2 = RColorBrewer::brewer.pal(8, "Dark2"),
    interp_size = .67
  )
  
  # data ----
  # validation
  if (!is.data.frame(data_file)) {
    stop("Data must be provided!")
  }
  
  message("=== BEGINNING ", toupper(band_name), " ANALYSIS ===")
  message("")
  
  # preprocessing 
  ss <- 
    data_file %>%
    filter(frequency %in% band, !is.na(psd_db)) %>%
    summarise(
      m = mean(psd_db), n = n(), 
      .by = c(ss, session, stim, block, eyes, electrode, task)
      ) %>%
    mutate(ss = str_extract(ss, "^\\d+(?=_)")) %>% # cleans ss number
    # turns certain cols to factors
    mutate(
      across(.cols = c(eyes, task, stim), .fns = ~factor(.x)),
      eyes = relevel(eyes, ref = "open"), 
      task = relevel(task, ref = "pre"),
      stim = relevel(stim, ref = "sham")
    )
  
  # Removing outliers
  if (thresh == 1) {
    cat("Outlier exclusion not performed...\n")
  } else{
    cat(sprintf("Outliers outside of %.0f%% threshold are excluded...\n", thresh*100))
  }
  tail <- ( 1- thresh ) / 2
  lower <- quantile(ss$m, tail, na.rm = TRUE)
  upper <- quantile(ss$m, 1-tail, na.rm = TRUE)
  
  # Removes outliers according to specified threshold ----
  # replaces outliers with NA
  ss_r <- ss %>% mutate(m = if_else(between(m, lower, upper), m, NA)) 
  ss_dropped <- ss_r %>% filter(is.na(m))
    
  # linear mixed modeling ----
  cat("Performing linear-mixed modeling...\n")
  
  # modeling
  mod <- 
    ss_r %>% 
    nest_by(eyes, electrode) %>%
    mutate(
      mod1 = list(
        lmer(m ~ 1 + session + block + stim*task + (1 | ss), data = data)
        )
      )
  
  ## model quality
  
  # grabs residuals
  mods <- mod$mod1 
  names(mods) <- interaction(mod$eyes, mod$electrode)
  resids <- 
    mods %>% 
    map(~tibble(resid = residuals(.))) %>% 
    list_rbind(names_to = "name") %>%
    separate(name, into = c("eyes", "electrode"))
  
  # omnibus model estimates
  cat("Extracting omnibus model statistics...\n")
  omni <- mod %>% reframe(broom::glance(mod1))
  
  # extracting coefficients
  cat("Extracting model estimates...\n")
  est <- mod %>% reframe(broom::tidy(mod1, conf.int = TRUE)) # see if this takes too long
  
  # cleans up fixed effects
  ests_fixed <- 
    est %>% 
    filter(effect == "fixed") %>% 
    group_by(eyes, term) %>%
    mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% # FDR HAPPENS REGARDLESS
    ungroup() %>%
    mutate(
      across(
        .cols = starts_with("p"), 
        .fns = ~if_else(.x < .05, "p<.05", "ns"),
        .names = "{.col}.sig"
        )
      )
  
  # for plotting function (do not remove)
  est_plot_args <- 
    expand.grid(term = unique(ests_fixed$term), eyes = unique(ests_fixed$eyes))
  
  # visualization ----
  
  # first collapse across blocks within pre and post
  ss_task <- 
    ss_r %>%
    filter(!is.na(m)) %>% # removes missing data
    summarise(psd = mean(m), len = n(), .by = c(ss, stim, eyes, task, electrode))
  
  # study summary
  study_sum <- 
    ss_task %>%
    summarise(
      M = mean(psd), 
      SD = sd(psd), 
      N = n(), 
      SEM = SD/sqrt(N), 
      MOE = qt(.975, df = N - 1) * SEM,
      .by = c(stim, eyes, task, electrode)
    )
  
  # to get shapes according to p-values
  ests_int <-
    ests_fixed %>% 
    filter(grepl("\\:", term)) %>% # filters only interaction
    mutate(
      stim = case_when(
        grepl("tacs", term) ~ "tacs",
        grepl("tdcs", term) ~ "tdcs",
        grepl("trns", term) ~ "trns",
        .default = NA
      )
      )
  
  # creates subtraction
  study_sum_sub <- 
    study_sum %>% 
    select(stim:M, N) %>%
    mutate(M = diff(M), .by = c(stim, eyes, electrode)) %>% # subtraction
    filter(task == "post") %>% # removes pre as these are the same values as post
    mutate(task = "post > pre") %>% # renames so that task is accurate
    left_join(., ests_int, by = c("stim", "eyes", "electrode"))
  
  # interpolates data
  interp <- 
    study_sum %>%
    split(interaction(.$stim, .$eyes, .$task, sep = "_")) %>%
    map_dfr(
      ~topo_interp(
        data = .x, meas = "M", gridRes = 100, size = CONFIG$interp_size, 
        elec_loc_path = "doc/chan-locs.rda"
        ), 
      .id = "name"
    ) %>%
    separate(name, into = c("stim", "eyes", "task")) %>%
    as_tibble() %>%
    mutate(task = factor(task), task = fct_relevel(task, c("pre", "post")))
  
  # interpolates subtraction
  interp_sub <-
    study_sum_sub %>%
    split(interaction(.$stim, .$eyes, .$task, sep = "_")) %>%
    map_dfr(
      ~topo_interp(
        data = .x, meas = "M", gridRes = 100, size = CONFIG$interp_size,
        elec_loc_path = "doc/chan-locs.rda"
        ), 
      .id = "name"
    ) %>%
    separate(name, into = c("stim", "eyes", "task"), sep = "\\_") %>%
    as_tibble()
  
  # Double subtraction (i.e., the interaction) ----
  
  # preallocation
  stims <- c("tacs", "tdcs", "trns")
  dub_sub <- vector("list", length = length(stims))
  names(dub_sub) <- stims
  
  # isolates sham condition as everything will have sham subtracted from it
  sham_sub <- study_sum_sub %>% filter(stim == "sham") %>% mutate(M = M*-1)
  
  for (i in 1:length(stims)) {
    this_sub <- study_sum_sub %>% filter(stim == stims[i])
    dub_sub[[i]] <- 
      bind_rows(sham_sub, this_sub) %>% 
      group_by(eyes, electrode) %>% 
      summarise(M_sub = sum(M)) %>% # this is a subtraction bc sham was *-1 above
      ungroup() %>%
      mutate(task = paste0(stims[i], " > sham"), stim = stims[i]) %>%
      left_join(., ests_int, by = c("eyes", "stim", "electrode")) # adds significance
  }
  
  dub_sub_df <- dub_sub %>% list_rbind() # combines into one df
  
  # interpolates data
  dub_sub_interp <- 
    dub_sub_df %>%
    split(interaction(.$stim, .$eyes, sep = "_")) %>%
    map_dfr(
      ~topo_interp(
        data = .x, meas = "M_sub", gridRes = 100, size = CONFIG$interp_size, 
        elec_loc_path = "doc/chan-locs.rda"
        ),
      .id = "name"
    ) %>%
    separate(name, into = c("stim", "eyes"), sep = "\\_") %>%
    as_tibble() %>%
    mutate(task = paste0(stim, " > sham")) # for plotting
  
  ## Extract significant interaction data for line plots ----
  cat("Extracting significant interactions for line graphs...\n")
  
  # Extract significant electrodes for UNCORRECTED p-values
  sig_list_uncorrected <- ests_int %>% 
    filter(p.value < .05) %>% 
    split(interaction(.$stim, .$eyes))
  
  sig_data_uncorrected <- vector("list", length = length(sig_list_uncorrected))
  names(sig_data_uncorrected) <- names(sig_list_uncorrected)
  
  for (i in seq_along(sig_list_uncorrected)) {
    this_name <- names(sig_list_uncorrected)[i]  # Get the name
    this_int <- str_split(this_name, "\\.", simplify = TRUE) # splits interaction
    this_stim <- this_int[1,1] # grabs stim
    this_eyes <- this_int[1,2] # grabs eyes
    
    # filters data
    res <- 
      study_sum %>% 
      filter(
        stim %in% c("sham", this_stim),
        eyes %in% this_eyes,
        electrode %in% sig_list_uncorrected[[i]]$electrode
      )
    
    # sets empty results to null
    if (nrow(res) == 0) {
      cat(paste0("  No results for ", this_name, " (uncorrected)\n"))
      sig_data_uncorrected[[this_name]] <- NULL
    } else{
      sig_data_uncorrected[[this_name]] <- res
    }
  }
  
  # Extract significant electrodes for FDR-CORRECTED p-values
  sig_list_fdr <- ests_int %>% 
    filter(p.fdr < .05) %>% 
    split(interaction(.$stim, .$eyes))
  
  sig_data_fdr <- vector("list", length = length(sig_list_fdr))
  names(sig_data_fdr) <- names(sig_list_fdr)
  
  for (i in seq_along(sig_list_fdr)) {
    this_name <- names(sig_list_fdr)[i]  # Get the name
    this_int <- str_split(this_name, "\\.", simplify = TRUE) # splits interaction
    this_stim <- this_int[1,1] # grabs stim
    this_eyes <- this_int[1,2] # grabs eyes
    
    # filters data
    res <- 
      study_sum %>% 
      filter(
        stim %in% c("sham", this_stim),
        eyes %in% this_eyes,
        electrode %in% sig_list_fdr[[i]]$electrode
      )
    
    # sets empty results to null
    if (nrow(res) == 0) {
      cat(paste0("  No results for ", this_name, " (FDR-corrected)\n"))
      sig_data_fdr[[this_name]] <- NULL
    } else{
      sig_data_fdr[[this_name]] <- res
    }
  }
  
  # Saving ----
  cat("Exporting results...\n")
  all_objects <- list(
    CONFIG = CONFIG,
    band = band,
    band_name = band_name,
    unit = this_unit,
    exclusion_threshold = list(thresh = thresh, lower = lower, upper = upper),
    #data = ss_r, # not saving the original data as it's too memory intensive
    data_dropped = ss_dropped,
    
    # SAVE PLOT DATA INSTEAD OF PLOTS
    power_histogram_data = ss,  # data for histogram
    histogram_thresholds = list(lower = lower, upper = upper),
    #power_histogram = plot1_all,  # REMOVE
    
    #models = mod, # KEEP COMMENTED - models are memory intensive
    #qq_plots = qq_plots,  # REMOVE - can regenerate from resids
    resids = resids,  # KEEP - small and useful for QQ plots
    
    omni = omni,
    ests = ests_fixed,
    
    #est_plots = est_plots,  # REMOVE
    est_plot_args = est_plot_args,  # KEEP - to know which plots to regenerate
    
    # SAVE DATA FOR TOPOGRAPHY PLOTS INSTEAD OF PLOT OBJECTS
    #eyes_closed_topos = eyes_closed_topos,  # REMOVE
    #eyes_open_topos = eyes_open_topos,  # REMOVE
    study_sum = study_sum,
    study_sum_sub = study_sum_sub,
    interp = interp,
    interp_sub = interp_sub,
    
    #interaction_topos = dub_sub_topos,  # REMOVE
    dub_sub_df = dub_sub_df,
    dub_sub_interp = dub_sub_interp,
    
    # SAVE BOTH CORRECTED AND UNCORRECTED SIGNIFICANT DATA for line plots
    sig_data_uncorrected = sig_data_uncorrected,
    sig_data_fdr = sig_data_fdr
  )
  return(all_objects)
}

