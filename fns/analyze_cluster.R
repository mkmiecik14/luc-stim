# analyze_cluster.R
# Matt Kmiecik
# Started 01 June 2024

# Purpose: function to analyze the mean of a cluster of electrodes

analyze_cluster <- function(clust = NULL, hz = seq(8, 12, .25), psd_cut = 100){
 
  # checks if cluster of elects waS added
  stopifnot("Cluster of electrodes was not added!" = !is.null(clust))
  
  # libraries ----
  library(tidyverse); library(readxl); library(patchwork)
  library(lme4); library(lmerTest); library(broom.mixed)
  
  # data ----
  files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
  walk(files, ~load(.x, .GlobalEnv)) # loads all files
  
  # bad data
  bd <-
    read_excel("../doc/ss-info.xlsx", sheet = "bad_data") %>% 
    mutate(ss = as.character(ss), session = as.character(session))
  
  # preps data
  ss <- 
    psd_res %>%
    filter(freq %in% hz, !is.na(psd)) %>%
    mutate(
      block = block - 1, # adjusts as IAF was not collected,
      task = if_else(block > 4, "post", "pre") # inserts task
    ) %>%
    group_by(ss, session, stim, block, eyes, elec, task) %>%
    summarise(m = mean(psd), n = n()) %>% # takes mean PSD across freq band
    ungroup() %>%
    left_join(., bd, by = c("ss", "session", "stim", "block", "eyes", "task")) %>%
    filter(is.na(drop)) %>% # drops bad data here
    # turns certain cols to factors
    mutate(
      across(.cols = c(eyes, task, stim), .fns = ~factor(.x)),
      eyes = relevel(eyes, ref = "open"), 
      task = relevel(task, ref = "pre"),
      stim = relevel(stim, ref = "sham")
    )
  
  # consider here replacing values with PSD > certain threshold?
  plot1_all <- 
    ggplot(ss, aes(m)) + 
    geom_histogram(binwidth = 5) +
    labs(x = "PSD", y = "Count", caption = "Bindwidth = 5 PSD.") +
    theme_bw()
  plot1_all
  
  # psd cutoff (default == 100 PSD) ----
  # replaces outliers with NA
  ss_r <- ss %>% mutate(m = if_else(m > psd_cut, NA, m)) 
  plot1_drop <-
    ggplot(ss_r, aes(m)) + 
    geom_histogram(binwidth = 5) +
    labs(x = "PSD", y = "Count", caption = "Bindwidth = 5 PSD.") +
    theme_bw()
  plot1_drop
  
  # Computing cluster!
  ss_r_c <- 
    ss_r %>%
    filter(!is.na(m), elec %in% clust) %>% # keeps only elecs in the cluster
    group_by(ss, session, stim, block, eyes, task) %>%
    summarise(M = mean(m), N = n()) %>% # computes mean PSD of cluster
    ungroup() %>%
    rename(m = M, n = N) # renames to be consistent with other function
  
  # linear mixed modeling ----
  
  # modeling
  mod <- 
    ss_r_c %>% 
    nest_by(eyes) %>%
    mutate(mod1 = list(lmer(m ~ 1 + stim*task + (1 | ss), data = data)))
  
  ## model quality
  
  # grabs residuals
  mods <- mod$mod1 
  names(mods) <- mod$eyes
  resids <- 
    mods %>% 
    map(~tibble(resid = residuals(.))) %>% 
    list_rbind(names_to = "eyes")
  
  # omnibus model estimates
  omni <- mod %>% reframe(broom::glance(mod1))
  
  # extracting coefficients
  est <- mod %>% reframe(broom::tidy(mod1))
  
  # cleans up fixed effects
  ests_fixed <- est %>% filter(effect == "fixed")
  
  # Visualization ----
  
  # participant-wise summary
  ss_sum <- 
    ss_r_c %>% 
    summarise(
      M = mean(m), N = n(),
      .by = c(ss, session, stim, eyes, task)
      ) %>%
    rename(m = M, n = N)
  
  # study-wise summary
  study_sum <- 
    ss_sum %>% 
    summarise(
      M = mean(m), N = n(), SD = sd(m), SEM = SD/sqrt(N),
      .by = c(stim, eyes, task)
      )
  
  # plots
  cpal <- palette.colors(palette = "Okabe-Ito") #  colors
  pd <- position_dodge(width = .2) # pos d
  pj <- position_jitter(width = .2) # pos j
  
  # includes underlying subjects
  p2 <- 
    ggplot(study_sum, aes(task, M, group = stim, color = stim)) +
    geom_point(data = ss_sum, aes(y = m), alpha = 1/3, position = pj, shape = 16, size = 1) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = M-SEM, ymax = M+SEM), width = .2, position = pd) +
    scale_color_manual(values = cpal) +
    geom_path(position = pd) +
    labs(
      x = "Time", y = "Mean Power Spectral Density", caption = "SEM error bars."
      ) +
    theme_bw() +
    facet_wrap(~eyes, ncol = 2, scales = "free")
  
  # just the means
  p3 <- 
    ggplot(study_sum, aes(task, M, group = stim, color = stim)) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = M-SEM, ymax = M+SEM), width = .2, position = pd) +
    scale_color_manual(values = cpal) +
    geom_path(position = pd) +
    labs(
      x = "Time", y = "Mean Power Spectral Density", caption = "SEM error bars."
    ) +
    theme_bw() +
    facet_wrap(~eyes, ncol = 2, scales = "free")
  
  # returns list
  rr <- 
    list(
      p1_all = plot1_all,
      p2_drop = plot1_drop,
      mods = mods,
      est = est,
      omni = omni,
      p2 = p2,
      p3 = p3
    )
  return(rr)
   
}

tc <- c("CP1", "Cz", "C2", "P1", "Pz", "P2", "POz")
res <- analyze_cluster(clust = tc)
