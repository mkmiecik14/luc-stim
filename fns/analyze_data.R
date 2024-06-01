# analyze_data.R
# Matt Kmiecik
# Started 01 June 2024

# Purpose: function to analyze EEG data with options

analyze_data <- function(hz = seq(8, 12, .25), psd_cut = 100){
  
  # libraries ----
  library(tidyverse); library(readxl); library(patchwork)
  library(lme4); library(lmerTest); library(broom.mixed)
  
  # functions ----
  source("fns/topo_interp.R"); source("fns/topo_plot.R")
  
  # data ----
  files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
  walk(files, ~load(.x, .GlobalEnv)) # loads all files
  
  # bad data
  bd <-
    read_excel("../doc/ss-info.xlsx", sheet = "bad_data") %>% 
    mutate(ss = as.character(ss), session = as.character(session))
  
  # preps data
  thz <- hz # freq range
  ss <- 
    psd_res %>%
    filter(freq %in% thz, !is.na(psd)) %>%
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
  
  # linear mixed modeling ----
  
  # modeling
  mod <- 
    ss_r %>% 
    nest_by(eyes, elec) %>%
    mutate(mod1 = list(lmer(m ~ 1 + stim*task + (1 | ss), data = data)))
  
  ## model quality
  
  # grabs residuals
  mods <- mod$mod1 
  names(mods) <- interaction(mod$eyes, mod$elec)
  resids <- 
    mods %>% 
    map(~tibble(resid = residuals(.))) %>% 
    list_rbind(names_to = "name") %>%
    separate(name, into = c("eyes", "elec"))
  
  # omnibus model estimates
  omni <- mod %>% reframe(broom::glance(mod1))
  
  # extracting coefficients
  est <- mod %>% reframe(broom::tidy(mod1))
  
  # cleans up fixed effects
  ests_fixed <- 
    est %>% 
    filter(effect == "fixed") %>% 
    group_by(eyes, term) %>%
    mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>%
    ungroup() %>%
    mutate(
      across(
        .cols = starts_with("p"), 
        .fns = ~if_else(.x < .05, "p<.05", "ns"),
        .names = "{.col}.sig"
      )
    )
  
  # visualization ----
  
  # first collapse across blocks within pre and post
  ss_task <- 
    ss_r %>%
    filter(!is.na(m)) %>% # removes missing data
    group_by(ss, stim, eyes, task, elec) %>%
    summarise(psd = mean(m), len = n()) %>%
    ungroup()
  
  # study summary
  study_sum <- 
    ss_task %>%
    group_by(stim, eyes, task, elec) %>%
    summarise(
      M = mean(psd), 
      SD = sd(psd), 
      N = n(), 
      SEM = SD/sqrt(N), 
      MOE = qt(.975, df = N - 1) * SEM
    ) %>%
    ungroup()
  
  # to get shapes accordingly to p-values
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
    group_by(stim, eyes, elec) %>%
    mutate(M = diff(M)) %>% # subtraction
    ungroup() %>%
    filter(task == "post") %>% # removes pre as these are the same values as post
    mutate(task = "post > pre") %>% # renames so that task is accurate
    left_join(., ests_int, by = c("stim", "eyes", "elec"))
  
  ## TOPO CONTSTANTS ----
  interp_size <- .67
  nose_adj <- .02
  dia <- 1.45
  
  # interpolates data
  interp <- 
    study_sum %>%
    split(interaction(.$stim, .$eyes, .$task, sep = "_")) %>%
    map_dfr(
      ~topo_interp(data = .x, meas = "M", gridRes = 100, size = interp_size), .id = "name"
    ) %>%
    separate(name, into = c("stim", "eyes", "task")) %>%
    as_tibble() %>%
    mutate(task = factor(task), task = fct_relevel(task, c("pre", "post")))
  
  # interpolates subtraction
  interp_sub <-
    study_sum_sub %>%
    split(interaction(.$stim, .$eyes, .$task, sep = "_")) %>%
    map_dfr(
      ~topo_interp(data = .x, meas = "M", gridRes = 100, size = interp_size), .id = "name"
    ) %>%
    separate(name, into = c("stim", "eyes", "task"), sep = "\\_") %>%
    as_tibble()
  
  # plot eyes closed
  this_orig <- study_sum %>% filter(eyes == "closed")
  this_interp <- interp %>% filter(eyes == "closed")
  this_orig %>% 
    group_by(stim, eyes, task) %>%
    summarise(mean = mean(M), min = min(M), max = max(M))
  
  psd_closed_plot <- 
    topo_plot(
      orig_data = this_orig, 
      interp_data = this_interp, 
      dv = M,
      color_pal_limits = c(0, 30),
      color_pal_breaks = seq(0, 30, 5),
      elec_shape_col = NULL,
      elec_shapes = 19,
      bwidth = 1.5, # width of colorbar
      bheight = .2, # height of colorbar
      d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
      contour_alpha = 1/3, # alpha level of contour lines
      contour_color = "black", # color of contour lines
      headshape_size = .5, # headshape size
      electrode_size = .25, # size of electrode points
      nose_size = .5, # size of nose shape,
      nose_adj = nose_adj, # adjusts position of nose,
      size_maskRing = 6,
      legend_name = "PSD (uV^2/Hz)"
    ) + 
    facet_grid(stim~task) +
    labs(title = "Eyes Closed")
  
  psd_closed_plot
  
  # subtraction plot for eyes closed
  this_orig <-  study_sum_sub %>% filter(eyes == "closed")
  this_interp <-  interp_sub %>% filter(eyes == "closed")
  this_orig %>% 
    group_by(stim, eyes, task) %>%
    summarise(mean = mean(M), min = min(M), max = max(M))
  
  psd_closed_sub_plot <- 
    topo_plot(
      orig_data = this_orig, 
      interp_data = this_interp, 
      dv = M,
      color_pal_limits = c(-5, 5),
      color_pal_breaks = seq(-5, 5, 2.5),
      color_pal = rev(brewer.pal(11, "RdBu")),
      elec_shape_col = NULL,
      elec_shapes = 19,
      bwidth = 1.5, # width of colorbar
      bheight = .2, # height of colorbar
      d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
      contour_alpha = 1/3, # alpha level of contour lines
      contour_color = "black", # color of contour lines
      headshape_size = .5, # headshape size
      electrode_size = .25, # size of electrode points
      nose_size = .5, # size of nose shape,
      nose_adj = nose_adj, # adjusts position of nose,
      legend_name = "PSD Difference \n (uV^2/Hz)"
    ) + 
    facet_grid(stim~task)
  
  psd_closed_sub_plot
  
  # puts them together
  psd_closed_fplot <- psd_closed_plot | psd_closed_sub_plot # final plot
  psd_closed_fplot
  
  ## eyes open
  this_orig <- study_sum %>% filter(eyes == "open")
  this_interp <-  interp %>% filter(eyes == "open")
  this_orig %>% 
    group_by(stim, eyes, task) %>%
    summarise(mean = mean(M), min = min(M), max = max(M))
  
  psd_open_plot <- 
    topo_plot(
      orig_data = this_orig, 
      interp_data = this_interp, 
      dv = M,
      color_pal_limits = c(0, 10),
      color_pal_breaks = seq(0, 10, 2),
      elec_shape_col = NULL,
      elec_shapes = 19,
      bwidth = 1.5, # width of colorbar
      bheight = .2, # height of colorbar
      d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
      contour_alpha = 1/3, # alpha level of contour lines
      contour_color = "black", # color of contour lines
      headshape_size = .5, # headshape size
      electrode_size = .25, # size of electrode points
      nose_size = .5, # size of nose shape,
      nose_adj = nose_adj, # adjusts position of nose,
      legend_name = "PSD (uV^2/Hz)"
    ) + 
    facet_grid(stim~task) +
    labs(title = "Eyes Open")
  
  psd_open_plot
  
  # subtraction plot for eyes open
  this_orig <-  study_sum_sub %>% filter(eyes == "open")
  this_interp <-  interp_sub %>% filter(eyes == "open")
  this_orig %>% 
    group_by(stim, eyes, task) %>%
    summarise(mean = mean(M), min = min(M), max = max(M))
  
  psd_open_sub_plot <- 
    topo_plot(
      orig_data = this_orig, 
      interp_data = this_interp, 
      dv = M,
      color_pal_limits = c(-4, 4),
      color_pal_breaks = seq(-4, 4, 2),
      color_pal = rev(brewer.pal(11, "RdBu")),
      elec_shape_col = NULL,
      elec_shapes = 19,
      bwidth = 1.5, # width of colorbar
      bheight = .2, # height of colorbar
      d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
      contour_alpha = 1/3, # alpha level of contour lines
      contour_color = "black", # color of contour lines
      headshape_size = .5, # headshape size
      electrode_size = .25, # size of electrode points
      nose_size = .5, # size of nose shape,
      nose_adj = nose_adj, # adjusts position of nose,
      legend_name = "PSD Difference \n (uV^2/Hz)"
    ) + 
    facet_grid(stim~task)
  
  psd_open_sub_plot
  
  # puts them together
  psd_open_fplot <- psd_open_plot | psd_open_sub_plot # final plot
  psd_open_fplot
  
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
      group_by(eyes, elec) %>% 
      summarise(M_sub = sum(M)) %>% # this is a subtraction bc sham was *-1 above
      ungroup() %>%
      mutate(task = paste0(stims[i], " > sham"), stim = stims[i]) %>%
      left_join(., ests_int, by = c("eyes", "stim", "elec")) # adds significance
  }
  
  dub_sub_df <- dub_sub %>% list_rbind() # combines into one df
  
  # interpolates data
  dub_sub_interp <- 
    dub_sub_df %>%
    split(interaction(.$stim, .$eyes, sep = "_")) %>%
    map_dfr(
      ~topo_interp(data = .x, meas = "M_sub", gridRes = 100, size = interp_size), 
      .id = "name"
    ) %>%
    separate(name, into = c("stim", "eyes"), sep = "\\_") %>%
    as_tibble() %>%
    mutate(task = paste0(stim, " > sham")) # for plotting
  
  # plotting eyes open interaction
  this_orig <- dub_sub_df %>% filter(eyes == "open")
  this_interp <- dub_sub_interp %>% filter(eyes == "open")
  this_orig %>% 
    group_by(eyes, stim) %>%
    summarise(mean = mean(M_sub), min = min(M_sub), max = max(M_sub))
  
  psd_open_int_plot <- 
    topo_plot(
      orig_data = this_orig, 
      interp_data = this_interp, 
      dv = M_sub,
      color_pal_limits = c(-3, 3),
      color_pal_breaks = seq(-3, 3, 1),
      color_pal = rev(brewer.pal(11, "RdBu")),
      elec_shape_col = p.value.sig,
      elec_shapes = c(1, 19),
      bwidth = 1.5, # width of colorbar
      bheight = .2, # height of colorbar
      d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
      contour_alpha = 1/3, # alpha level of contour lines
      contour_color = "black", # color of contour lines
      headshape_size = .5, # headshape size
      electrode_size = 1, # size of electrode points
      nose_size = .5, # size of nose shape,
      nose_adj = nose_adj, # adjusts position of nose,
      legend_name = "PSD Difference \n (uV^2/Hz)"
    ) + 
    facet_wrap(~task) +
    labs(
      title = "Stim * Task (pre vs. post) Interaction\nEyes Open", 
      caption = "Closed cirlces p<.05 (uncorrected)"
    )
  
  psd_open_int_plot
  
  # plotting eyes closed interaction
  this_orig <- dub_sub_df %>% filter(eyes == "closed")
  this_interp <- dub_sub_interp %>% filter(eyes == "closed")
  this_orig %>% 
    group_by(eyes, stim) %>%
    summarise(mean = mean(M_sub), min = min(M_sub), max = max(M_sub))
  
  psd_closed_int_plot <- 
    topo_plot(
      orig_data = this_orig, 
      interp_data = this_interp, 
      dv = M_sub,
      color_pal_limits = c(-5, 5),
      color_pal_breaks = seq(-5, 5, 2.5),
      color_pal = rev(brewer.pal(11, "RdBu")),
      elec_shape_col = p.value.sig,
      elec_shapes = c(1, 19),
      bwidth = 1.5, # width of colorbar
      bheight = .2, # height of colorbar
      d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
      contour_alpha = 1/3, # alpha level of contour lines
      contour_color = "black", # color of contour lines
      headshape_size = .5, # headshape size
      electrode_size = 1, # size of electrode points
      nose_size = .5, # size of nose shape,
      nose_adj = nose_adj, # adjusts position of nose,
      legend_name = "PSD Difference \n (uV^2/Hz)"
    ) + 
    facet_wrap(~task) +
    labs(
      title = "Stim * Task (pre vs. post) Interaction\nEyes Closed", 
      caption = "Closed cirlces p<.05 (uncorrected)"
    )
  
  psd_closed_int_plot
  
  ## line plots to vizualize the tdcs interactions
  # interaction plots
  
  # grabs significant elecs
  sig_elecs <-  
    ests_int %>% 
    filter(eyes == "closed", term == "stimtdcs:taskpost", p.value < .05)
  
  # narrows data to plot
  this_data <-
    study_sum %>% 
    filter(
      stim %in% c("sham", "tdcs"), 
      elec %in% sig_elecs$elec
    )
  
  # plots tdcs * sham interaction 
  rg <- RColorBrewer::brewer.pal(11, "RdGy")
  pd <- position_dodge(width = .2)
  tdcs_sham_closed_int_plot <-
    ggplot(
      this_data %>% filter(eyes == "closed"), 
      aes(task, M, group = stim, color = stim)
    ) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = M-SEM, ymax = M+SEM, width = .2), position = pd) +
    geom_path(position = pd) +
    scale_color_manual(values = c(rg[3], rg[10])) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 30)) +
    labs(
      x = "Task", 
      y = "Mean PSD", 
      title = "tdcs (vs. sham) * Task (pre vs. post) Interaction\nEyes Closed"
    ) +
    theme(legend.position = "inside", legend.position.inside = c(.85, .25)) +
    facet_wrap(~elec)
  tdcs_sham_closed_int_plot
  
  # same plot as above but eyes open
  tdcs_sham_open_int_plot <-
    ggplot(
      this_data %>% filter(eyes == "open"), 
      aes(task, M, group = stim, color = stim)
    ) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = M-SEM, ymax = M+SEM, width = .2), position = pd) +
    geom_path(position = pd) +
    scale_color_manual(values = c(rg[3], rg[10])) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 10)) +
    labs(
      x = "Task", 
      y = "Mean PSD", 
      title = "tdcs (vs. sham) * Task (pre vs. post) Interaction\nEyes Open"
    ) +
    theme(legend.position = "inside", legend.position.inside = c(.85, .25)) +
    facet_wrap(~elec)
  tdcs_sham_open_int_plot
  
  # returns list to user ----
  rr <- 
    list(
      p1_all = plot1_all, 
      p1_drop = plot1_drop,
      omni = omni,
      est = est,
      mods = mods,
      ec_topo = psd_closed_fplot,
      eo_topo = psd_closed_fplot,
      ec_topo_i = psd_closed_int_plot,
      eo_topo_i = psd_open_int_plot
      )
  return(rr)
  
}
