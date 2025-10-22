# bandwise_analysis.R
# Matt Kmiecik
# Purpose: analyze the PSD data using linear mixed models after QC

bandwise_analysis <- function(data_file = NULL,
    this_unit = "dB", band_name = "alpha", band = seq(8, 12, .25), p_cor = TRUE, 
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
  plot1_all <- 
    ggplot(ss, aes(m)) + 
    geom_histogram() +
    geom_vline(xintercept = c(lower, upper), linetype = 2, color = CONFIG$rdgy[3]) +
    labs(
      x = paste0("PSD (", this_unit, ")"), 
      y = "Count", 
      caption = paste0("Lines depict ", 100*thresh, "% cutoffs.")
      ) +
    ggtitle("Electrodes x Blocks x Conditions") +
    theme_bw()
  plot1_all
  
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
  
  # QQ plot
  elec_qqplot <- function(df, teyes, velec){
    elecs <- unique(df$electrode)
    p <-
      ggplot(
        resids %>% filter(eyes == teyes, electrode %in% elecs[velec]), 
        aes(sample = resid)
      ) + 
        stat_qq(size = .5) + 
        stat_qq_line(color = "red") +
        theme_bw() +
        labs(
          x = "Theoretical Quantiles", 
          y = "Sample Quantiles", 
          title = paste0("Eyes ", teyes)
          ) +
        facet_wrap(~electrode, ncol = 8)
    return(p)
  }
  
  # PLOTS QQ PLOTS
  cat("Plotting QQ plots...\n")
  qq_plots <- list(
    eyes_open = elec_qqplot(df = resids, teyes = "open", velec = 1:64),
    eyes_closed = elec_qqplot(df = resids, teyes = "closed", velec = 1:64)
  )
  
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
  
  # function to help plot fixed effects
  plot_fixed <- function(data, tterm, teyes, p_fdr = TRUE){
    
    # Select column based on p_fdr flag
    sig_col <- if (p_fdr) {
      sym("p.fdr.sig")  # FDR-corrected p-values
    } else {
      sym("p.value.sig")  # Uncorrected p-values
    }
    
    this_est <- data %>% filter(term == tterm, eyes == teyes)
    tcol <- c("ns" = CONFIG$rdgy[8], "p<.05" = CONFIG$rdgy[3])
    this_plot<- 
      ggplot(
        this_est, 
        aes(estimate, reorder(electrode, estimate), color = !!sig_col)
        ) + 
      geom_vline(xintercept = 0, linetype = 2) +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = .2) +
      geom_point() +
      labs(
        x = "Estimate", 
        y = "Electrode", 
        title = paste("Term:", tterm, "\nEyes:", teyes),
        caption = "95% CI error bars."
        ) +
      scale_color_manual(values = tcol) +
      theme_classic()
    return(this_plot)
  }
  
  # PLOTS AND STORES FIXED-EFFECT ESTIMATES
  cat("Plotting fixed-effect estimates...\n")
  est_plot_args <- 
    expand.grid(term = unique(ests_fixed$term), eyes = unique(ests_fixed$eyes))
  est_plots <- vector("list", length = nrow(est_plot_args))
  names(est_plots) <- paste(est_plot_args$term, est_plot_args$eyes, sep = ".")
  for (i in 1:nrow(est_plot_args)) {
    est_plots[[i]] <- 
    plot_fixed(
      ests_fixed, # data 
      est_plot_args[["term"]][i], # term 
      est_plot_args[["eyes"]][i], # eyes
      p_fdr = p_cor
      )
  }
  
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
  
  # helper function for plotting etc.
  plot_topos <- function(
      data_orig, data_interp, data_sub_orig, data_sub_interp, cond = "closed", 
      unit = "dB"
  ){
    
    # TOPO CONTSTANTS
    nose_adj <- .02
    dia <- 1.45
    
    # filters data by condition
    this_orig <- data_orig %>% filter(eyes == cond)
    this_interp <- data_interp %>% filter(eyes == cond)
    this_minmax <- 
      this_orig %>% 
      group_by(stim, eyes, task) %>%
      summarise(mean = mean(M), min = min(M), max = max(M))
    this_min <- round(min(this_minmax$min))
    this_max <- round(max(this_minmax$max))
    this_color_pal_breaks <- seq(this_min, this_max, length.out = 5)
    
    # main effect of condition (eyes)
    eyes_plot <- 
      topo_plot(
        orig_data = this_orig, 
        interp_data = this_interp, 
        dv = M,
        color_pal_limits = c(min(this_color_pal_breaks), max(this_color_pal_breaks)),
        color_pal_breaks = this_color_pal_breaks,
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
        legend_name = paste0("PSD (", unit, ")"),
        elec_loc_path = "doc/chan-locs.rda"
      ) + 
      facet_grid(stim~task) +
      labs(title = paste0("eyes ", cond))
    
    # subtraction plot 
    this_sub_orig <-  data_sub_orig %>% filter(eyes == cond)
    this_sub_interp <- data_sub_interp %>% filter(eyes == cond)
    this_sub_minmax <- 
      this_sub_orig %>% 
      group_by(stim, eyes, task) %>%
      summarise(mean = mean(M), min = min(M), max = max(M))
    this_sub_min <- min(this_sub_minmax$min)
    this_sub_max <- max(this_sub_minmax$max)
    # use this value as the extreme on the diverging scale
    nn <- 
      if_else(
        abs(this_sub_min) > abs(this_sub_max), 
        abs(this_sub_min), 
        abs(this_sub_max)
      )
    # helper function for finding breaks
    symmetric_vector <- function(magnitude, num_each_side = 3) {
      total_length <- 2 * num_each_side + 1  # +1 for zero
      return(seq(-magnitude, magnitude, length.out = total_length))
    }
    this_sub_color_pal_breaks <- round(symmetric_vector(nn))
    
    sub_plot <- 
      topo_plot(
        orig_data = this_sub_orig, 
        interp_data = this_sub_interp, 
        dv = M,
        color_pal_limits = c(min(this_sub_color_pal_breaks), max(this_sub_color_pal_breaks)),
        color_pal_breaks = this_sub_color_pal_breaks,
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
        legend_name = paste0("PSD Difference \n (", unit, ")"),
        elec_loc_path = "doc/chan-locs.rda"
      ) + 
      facet_grid(stim~task)
    
    # puts them together
    comb_plot <- eyes_plot | sub_plot # final plot
    
    # creating return object
    res <- 
      list(
        eyes_plot = eyes_plot,
        sub_plot = sub_plot,
        comb_plot = comb_plot
      )
    return(res)
  }
  
  cat("Plotting topographies...\n")
  eyes_closed_topos <- 
    plot_topos(
      data_orig = study_sum, 
      data_interp = interp, 
      data_sub_orig = study_sum_sub, 
      data_sub_interp = interp_sub,
      cond = "closed",
      unit = this_unit
    )
  
  
  eyes_open_topos <- 
    plot_topos(
      data_orig = study_sum, 
      data_interp = interp, 
      data_sub_orig = study_sum_sub, 
      data_sub_interp = interp_sub,
      cond = "open",
      unit = this_unit
    )
  
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
  
  # helper function for double subtraction
  dub_sub_topo_plot <- function(
      data_orig, data_interp, cond = "closed", p_fdr = TRUE, unit = "dB"){
    
    # plotting eyes open interaction
    this_orig <- data_orig %>% filter(eyes == cond)
    this_interp <- data_interp %>% filter(eyes == cond)
    this_min_max <- 
      this_orig %>% 
      group_by(eyes, stim) %>%
      summarise(mean = mean(M_sub), min = min(M_sub), max = max(M_sub))
    this_min <- round(min(this_min_max$min))
    this_max <- round(max(this_min_max$max))
    # use this value as the extreme on the diverging scale
    nn <- 
      if_else(
        abs(this_min) > abs(this_max), 
        abs(this_min), 
        abs(this_max)
      )
    # helper function for finding breaks
    symmetric_vector <- function(magnitude, num_each_side = 3) {
      total_length <- 2 * num_each_side + 1  # +1 for zero
      return(seq(-magnitude, magnitude, length.out = total_length))
    }
    this_sub_color_pal_breaks <- round(symmetric_vector(nn), 1)
    
    # Handling pvalue correction flag ----
    
    # Select column based on p_fdr flag
    shape_col <- if (p_fdr) {
      sym("p.fdr.sig")  # FDR-corrected p-values
    } else {
      sym("p.value.sig")  # Uncorrected p-values
    }
    
    # dynamic figure caption
    caption_text <- if (p_fdr) {
      "Closed circles p<.05 (FDR-corrected)"
    } else {
      "Closed circles p<.05 (uncorrected)"
    }
    p_value_shapes <- c("ns" = 1, "p<.05" = 19) # shapes
    
    # TOPO CONTSTANTS
    nose_adj <- .02
    dia <- 1.45
    
    # the plot
    p <- 
      topo_plot(
        orig_data = this_orig, 
        interp_data = this_interp, 
        dv = M_sub,
        color_pal_limits = c(min(this_sub_color_pal_breaks), max(this_sub_color_pal_breaks)),
        color_pal_breaks = this_sub_color_pal_breaks,
        color_pal = rev(brewer.pal(11, "RdBu")),
        elec_shape_col = !!shape_col,
        elec_shapes = p_value_shapes,
        bwidth = 1.5, # width of colorbar
        bheight = .2, # height of colorbar
        d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
        contour_alpha = 1/3, # alpha level of contour lines
        contour_color = "black", # color of contour lines
        headshape_size = .5, # headshape size
        electrode_size = 1, # size of electrode points
        nose_size = .5, # size of nose shape,
        nose_adj = nose_adj, # adjusts position of nose,
        legend_name = paste0("PSD (", unit, ")"),
        elec_loc_path = "doc/chan-locs.rda"
      ) + 
      facet_wrap(~task) +
      labs(
        title = "Stim * Task (pre vs. post) Interaction\nEyes Open", 
        caption = caption_text
      )
    
    return(p)
  }
  
  # Eyes open double subtraction
  conds <- c("open", "closed")
  dub_sub_topos <- as.list(set_names(conds))
  for (i in seq_along(conds)) {
    dub_sub_topos[[i]] <- 
      dub_sub_topo_plot(
        dub_sub_df, dub_sub_interp, cond = conds[i], p_fdr = p_cor, unit = this_unit
        )
  }
  
  ## line plots to visualize any interactions
  # interaction plots
  cat("Plotting significant interactions as line graphs...\n")
  
  # grabs significant elecs
  sig_list <- ests_int %>% filter(p.value < .05) %>% split(interaction(.$stim, .$eyes))
  sig_data <- vector("list", length = length(sig_list))
  names(sig_data) <- names(sig_list)
  pd <- position_dodge(width = .2)
  stim_colors <- c(
    "sham" = CONFIG$dark2[8], 
    "tacs" = CONFIG$dark2[1], 
    "tdcs" = CONFIG$dark2[2], 
    "trns" = CONFIG$dark2[3]
  )
  
  for (i in seq_along(sig_list)) {
    this_name <- names(sig_list)[i]  # Get the name
    this_int <- str_split(this_name, "\\.", simplify = TRUE) # splits interaction
    this_stim <- this_int[1,1] # grabs stim
    this_eyes <- this_int[1,2] # grabs eyes
    
    # filters data
    res <- 
      study_sum %>% 
      filter(
        stim %in% c("sham", this_stim),
        eyes %in% this_eyes,
        electrode %in% sig_list[[i]]$electrode
        )
    
    # sets empty results to null
    if (nrow(res) == 0) {
      print(paste0("no results for ", this_name))
      sig_data[[this_name]] <- NULL
    } else{
      p <- 
        ggplot(res, aes(task, M, group = stim, color = stim)) +
        geom_point(position = pd) +
        geom_errorbar(aes(ymin = M-SEM, ymax = M+SEM, width = .1), position = pd) +
        geom_path(position = pd) +
        scale_color_manual(values = stim_colors) +
        theme_bw() +
        labs(
          x = "Task", 
          y = "Mean PSD", 
          title = paste0("Sham vs. ", this_stim, " * Task (pre vs. post) Interaction\nEyes ", this_eyes),
          caption = "Error bars are SEM."
        ) +
        theme(legend.position = "bottom") +
        facet_wrap(~electrode)
      sig_data[[this_name]] <- list(res = res, plot = p)
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
    power_histogram = plot1_all,
    #models = mod, # also dropping this as likely memory intensive
    qq_plots = qq_plots,
    omni = omni,
    ests = ests_fixed,
    est_plots = est_plots,
    eyes_closed_topos = eyes_closed_topos,
    eyes_open_topos = eyes_open_topos,
    interaction_topos = dub_sub_topos,
    interaction_line_graphs = sig_data
  )
  return(all_objects)
}

