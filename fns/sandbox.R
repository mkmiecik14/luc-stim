# sandbox.R
# Matt Kmiecik
# Functions that I wrote but were re-purposed; keeping here in case I need them

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