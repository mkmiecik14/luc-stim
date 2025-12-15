# analysis-viz-mediation-results.R
# Matt Kmiecik
# Purpose: vizualize the results of the mediation analysis

# LIBRARIES ====================================================================
library(tidyverse); library(glue); library(patchwork)

# CONSTANTS ====================================================================
CPAL <- RColorBrewer::brewer.pal(11, "RdGy")

# DATA =========================================================================
f <- file.path(
  "output", "r-analysis", "aut-psd-OS-SEMDIS-originality-med-res.rds"
  )
dd <- read_rds(f) # loads data

# isolates data and does some computation
# ( idea ): can add post-hoc correction
med_res <- dd$mediation_results %>% mutate(sig = p_value < .05)

# helper function to plot the estimates from mediation analysis
plot_med_res <- function(data, eye, eff, stim, color_pal = CPAL){
  
  tcols <- c("TRUE" = color_pal[3], "FALSE" = color_pal[8])
  tdata <- data %>% filter(eyes == eye, effect == eff, stim_condition == stim)
  # title constructor
  title_const <- glue("Eyes: {eye}, Stim: {stim}, Effect: {eff}")
  p <- 
    ggplot(tdata, aes(estimate, reorder(electrode, estimate), color = sig)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = .2) +
    geom_point() +
    scale_color_manual(values = tcols) +
    labs(
      x = "Estimate", 
      y = "Electrode", 
      caption = "95% CI error bars.", 
      color = "p<.05",
      title = title_const
      ) +
    theme_bw() +
    facet_wrap(~stim_condition, labeller = label_both)
  return(p)
  
}

# wrapper function to plot all three stims
wrap_stims <- function(stims, eye, effect){
  l <- stims %>% map(~plot_med_res(med_res, eye, effect, .x))
  return(wrap_plots(l))
}

stims <- unique(med_res$stim_condition) # vector of stims

# TOTAL EFFECT =================================================================
# The overall effect of stimulation on cognition, regardless of mechanism
wrap_stims(stims, "open", "Total Effect")
wrap_stims(stims, "closed", "Total Effect")

# ACME =========================================================================
# average causal mediation effect (indirect effect that operates thru mediator)
# The effect of stim on cognition that operates through alpha power changes
wrap_stims(stims, "open", "ACME")
wrap_stims(stims, "closed", "ACME")

# ADE ==========================================================================
# average direct effect
# the effect of stim on cognition that does NOT go through alpha power changes
wrap_stims(stims, "open", "ADE")
wrap_stims(stims, "closed", "ADE")

# PROP. MEDIATED ===============================================================
# What proportion of the total effect goes through the mediator?
wrap_stims(stims, "open", "Prop. Mediated")
wrap_stims(stims, "closed", "Prop. Mediated")


