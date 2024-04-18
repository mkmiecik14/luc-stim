# analysis-4.R
# Matt Kmiecik
# Started 06 April 2024

# Purpose: preliminary viz and analysis of LUC stim data

# source("prepro-to-r.r") # not run

# libraries ----
library(tidyverse)

# functions ----
source("fns/topo_interp.R"); source("fns/topo_plot.R")

# data ----
files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
walk(files, ~load(.x, .GlobalEnv))

alpha <- seq(8, 12, .25)
ss <- 
  psd_res %>%
  filter(freq %in% alpha, !is.na(psd)) %>%
  # determines stimulation
  mutate(
    stim_type = case_when(
      stim_type == 1 ~ "tdcs", 
      stim_type == 3 ~ "tacs", 
      stim_type == 4 ~ "trns"
    ),
    stim = if_else(stim_version == "B", "sham", stim_type),
    stim = factor(stim), stim = relevel(stim, ref = "sham")
  ) %>%
  mutate(
    block = block - 1, # adjusts as IAF was not collected
    ss = sub("\\_.*", "", ss),
    task = if_else(block > 4, "post", "pre")
  ) %>%
  # turns certain cols to factors
  mutate(
    across(.cols = c(eyes, task), .fns = ~factor(.x)),
    eyes = relevel(eyes, ref = "open"), 
    task = relevel(task, ref = "pre")
  )

# examines contrasts
contrasts(ss$stim)
contrasts(ss$eyes)
contrasts(ss$task)

# subject-wise summary
ss_sum <- 
  ss %>%
  group_by(ss, stim, eyes, task, elec) %>%
  summarise(m = mean(psd), n = n()) %>%
  ungroup()

library(scales)
ggplot(ss_sum %>% filter(eyes == "open"), aes(m)) + 
  geom_density() + 
  scale_x_continuous(trans = "log2") +
  facet_grid(stim~task, scales = "free")

# study summary
study_sum <- 
  ss_sum %>% 
  group_by(stim, eyes, task, elec) %>%
  summarise(
    M = mean(m), 
    SD = sd(m), 
    N = n(), 
    SEM = SD/sqrt(N), 
    MOE = qt(.975, df = N - 1) * SEM
  ) %>%
  ungroup()

ggplot(study_sum %>% filter(eyes == "closed"), aes(M)) + 
  geom_density() +
  scale_y_continuous(trans = "log10") +
  facet_grid(stim~task, scales = "free")

# interpolates data
interp <- 
  study_sum %>%
  split(interaction(.$stim, .$eyes, .$task, sep = "_")) %>%
  map_dfr(
    ~topo_interp(data = .x, meas = "M", gridRes = 100, size = .6), .id = "name"
  ) %>%
  separate(name, into = c("stim", "eyes", "task")) %>%
  as_tibble() %>%
  mutate(task = factor(task), task = fct_relevel(task, c("pre", "post")))


quick_plot <- function(orig, interp, this_stim, this_eyes, min = NULL, max = NULL, cust_seq = NULL){
  
  
  this_orig <- orig %>% filter(stim == this_stim, eyes == this_eyes)
  this_interp <-  interp %>% filter(stim == this_stim, eyes == this_eyes)
  
  ranges <- 
    this_orig %>% 
    group_by(stim, eyes, task) %>% 
    summarise(min = min(M), max = max(M)) %>%
    ungroup()
  print(ranges)
  
  if (is.null(min)) {
    this_min <- round(min(ranges$min))
  } else{ this_min <- min }
  
  if (is.null(min)) {
    this_max <- round(max(ranges$max))
  } else{ this_max <- max }
  
  if (is.null(cust_seq)) {
    this_seq <- seq(this_min, this_max, round((this_max - this_min)/3, 1))
  } else{ this_seq <- cust_seq }
  
  print(this_min); print(this_max)
  
  topo_plot(
    orig_data = this_orig, 
    interp_data = this_interp, 
    dv = M,
    color_pal_limits = c(min(ranges$min), max(ranges$max)),
    color_pal_breaks = this_seq,
    elec_shape_col = NULL,
    elec_shapes = 19,
    bwidth = 1, # width of colorbar
    bheight = .2, # height of colorbar
    d = 1.2, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
    contour_alpha = 1/3, # alpha level of contour lines
    contour_color = "black", # color of contour lines
    headshape_size = .5, # headshape size
    electrode_size = 1, # size of electrode points
    nose_size = .5, # size of nose shape,
    nose_adj = -.12, # adjusts position of nose,
    legend_name = "PSD (uV^2/Hz)"
  ) + 
    facet_grid(stim~task) +
    labs(title = paste0("Eyes ", this_eyes))
}

quick_plot(study_sum, interp, "sham", "open", 0, 10, seq(0, 10, 2)) # good
ggsave(filename = "../output/psd-open-sham.png", width = 5, height = 4, units = "in", bg = "white")
quick_plot(study_sum, interp, "sham", "closed") # crazy

quick_plot(study_sum, interp, "tacs", "open") # crazy
quick_plot(study_sum, interp, "tacs", "closed", 0, 40, seq(0, 40, 10)) # good

quick_plot(study_sum, interp, "tdcs", "open", 0, 35, seq(0, 40, 10)) # good
ggsave(filename = "../output/psd-open-tdcs.png", width = 5, height = 4, units = "in", bg = "white")
quick_plot(study_sum, interp, "tdcs", "closed") # crazy

quick_plot(study_sum, interp, "trns", "open", 0, 20, seq(0, 20, 5)) # good
ggsave(filename = "../output/psd-open-trns.png", width = 5, height = 4, units = "in", bg = "white")
quick_plot(study_sum, interp, "trns", "closed", 0, 40, seq(0, 40, 10)) # good



# plot eyes closed
this_orig <- study_sum %>% filter(stim == "tacs", eyes == "closed")
this_interp <-  interp %>% filter(stim == "tacs", eyes == "closed")
this_orig %>% group_by(stim, eyes, task) %>% summarise(min = min(M), max = max(M))
topo_plot(
  orig_data = this_orig, 
  interp_data = this_interp, 
  dv = M,
  color_pal_limits = c(0, 40),
  color_pal_breaks = seq(0, 40, 10),
  elec_shape_col = NULL,
  elec_shapes = 19,
  bwidth = 1, # width of colorbar
  bheight = .2, # height of colorbar
  d = 1.2, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
  contour_alpha = 1/3, # alpha level of contour lines
  contour_color = "black", # color of contour lines
  headshape_size = .5, # headshape size
  electrode_size = 1, # size of electrode points
  nose_size = .5, # size of nose shape,
  nose_adj = -.12, # adjusts position of nose,
  legend_name = "PSD (uV^2/Hz)"
) + 
  facet_grid(stim~task) +
  labs(title = "Eyes Closed")

this_orig <- study_sum %>% filter(stim == "tacs", eyes == "open")
this_interp <-  interp %>% filter(stim == "tacs", eyes == "open")
this_orig %>% group_by(stim, eyes, task) %>% summarise(min = min(M), max = max(M))
topo_plot(
  orig_data = this_orig, 
  interp_data = this_interp, 
  dv = M,
  color_pal_limits = c(0, 40),
  color_pal_breaks = seq(0, 40, 10),
  elec_shape_col = NULL,
  elec_shapes = 19,
  bwidth = 1, # width of colorbar
  bheight = .2, # height of colorbar
  d = 1.2, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
  contour_alpha = 1/3, # alpha level of contour lines
  contour_color = "black", # color of contour lines
  headshape_size = .5, # headshape size
  electrode_size = 1, # size of electrode points
  nose_size = .5, # size of nose shape,
  nose_adj = -.12, # adjusts position of nose,
  legend_name = "PSD (uV^2/Hz)"
) + 
  facet_grid(stim~task) +
  labs(title = "Eyes Closed")
