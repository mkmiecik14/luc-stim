# analysis-3.R
# Matt Kmiecik
# Started 15 April 2024

# Purpose: preliminary viz and analysis of LUC stim data

# source("prepro-to-r.r") # not run

# libraries ----
library(tidyverse); library(purrr)

# functions ----
source("fns/topo_interp.R"); source("fns/topo_plot.R")

# data ----
files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
walk(files, ~load(.x, .GlobalEnv))

# stimulation types:
# 1 == tdcs
# 3 == tacs
# 4 == trns

# Starting with COG results ----
ss <- 
  paf_res %>%
  filter(!task %in% c("iaf")) %>% # removes IAF data
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
    #ss = sub("^(\\d+)_\\d+$", "\\1", ss) # removes session info
    ss = sub("\\_.*", "", ss)
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
  filter(!is.na(paf)) %>% # removes missing data 
  group_by(ss, stim, eyes, task, elec) %>%
  summarise(m = mean(paf), n = n()) %>%
  ungroup()

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


# interpolates data
interp <- 
  study_sum %>%
  split(interaction(.$stim, .$eyes, .$task, sep = "_")) %>%
  map_dfr(
    ~topo_interp(data = .x, meas = "M", gridRes = 100, size = .6), .id = "name"
  ) %>%
  separate(name, into = c("stim", "eyes", "task"))

# plot eyes closed
this_orig <- study_sum %>% filter(eyes == "closed")
this_interp <-  interp %>% filter(eyes == "closed")
topo_plot(
  orig_data = this_orig, 
  interp_data = this_interp, 
  dv = M,
  color_pal_limits = c(8, 12),
  color_pal_breaks = seq(8, 12, 1),
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
  legend_name = "PAF (Hz)"
) + 
  facet_grid(task~stim) +
  labs(title = "Eyes Closed")

# plot eyes open
this_orig <- study_sum %>% filter(eyes == "open")
this_interp <-  interp %>% filter(eyes == "open")
topo_plot(
  orig_data = this_orig, 
  interp_data = this_interp, 
  dv = M,
  color_pal_limits = c(8, 12),
  color_pal_breaks = seq(8, 12, 1),
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
  legend_name = "PAF (Hz)"
) + 
  facet_grid(task~stim) +
  labs(title = "Eyes Open")

# difference plots ----
post_g_pre <- 
  study_sum %>% 
  group_by(stim, eyes, elec) %>% 
  mutate(diff = diff(M)) %>% 
  ungroup() %>% 
  filter(task == "post")

# interpolates data
post_g_pre_interp <- 
  post_g_pre %>%
  split(interaction(.$stim, .$eyes, sep = "_")) %>%
  map_dfr(
    ~topo_interp(data = .x, meas = "diff", gridRes = 100, size = .6), .id = "name"
  ) %>%
  separate(name, into = c("stim", "eyes"))

# plot eyes closed
this_orig <- post_g_pre %>% filter(eyes == "closed")
this_interp <-  post_g_pre_interp %>% filter(eyes == "closed")
plot1 <- 
  topo_plot(
    orig_data = this_orig, 
    interp_data = this_interp, 
    dv = diff,
    color_pal_limits = c(-.6, .6),
    color_pal_breaks = seq(-.6, .6, .3),
    color_pal = rev(brewer.pal(11, "RdBu")),
    elec_shape_col = NULL,
    elec_shapes = 19,
    bwidth = 1.25, # width of colorbar
    bheight = .2, # height of colorbar
    d = 1.2, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
    contour_alpha = 1/3, # alpha level of contour lines
    contour_color = "black", # color of contour lines
    headshape_size = .5, # headshape size
    electrode_size = .25, # size of electrode points
    nose_size = .5, # size of nose shape,
    nose_adj = -.12, # adjusts position of nose,
    legend_name = "PAF (Hz)"
  ) + 
  facet_wrap(~stim, nrow = 1) +
  labs(title = "Eyes Closed: Post > Pre")
plot1

# saves out
ggsave(
  filename = "../output/paf-closed-diff-topo.png", 
  plot = plot1, 
  width = 5, 
  height = 4.5, 
  units = "in",
  bg = "white"
)




