# analysis-4.R
# Matt Kmiecik
# Started 06 April 2024

# Purpose: preliminary viz and analysis of LUC stim data

# source("prepro-to-r.r") # not run

# libraries ----
library(tidyverse); library(scales)

# functions ----
source("fns/topo_interp.R"); source("fns/topo_plot.R")

# data ----
files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
walk(files, ~load(.x, .GlobalEnv))

alpha <- seq(8, 12, .25)
ss <- 
  psd_res %>%
  filter(freq %in% alpha, !is.na(psd)) %>%
  mutate(
    block = block - 1, # adjusts as IAF was not collected,
    task = if_else(block > 4, "post", "pre") # inserts task
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


ggplot(ss_sum %>% filter(eyes == "open"), aes(m)) + 
  geom_density() + 
  scale_x_continuous(trans = "log2") +
  facet_grid(stim~task, scales = "free")

## summarises across electrode
ss_sum_elec <- 
  ss_sum %>% 
  group_by(ss, stim, eyes, task) %>%
  summarise(M = mean(m), SD = sd(m), N = n()) %>%
  ungroup() %>%
  mutate(exclude = if_else(M > 50, 1, 0))

# study summary
study_sum <- 
  ss_sum %>% 
  left_join(
    ., 
    ss_sum_elec %>% select(ss:task, exclude), 
    by = c("ss", "stim", "eyes", "task")
    ) %>%
  filter(exclude == 0) %>%
  group_by(stim, eyes, task, elec) %>%
  summarise(
    M = mean(m), 
    SD = sd(m), 
    N = n(), 
    SEM = SD/sqrt(N), 
    MOE = qt(.975, df = N - 1) * SEM
  ) %>%
  ungroup()

## creates subtraction
study_sum_sub <- 
  study_sum %>% 
  select(stim:M, N) %>%
  group_by(stim, eyes, elec) %>%
  mutate(M = diff(M)) %>% # subtraction
  ungroup() %>%
  filter(task == "post") %>% # removes pre as these are the same values as post
  mutate(task = "post > pre") # renames so that task is accurate

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
this_orig <-  study_sum %>% filter(eyes == "closed")
this_interp <-  interp %>% filter(eyes == "closed")
this_orig %>% 
  group_by(stim, eyes, task) %>%
  summarise(mean = mean(M), min = min(M), max = max(M))

psd_closed_plot <- 
  topo_plot(
  orig_data = this_orig, 
  interp_data = this_interp, 
  dv = M,
  color_pal_limits = c(0, 50),
  color_pal_breaks = seq(0, 50, 10),
  elec_shape_col = NULL,
  elec_shapes = 19,
  bwidth = 1.5, # width of colorbar
  bheight = .2, # height of colorbar
  d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
  contour_alpha = 1/3, # alpha level of contour lines
  contour_color = "black", # color of contour lines
  headshape_size = .5, # headshape size
  electrode_size = .5, # size of electrode points
  nose_size = .5, # size of nose shape,
  nose_adj = nose_adj, # adjusts position of nose,
  size_maskRing = 6,
  legend_name = "PSD (uV^2/Hz)"
) + 
  facet_grid(stim~task) +
  labs(title = "Eyes Closed")

psd_closed_plot

# saves out
# ggsave(
#   filename = "../output/psd-closed.png", 
#   width = 6.5, 
#   height = 5.5, 
#   units = "in",
#   bg = "white"
#   )

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
  color_pal_limits = c(-14, 14),
  color_pal_breaks = seq(-14, 14, 7),
  color_pal = rev(brewer.pal(11, "RdBu")),
  elec_shape_col = NULL,
  elec_shapes = 19,
  bwidth = 1.5, # width of colorbar
  bheight = .2, # height of colorbar
  d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
  contour_alpha = 1/3, # alpha level of contour lines
  contour_color = "black", # color of contour lines
  headshape_size = .5, # headshape size
  electrode_size = .5, # size of electrode points
  nose_size = .5, # size of nose shape,
  nose_adj = nose_adj, # adjusts position of nose,
  legend_name = "PSD Difference \n (uV^2/Hz)"
) + 
  facet_grid(stim~task)

psd_closed_sub_plot

# puts them together
library(patchwork)
psd_closed_fplot <- psd_closed_plot | psd_closed_sub_plot # final plot
psd_closed_fplot

# saves out
# ggsave(
#   filename = "../output/psd-closed-all.png", 
#   plot = psd_closed_fplot,
#   width = 7, height = 8, units = "in", bg = "white"
#   )

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
    color_pal_limits = c(0, 15),
    color_pal_breaks = seq(0, 15, 3),
    elec_shape_col = NULL,
    elec_shapes = 19,
    bwidth = 1.5, # width of colorbar
    bheight = .2, # height of colorbar
    d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
    contour_alpha = 1/3, # alpha level of contour lines
    contour_color = "black", # color of contour lines
    headshape_size = .5, # headshape size
    electrode_size = .5, # size of electrode points
    nose_size = .5, # size of nose shape,
    nose_adj = nose_adj, # adjusts position of nose,
    legend_name = "PSD (uV^2/Hz)"
    ) + 
  facet_grid(stim~task) +
  labs(title = "Eyes Open")

psd_open_plot

# ggsave(
#   filename = "../output/psd-open.png", 
#   width = 6.5, 
#   height = 5.5, 
#   units = "in",
#   bg = "white"
# )

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
    color_pal_limits = c(-10, 10),
    color_pal_breaks = seq(-10, 10, 5),
    color_pal = rev(brewer.pal(11, "RdBu")),
    elec_shape_col = NULL,
    elec_shapes = 19,
    bwidth = 1.5, # width of colorbar
    bheight = .2, # height of colorbar
    d = dia, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
    contour_alpha = 1/3, # alpha level of contour lines
    contour_color = "black", # color of contour lines
    headshape_size = .5, # headshape size
    electrode_size = .5, # size of electrode points
    nose_size = .5, # size of nose shape,
    nose_adj = nose_adj, # adjusts position of nose,
    legend_name = "PSD Difference \n (uV^2/Hz)"
  ) + 
  facet_grid(stim~task)

psd_open_sub_plot

# puts them together
psd_open_fplot <- psd_open_plot | psd_open_sub_plot # final plot
psd_open_fplot

# saves out
ggsave(
  filename = "../output/psd-open-all.png", 
  plot = psd_open_fplot,
  width = 7, height = 8, units = "in", bg = "white"
)


# testing out a quick plotting function

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
    color_pal_limits = c(this_min, this_max),
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

quick_plot(study_sum, interp, "sham", "open", 0, 40, seq(0, 40, 10)) # good
ggsave(filename = "../output/psd-open-sham.png", width = 5, height = 4, units = "in", bg = "white")
quick_plot(study_sum, interp, "sham", "closed") # crazy

quick_plot(study_sum, interp, "tacs", "open") # crazy
quick_plot(study_sum, interp, "tacs", "closed", 0, 40, seq(0, 40, 10)) # good

quick_plot(study_sum, interp, "tdcs", "open", 0, 40, seq(0, 40, 10)) # good
ggsave(filename = "../output/psd-open-tdcs.png", width = 5, height = 4, units = "in", bg = "white")
quick_plot(study_sum, interp, "tdcs", "closed") # crazy

quick_plot(study_sum, interp, "trns", "open", 0, 40, seq(0, 40, 10)) # good
ggsave(filename = "../output/psd-open-trns.png", width = 5, height = 4, units = "in", bg = "white")
quick_plot(study_sum, interp, "trns", "closed", 0, 40, seq(0, 40, 10)) # good




