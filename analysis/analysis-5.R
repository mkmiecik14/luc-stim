# analysis-5.R
# Matt Kmiecik
# Started 05 May 2024

# Purpose: analyze the PSD data using linear mixed models after QC

# libraries ----
library(tidyverse)

# functions ----
source("fns/topo_interp.R"); source("fns/topo_plot.R")

# data ----
files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
walk(files, ~load(.x, .GlobalEnv)) # loads all files

# bad data
library(readxl)
bd <-
  read_excel("../doc/ss-info.xlsx", sheet = "bad_data") %>% 
  mutate(ss = as.character(ss), session = as.character(session))

# preps data
alpha <- seq(8, 12, .25) # alpha range
ss <- 
  psd_res %>%
  filter(freq %in% alpha, !is.na(psd)) %>%
  mutate(
    block = block - 1, # adjusts as IAF was not collected,
    task = if_else(block > 4, "post", "pre") # inserts task
  ) %>%
  group_by(ss, session, stim, block, eyes, elec, task) %>%
  summarise(m = mean(psd), n = n()) %>% # takes mean PSD across alpha band
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
ggplot(ss, aes(m)) + geom_histogram(binwidth = 5)
ss %>% filter(m > 100) %>% arrange(-m)
ss %>% summarise(M = mean(m), sd = sd(m))

# let's try 100 PSD
ss_r <- ss %>% mutate(m = if_else(m > 100, NA, m)) # replaces outliers with NA
ggplot(ss_r, aes(m)) + geom_histogram(binwidth = 5)

# linear mixed modeling ----
library(lme4); library(lmerTest); library(broom.mixed) # pkgs

# modeling
mod <- 
  ss_r %>% 
  nest_by(eyes, elec) %>% # modulate block by 1 here in future
  mutate(mod1 = list(lmer(m ~ 1 + block + stim*task + (1 | ss), data = data)))

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

# function to help plot fixed effects
plot_fixed <- function(data, tterm, teyes){
  this_est <- data %>% filter(term == tterm, eyes == teyes)
  this_plot<- 
    ggplot(this_est, aes(estimate, reorder(elec, estimate), color = p.value.sig)) + 
    geom_point() +
    labs(x = "Estimate", y = "Electrode", title = paste("Term:", tterm, "\nEyes:", teyes)) +
    theme_bw()
  return(this_plot)
}

# interaction effects
plot_fixed(ests_fixed, "(Intercept)", "open")
plot_fixed(ests_fixed, "(Intercept)", "closed")

plot_fixed(ests_fixed, "stimtacs:taskpost", "open")
plot_fixed(ests_fixed, "stimtacs:taskpost", "closed")

plot_fixed(ests_fixed, "stimtdcs:taskpost", "open")
plot_fixed(ests_fixed, "stimtdcs:taskpost", "closed")

plot_fixed(ests_fixed, "stimtrns:taskpost", "open")
plot_fixed(ests_fixed, "stimtrns:taskpost", "closed")

# pre vs. post effect
plot_fixed(ests_fixed, "taskpost", "open")
plot_fixed(ests_fixed, "taskpost", "closed")

# block effect
plot_fixed(ests_fixed, "block", "open")
plot_fixed(ests_fixed, "block", "closed")

# stim effect
plot_fixed(ests_fixed, "stimtacs", "open")
plot_fixed(ests_fixed, "stimtdcs", "open")
plot_fixed(ests_fixed, "stimtrns", "open")

plot_fixed(ests_fixed, "stimtacs", "closed")
plot_fixed(ests_fixed, "stimtdcs", "closed")
plot_fixed(ests_fixed, "stimtrns", "closed")


# visualization ----

# first collapse across  blocks within pre and post
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
# interaction plots
ggplot(
  study_sum %>% filter(eyes == "closed", stim %in% c("sham", "tdcs")), 
  aes(task, M, group = stim, color = stim)
  ) +
  geom_point() +
  geom_path() +
  facet_wrap(~elec)

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
    ))

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
    electrode_size = .5, # size of electrode points
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
    elec_shape_col = p.value.sig, # null
    elec_shapes = c(1, 19), # 19
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
  facet_grid(stim~task)

psd_closed_sub_plot

# puts them together
library(patchwork)
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
    electrode_size = .5, # size of electrode points
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
    elec_shape_col = p.value.sig,
    elec_shapes = c(1, 19), #19
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
  facet_grid(stim~task)

psd_open_sub_plot

# puts them together
psd_open_fplot <- psd_open_plot | psd_open_sub_plot # final plot
psd_open_fplot

