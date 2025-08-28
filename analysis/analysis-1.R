# analysis-1.R
# Matt Kmiecik
# Started 06 April 2024

# Purpose: preliminary viz and analysis of LUC stim data

# source("prepro-to-r.r") # not run

# libraries ----
library(tidyverse)

# data ----
files <- as.list(dir(path = "../output/", pattern = "*.rda", full.names = TRUE))
walk(files, ~load(.x, .GlobalEnv))

# stimulation types:
# 1 == tdcs
# 3 == tacs
# 4 == trns

# Total Peak Alpha Frequency (PAF) and Center of Gravity (COG) ----

## subject-wise data
iaf_ss <- 
  iaf_res %>%  
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
contrasts(iaf_ss$stim)
contrasts(iaf_ss$eyes)
contrasts(iaf_ss$task)

## subject-wise summary
iaf_ss_sum <- 
  iaf_ss %>%
  pivot_longer(cols = c(paf, cog)) %>% 
  filter(!is.na(value)) %>% # removes missing data 
  group_by(ss, stim, eyes, task, name) %>%
  summarise(m = mean(value), n = n()) %>%
  ungroup()

# study summary
iaf_sum <- 
  iaf_ss_sum %>% 
  group_by(stim, eyes, task, name) %>%
  summarise(
    M = mean(m), 
    SD = sd(m), 
    N = n(), 
    SEM = SD/sqrt(N), 
    MOE = qt(.975, df = N - 1) * SEM
    ) %>%
  ungroup()

# plot distributions and mean plots here
ggplot(
  iaf_ss_sum %>% filter(name == "cog"), 
  aes(m, group = task, color = task, fill = task)
  ) +
  geom_density(alpha = 1/2) + 
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  labs(
    x = "Center of Gravity (Hz)", y = "Density", title = "Center of Gravity"
    ) +
  coord_cartesian(xlim = c(8, 13)) +
  facet_grid(stim~eyes) +
  theme_bw()

ggplot(
  iaf_ss_sum %>% filter(name == "paf"), 
  aes(m, group = task, color = task, fill = task)
) +
  geom_density(alpha = 1/2) + 
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  labs(
    x = "Peak Alpha Frequency (Hz)", y = "Density", title = "Peak Alpha Frequency"
    ) +
  coord_cartesian(xlim = c(8, 13)) +
  facet_grid(stim~eyes) +
  theme_bw()


# COG means
library(ggsci)
jco <- pal_jco("default")(10)

pd <- position_dodge(width = .4)
plot1 <- 
  ggplot(
  iaf_sum %>% filter(name == "cog"),
  aes(task, M, group = stim, color = stim)
  ) +
  geom_point(position = pd) +
  geom_path(position = pd) +
  geom_errorbar(aes(ymin = M-SEM, ymax = M+SEM), width = .2, position = pd) +
  labs(
    x = "Time", 
    y = "Mean Center of Gravity (Hz)", 
    caption = "SEM error bars.", 
    title = "Center of Gravity"
    ) +
  coord_cartesian(ylim = c(8, 12)) +
  facet_wrap(~eyes) +
  scale_color_jco() +
  theme_bw()
plot1
ggsave(
  filename = "../output/cog-mean-plot.png", 
  plot1, 
  width = 5, 
  height = 4.5, 
  units = "in"
  )

pd <- position_dodge(width = .4)
plot2 <- 
  ggplot(
  iaf_sum %>% filter(name == "paf"),
  aes(task, M, group = stim, color = stim)
) +
  geom_point(position = pd) +
  geom_path(position = pd) +
  geom_errorbar(aes(ymin = M-SEM, ymax = M+SEM), width = .2, position = pd) +
  labs(
    x = "Time", 
    y = "Mean Peak Alpha Frequency (Hz)", 
    caption = "SEM error bars.", 
    title = "Peak Alpha Frequency"
  ) +
  coord_cartesian(ylim = c(8, 12)) +
  facet_wrap(~eyes) +
  scale_color_jco() +
  theme_bw()
plot2
ggsave(
  filename = "../output/paf-mean-plot.png", 
  plot2, 
  width = 5, 
  height = 4.5, 
  units = "in"
)

# modeling ----
library(lme4); library(lmerTest) # libraries

## Eyes closed

# COG
cog_mod_closed <- 
  lmer(
    cog ~ 1 + block + task*stim + (1 | ss), 
    data = iaf_ss %>% filter(eyes == "closed")
    )
summary(cog_mod_closed)

# PAF
paf_mod_closed <- 
  lmer(
    paf ~ 1 + block + task*stim + (1 | ss), 
    data = iaf_ss %>% filter(eyes == "closed")
  )
summary(paf_mod_closed)

## Eyes open

# COG
cog_mod_open <- 
  lmer(
    cog ~ 1 + block + task*stim + (1 | ss), 
    data = iaf_ss %>% filter(eyes == "open")
  )
summary(cog_mod_open)

# PAF
paf_mod_open <- 
  lmer(
    paf ~ 1 + block + task*stim + (1 | ss), 
    data = iaf_ss %>% filter(eyes == "open")
  )
summary(paf_mod_open)