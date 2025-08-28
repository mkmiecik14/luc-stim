# psd-clip.R
# Matt Kmiecik
# 30 March 2025

# Purpose: reduces size of PSD file

# libraries ----
library(tidyverse)

load("../output/psd-res-2025-03-25.rda")
alpha <- seq(8, 12, .25) # alpha range
psd_res2 <- psd_res %>% filter(freq %in% alpha)

# writes out ----
saveRDS(psd_res2, file = "../output/psd-res-alpha-2025-03-25.rds")

