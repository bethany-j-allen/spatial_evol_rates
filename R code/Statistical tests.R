#Bethany Allen   21st October 2022
#Code to calculate test statistics on difference results

#setwd("#####")

library(tidyverse)
library(lsr)

differences <- read_csv("data/Simulation results/Sim_diffs_overall_main.csv")

differences <- filter(differences, sampling == "0") %>% filter(bin_size == "lat_band")

differences <- filter(differences, !is.nan(difference)) #not if paired

raw_orig <- filter(differences, method == "raw") %>% filter(rate == "origination")
raw_ext <- filter(differences, method == "raw") %>% filter(rate == "extinction")
BC_orig <- filter(differences, method == "boundary-crosser") %>% filter(rate == "origination")
BC_ext <- filter(differences, method == "boundary-crosser") %>% filter(rate == "extinction")
TT_orig <- filter(differences, method == "three-timer") %>% filter(rate == "origination")
TT_ext <- filter(differences, method == "three-timer") %>% filter(rate == "extinction")

t.test(x = raw_orig$difference, mu = 0)
t.test(x = raw_orig$difference, y = BC_orig$difference, paired = TRUE)

cohensD(x = raw_orig$difference, mu = 0)
cohensD(x = raw_orig$difference, y = BC_orig$difference, method = "paired")
