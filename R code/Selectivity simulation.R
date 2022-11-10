#Matthew Clapham

library(dplyr)
library(purrr)
library(ggplot2)
library(viridis)
library(patchwork)

#code to download brachiopod occurrences
#used as basis for simulated counts
#word_br <- read.csv("https://paleobiodb.org/data1.2/occs/list.txt?base_name=Rhynchonelliformea&interval=Wordian&show=class")
#genus_ct <- table(word_br$genus)
#genus_ct <- genus_ct[names(genus_ct) != ""]
#genus_ct <- sample(genus_ct, 400, replace=F)
#paste(sort(as.vector(genus_ct)), collapse=",")

#based on Wordian brachiopods from PBDB
genus_cts <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
               1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
               1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
               1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
               2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
               3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,
               5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,
               8,8,8,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,11,11,11,11,11,11,12,
               12,12,12,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,15,15,16,16,
               16,16,16,16,16,16,16,16,17,17,17,18,18,18,18,19,19,19,19,20,20,20,20,
               20,20,20,20,20,21,21,21,21,22,22,22,22,22,22,23,23,23,23,24,24,25,25,
               26,26,26,26,27,27,28,28,28,29,30,30,31,32,32,33,33,33,33,33,34,35,35,
               38,40,42,43,44,44,45,46,46,47,47,48,48,49,49,50,50,52,52,53,53,54,55,
               56,59,61,62,67,67,69,70,71,74,76,86,88,91,97,108,108,110,140,154)

#converts counts to proportion
genus_prob <- genus_cts / sum(genus_cts)

#sets up sequential letter names
genus_names_all <- c(LETTERS, c(t(outer(LETTERS, LETTERS, paste, sep = ""))))
genus_names <- genus_names_all[1:400]

#function to subsample extinction with varying selectivity
local.ext <- function(selection, size) {
  
  #creates vector based on selectivity values and uses Poisson distribution to convert to counts
  selectivity <- sapply(seq(0.69-selection, 0.69+selection, length.out=400), function(x) rpois(1, x))

  #discretizes Poisson counts and combines with genus names and counts
  genus_outcomes <- data.frame(genus_names, genus_cts, outcome = ifelse(selectivity > 0, "survive", "extinct"))
  
  #true extinction will be close to 50% but will vary, so is calculated here
  true_extinction <- mean(genus_outcomes$outcome == "extinct")
  
  #coefficient for log-odds ratio of the effect of abundance on extinction
  log_odds <- coef(glm(as.factor(outcome) ~ genus_cts, family="binomial", data=genus_outcomes))[2]
  
  #expand counts into list of occurrences for subsampling
  genus_ext <- genus_outcomes[rep(seq(nrow(genus_outcomes)), times=genus_outcomes$genus_cts),]
  
  #selects number of occurrences to sample, by converting size proportion into integer
  size_samp <- round(nrow(genus_ext) * size, 0)
  
  #selects occurrences and converts into counts of extinct and survive
  genus_res <- genus_ext %>% 
    slice_sample(n=size_samp) %>% 
    distinct(genus_names, .keep_all = T) %>% 
    count(outcome)
  
  #calculates extinction proportion (number of extinct / total number of taxa)
  extinction <- genus_res$n[1] / sum(genus_res$n)
  
  #compiles results into vector
  temp_res <- c(selectivity = log_odds, extinction_bias = extinction - true_extinction, size = size_samp)
  
  temp_res
}

#Subsampling and boundary-crosser extinction
#varies sampling of t0 (pre_size), t1 (size), t2 (post_size), and selectivity
bc.sampling <- function(size, pre_size, post_size, selection) {
  
  #creates vector based on selectivity values and uses Poisson distribution to convert to counts
  selectivity <- sapply(seq(0.69-selection, 0.69+selection, length.out=400), function(x) rpois(1, x))
  
  #discretizes Poisson counts and combines with genus names and counts
  genus_outcomes <- data.frame(genus_names, genus_cts, outcome = ifelse(selectivity > 0, "survive", "extinct"))
  
  #true extinction will be close to 50% but will vary, so is calculated here
  true_extinction <- mean(genus_outcomes$outcome == "extinct")
  
  #coefficient for log-odds ratio of the effect of abundance on extinction
  log_odds <- coef(glm(as.factor(outcome) ~ genus_cts, family="binomial", data=genus_outcomes))[2]

  #populates pre-extinction time bin t0
  pre_extinction <- genus_outcomes
  
  #populates post-extinction time bin t2
  post_extinction <- genus_outcomes[genus_outcomes$outcome == "survive",]

  #expand counts into list of occurrences for subsampling
  genus_ext <- genus_outcomes[rep(seq(nrow(genus_outcomes)), times=genus_outcomes$genus_cts),]
  
  #expands post-extinction counts into list of occurrences for subsampling
  post_sample <- post_extinction[rep(seq(nrow(post_extinction)), times=post_extinction$genus_cts),]
  
  #selects number of t1 occurrences to sample, by converting size proportion into integer
  size_samp <- round(nrow(genus_ext) * size, 0)
  
  #selects number of t0 occurrences to sample, by converting size proportion into integer
  pre_size <- round(nrow(genus_ext) * pre_size, 0)
  
  #selects number of t2 occurrences to sample, by converting size proportion into integer
  post_size_samp <- round(nrow(post_sample) * post_size, 0)
  
  #subsets t0 bin
  pre_taxa <- genus_ext %>% 
    slice_sample(n=pre_size) %>% 
    distinct(genus_names)
  
  #subsets t2 bin
  post_taxa <- post_sample %>% 
    slice_sample(n=post_size_samp) %>% 
    distinct(genus_names)
  
  #assigns extinct outcome to t1 list where names not in t2
  genus_ext$outcome[!genus_ext$genus_names %in% post_taxa$genus_names] <- "extinct"
  
  #selects boundary-crosser occurrences and converts into counts of extinct and survive
  genus_res <- genus_ext %>% 
    slice_sample(n=size_samp) %>% 
    filter(genus_names %in% pre_taxa$genus_names) %>% 
    distinct(genus_names, .keep_all = T) %>% 
    count(outcome)
  
  #calculates extinction proportion
  extinction <- genus_res$n[1] / sum(genus_res$n)
  
  #compiles results into vector
  temp_res <- c(selectivity = log_odds, extinction_bias = extinction - true_extinction, size = size_samp,
                pre_size = pre_size, post_size = post_size)
  
  temp_res
  
}

#Subsampling and three-timer extinction
#varies sampling of t0 (pre_size), t1 (size), t2 (post_size), and selectivity
tt.sampling <- function(size, pre_size, post_size, selection) {
  
  #creates vector based on selectivity values and uses Poisson distribution to convert to counts
  selectivity <- sapply(seq(0.69-selection, 0.69+selection, length.out=400), function(x) rpois(1, x))
  
  #discretizes Poisson counts and combines with genus names and counts
  genus_outcomes <- data.frame(genus_names, genus_cts, outcome = ifelse(selectivity > 0, "survive", "extinct"))
  
  #true extinction will be close to 50% but will vary, so is calculated here
  true_extinction <- mean(genus_outcomes$outcome == "extinct")
  
  #coefficient for log-odds ratio of the effect of abundance on extinction
  log_odds <- coef(glm(as.factor(outcome) ~ genus_cts, family="binomial", data=genus_outcomes))[2]
  
  #populates pre-extinction time bin t0
  pre_extinction <- genus_outcomes
  
  #populates post-extinction time bin t2
  post_extinction <- genus_outcomes[genus_outcomes$outcome == "survive",]
  
  #expand counts into list of occurrences for subsampling
  genus_ext <- genus_outcomes[rep(seq(nrow(genus_outcomes)), times=genus_outcomes$genus_cts),]
  
  #expands post-extinction counts into list of occurrences for subsampling
  post_sample <- post_extinction[rep(seq(nrow(post_extinction)), times=post_extinction$genus_cts),]
  
  #selects number of t1 occurrences to sample, by converting size proportion into integer
  size_samp <- round(nrow(genus_ext) * size, 0)
  
  #selects number of t0 occurrences to sample, by converting size proportion into integer
  pre_size <- round(nrow(genus_ext) * pre_size, 0)
  
  #selects number of t2 occurrences to sample, by converting size proportion into integer
  post_size_samp <- round(nrow(post_sample) * post_size, 0)
  
  #subsets t0 bin
  pre_taxa <- genus_ext %>% 
    slice_sample(n=pre_size) %>% 
    distinct(genus_names)
  
  #subsets t2 bin
  post_taxa <- post_sample %>% 
    slice_sample(n=post_size_samp) %>% 
    distinct(genus_names)
  
  #assigns extinct outcome to t1 list where names not in t2
  genus_ext$outcome[!genus_ext$genus_names %in% post_taxa$genus_names] <- "extinct"
  
  #selects two-timer/three-timer cohort taxa
  genus_focal <- genus_ext %>% 
    slice_sample(n=size_samp) %>% 
    filter(genus_names %in% pre_taxa$genus_names) %>% 
    distinct(genus_names, .keep_all = T) 
  
  #converts into counts of extinct and survive
  genus_res <- genus_focal %>% 
    count(outcome)
  
  #counts two-timer taxa
  ct_2t <- sum(genus_res$n)
  
  #counts three-timer taxa
  ct_3t <- genus_res$n[2]
  
  #three-timer taxa for t2 bin
  post_3t <- post_sample %>% 
    filter(genus_names %in% genus_focal$genus_names,
           genus_names %in% post_taxa$genus_names)
  
  #counts three-timers for t1, t2, t3
  ct_3t2 <- length(unique(post_3t$genus_names))
  
  #selects part-timer taxa
  part_timers <- post_sample %>% 
    filter(!genus_names %in% post_taxa$genus_names) %>% 
    filter(genus_names %in% genus_focal$genus_names)
  
  #counts part-timers
  ct_pt <- length(unique(part_timers$genus_names))
  
  #calculates probability of sampling
  ps <- ct_3t2 / (ct_3t2 + ct_pt)
  
  #calculates extinction proportion
  extinction <- 1 - ct_3t/(ps * ct_2t)
  
  #compiles results into vector
  temp_res <- c(selectivity=log_odds, extinction_bias=extinction - true_extinction, size=size_samp,
                pre_size=pre_size, post_size=post_size)
  
  temp_res
  
}

#Subsampling and range-through extinction
#varies sampling of t0 (pre_size), t1 (size), t2 (post_size), and selectivity
rt.sampling <- function(size, pre_size, post_size, selection) {
  
  #creates vector based on selectivity values and uses Poisson distribution to convert to counts
  selectivity <- sapply(seq(0.69-selection, 0.69+selection, length.out=400), function(x) rpois(1, x))
  
  #discretizes Poisson counts and combines with genus names and counts
  genus_outcomes <- data.frame(genus_names, genus_cts, outcome = ifelse(selectivity > 0, "survive", "extinct"))
  
  #true extinction will be close to 50% but will vary, so is calculated here
  true_extinction <- mean(genus_outcomes$outcome == "extinct")
  
  #coefficient for log-odds ratio of the effect of abundance on extinction
  log_odds <- coef(glm(as.factor(outcome) ~ genus_cts, family="binomial", data=genus_outcomes))[2]
  
  #populates post-extinction time bin t2
  post_extinction <- genus_outcomes[genus_outcomes$outcome == "survive",]
  
  #expand counts into list of occurrences for subsampling
  genus_ext <- genus_outcomes[rep(seq(nrow(genus_outcomes)), times=genus_outcomes$genus_cts),]
  
  #selects number of t1 occurrences to sample, by converting size proportion into integer
  size_samp <- round(nrow(genus_ext) * size, 0)
  
  #expands post-extinction counts into list of occurrences for subsampling
  post_sample <- post_extinction[rep(seq(nrow(post_extinction)), times=post_extinction$genus_cts),]
  
  #selects number of t0 occurrences to sample, by converting size proportion into integer
  pre_size <- round(nrow(genus_ext) * pre_size, 0)
  
  #selects number of t2 occurrences to sample, by converting size proportion into integer
  post_size_samp <- round(nrow(post_sample) * post_size, 0)
  
  #subsets t2 bin
  post_taxa <- post_sample %>% 
    slice_sample(n=post_size_samp) %>% 
    distinct(genus_names)
  
  #assigns extinct outcome to t1 list where names not in t2
  genus_ext$outcome[!genus_ext$genus_names %in% post_taxa$genus_names] <- "extinct"
  
  #selects occurrences and converts into counts of extinct and survive
  genus_res <- genus_ext %>% 
    slice_sample(n=size_samp) %>% 
    distinct(genus_names, .keep_all = T) %>% 
    count(outcome)
  
  #calculates extinction proportion
  extinction <- genus_res$n[1] / sum(genus_res$n)
  
  #compiles results into vector
  temp_res <- c(selectivity=log_odds, extinction_bias=extinction - true_extinction, size=size_samp,
                pre_size=pre_size, post_size=post_size)
  
  temp_res
  
}



#SUBSAMPLING OF A SINGLE INTERVAL

#Sets up test conditions of size and selectivity
test_conditions <- expand.grid(size = seq(0.1, 1, length.out=5), selection = seq(-0.6, 0.6, by=0.001))

#applies local extinction function to the range of test conditions
outcomes_syn <- map2_dfr(.x=test_conditions$selection, .y=test_conditions$size,  .f=local.ext)

#plots results
#uses the GAM smoothing for most results, except for perfect sampling where it doesn't work
#uses linear model for perfect sampling
ggplot(outcomes_syn, aes(selectivity.genus_cts, extinction_bias, color=as.factor(size))) +
  geom_point() +
  stat_smooth(data = filter(outcomes_syn, size == max(size)), method="lm", se=F, size=2) +
  stat_smooth(data = filter(outcomes_syn, size != max(size)), se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  annotate(geom = "text", x = max(outcomes_syn$selectivity.genus_cts) - 0.01,
           y = min(outcomes_syn$extinction_bias) - 0.01, label = "Common survive") +
  annotate(geom = "text", x = min(outcomes_syn$selectivity.genus_cts) + 0.01,
           y = min(outcomes_syn$extinction_bias) - 0.01, label = "Rare survive") +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  scale_color_manual(values = c("#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32"), 
                     labels = c("Smallest", "Small", "Medium", "Large", "Perfect"),
                     name="Bin sampling") +
  theme_classic()



#VARIABLE SAMPLING OF T0 BIN
test_conditions <- expand.grid(pre_size = seq(0.1, 1, length.out=5), selection = seq(-0.6, 0.6, by=0.001))

#runs for boundary-crosser extinction
outcomes_pre_bc <- map2_dfr(.x=test_conditions$pre_size, .y=test_conditions$selection, .f=bc.sampling, size = 1, post_size = 1)

#runs for three-timer extinction
outcomes_pre_tt <- map2_dfr(.x=test_conditions$pre_size, .y=test_conditions$selection, .f=tt.sampling, size = 1, post_size = 1)

#runs for range-through extinction
outcomes_pre_rt <- map2_dfr(.x=test_conditions$pre_size, .y=test_conditions$selection, .f=rt.sampling, size = 1, post_size = 1)

#boundary-crosser plot
bc_t0 <- ggplot(outcomes_pre_bc, aes(selectivity.genus_cts, extinction_bias, color=as.factor(pre_size))) +
  geom_point() +
  stat_smooth(data = filter(outcomes_pre_bc, pre_size == max(pre_size)), method="lm", se=F, size=2) +
  stat_smooth(data = filter(outcomes_pre_bc, pre_size != max(pre_size)), se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  annotate(geom = "text", x = max(outcomes_pre_bc$selectivity.genus_cts) - 0.01,
           y = -0.2, label = "Common survive") +
  annotate(geom = "text", x = min(outcomes_pre_bc$selectivity.genus_cts) + 0.01,
           y = -0.2, label = "Rare survive") +
  annotate(geom = "text", x = min(outcomes_pre_bc$selectivity.genus_cts) + 0.02,
           y = 0.2, label = "Boundary crosser") +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  ylim(c(-0.2, 0.2)) +
  scale_color_manual(values = c("#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#4a1486"),
                     labels = c("Smallest", "Small", "Medium", "Large", "Perfect"),
                     name = "Sampling t0") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))

#three-timer plot
tt_t0 <- ggplot(outcomes_pre_tt, aes(selectivity.genus_cts, extinction_bias, color=as.factor(pre_size))) +
  geom_point() +
  stat_smooth(data = filter(outcomes_pre_tt, pre_size == max(pre_size)), method="lm", se=F, size=2) +
  stat_smooth(data = filter(outcomes_pre_tt, pre_size != max(pre_size)), se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  annotate(geom = "text", x = max(outcomes_pre_tt$selectivity.genus_cts) - 0.01,
           y = -0.2, label = "Common survive") +
  annotate(geom = "text", x = min(outcomes_pre_tt$selectivity.genus_cts) + 0.01,
           y = -0.2, label = "Rare survive") +
  annotate(geom = "text", x = min(outcomes_pre_tt$selectivity.genus_cts) + 0.02,
           y = 0.2, label = "Three timer") +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  ylim(c(-0.2, 0.2)) +
  scale_color_manual(values = c("#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#8c2d04"),
                     labels = c("Smallest", "Small", "Medium", "Large", "Perfect"),
                     name = "Sampling t0") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))

#range-through plot
rt_t0 <- ggplot(outcomes_pre_rt, aes(selectivity.genus_cts, extinction_bias, color=as.factor(pre_size))) +
  geom_point() +
  stat_smooth(method="lm", se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  annotate(geom = "text", x = max(outcomes_pre_rt$selectivity.genus_cts) - 0.01,
           y = -0.2, label = "Common survive") +
  annotate(geom = "text", x = min(outcomes_pre_rt$selectivity.genus_cts) + 0.01,
           y = -0.2, label = "Rare survive") +
  annotate(geom = "text", x = min(outcomes_pre_rt$selectivity.genus_cts) + 0.02,
           y = 0.2, label = "Range through") +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  ylim(c(-0.2, 0.2)) +
  scale_color_manual(values = c("#fcc5c0", "#fa9fb5", "#f768a1", "#c51b8a", "#7a0177"),
                     labels = c("Smallest", "Small", "Medium", "Large", "Perfect"),
                     name = "Sampling t0") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))



#VARIABLE SAMPLING OF T2 BIN
test_conditions <- expand.grid(post_size = seq(0.1, 1, length.out=5), selection = seq(-0.6, 0.6, by=0.001))

#runs for boundary-crosser extinction
outcomes_post_bc <- map2_dfr(.x=test_conditions$post_size, .y=test_conditions$selection, .f=bc.sampling, size = 1, pre_size = 1)

#runs for three-timer extinction
outcomes_post_tt <- map2_dfr(.x=test_conditions$post_size, .y=test_conditions$selection, .f=tt.sampling, size = 1, pre_size = 1)

#runs for range-through extinction
outcomes_post_rt <- map2_dfr(.x=test_conditions$post_size, .y=test_conditions$selection, .f=rt.sampling, size = 1, pre_size = 1)

#boundary-crosser plot
bc_t2 <- ggplot(outcomes_post_bc, aes(selectivity.genus_cts, extinction_bias, color=as.factor(post_size))) +
  geom_point() +
  stat_smooth(data = filter(outcomes_post_bc, post_size == max(post_size)), method="lm", se=F, size=2) +
  stat_smooth(data = filter(outcomes_post_bc, post_size != max(post_size)), se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  annotate(geom = "text", x = max(outcomes_post_bc$selectivity.genus_cts) - 0.01,
           y = -0.018, label = "Common survive") +
  annotate(geom = "text", x = min(outcomes_post_bc$selectivity.genus_cts) + 0.01,
           y = -0.018, label = "Rare survive") +
  annotate(geom = "text", x = min(outcomes_post_bc$selectivity.genus_cts) + 0.02,
           y = 0.4, label = "Boundary crosser") +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  ylim(c(-0.02, 0.4)) +
  scale_color_manual(values = c("#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#4a1486"),
                     labels = c("Smallest", "Small", "Medium", "Large", "Perfect"),
                     name = "Sampling t2") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))

#three-timer plot
tt_t2 <- ggplot(outcomes_post_tt, aes(selectivity.genus_cts, extinction_bias, color=as.factor(post_size))) +
  geom_point() +
  stat_smooth(data = filter(outcomes_post_tt, post_size == max(post_size)), method="lm", se=F, size=2) +
  stat_smooth(data = filter(outcomes_post_tt, post_size != max(post_size)), se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  annotate(geom = "text", x = max(outcomes_post_tt$selectivity.genus_cts) - 0.01,
           y = -0.018, label = "Common survive") +
  annotate(geom = "text", x = min(outcomes_post_tt$selectivity.genus_cts) + 0.01,
           y = -0.018, label = "Rare survive") +
  annotate(geom = "text", x = min(outcomes_post_tt$selectivity.genus_cts) + 0.02,
           y = 0.4, label = "Three timer") +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  ylim(c(-0.02, 0.4)) +
  scale_color_manual(values = c("#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#8c2d04"),
                     labels = c("Smallest", "Small", "Medium", "Large", "Perfect"),
                     name = "Sampling t2") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))

#range-through plot
rt_t2 <- ggplot(outcomes_post_rt, aes(selectivity.genus_cts, extinction_bias, color=as.factor(post_size))) +
  geom_point() +
  stat_smooth(data = filter(outcomes_post_rt, post_size == max(post_size)), method="lm", se=F, size=2) +
  stat_smooth(data = filter(outcomes_post_rt, post_size != max(post_size)), se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  annotate(geom = "text", x = max(outcomes_post_rt$selectivity.genus_cts) - 0.01,
           y = -0.018, label = "Common survive") +
  annotate(geom = "text", x = min(outcomes_post_rt$selectivity.genus_cts) + 0.01,
           y = -0.018, label = "Rare survive") +
  annotate(geom = "text", x = min(outcomes_post_rt$selectivity.genus_cts) + 0.02,
           y = 0.4, label = "Range through") +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  ylim(c(-0.02, 0.4)) +
  scale_color_manual(values = c("#fcc5c0", "#fa9fb5", "#f768a1", "#c51b8a", "#7a0177"),
                     labels = c("Smallest", "Small", "Medium", "Large", "Perfect"),
                     name = "Sampling t2") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))


#arranges plots using patchwork package
(rt_t0 + rt_t2) /
(bc_t0 + bc_t2) /
(tt_t0 + tt_t2)





#VARIABLE SAMPLING OF T0, T1, T2 BIN, ALONG WITH VARYING SELECTIVITY 
test_conditions <- expand.grid(size = c(0.3, 0.6), post_size = c(0.3, 0.6, 1), pre_size = c(0.3, 0.6, 1), selection = seq(-0.6, 0.6, by=0.001))

#runs for boundary-crosser extinction
outcomes_bc <- pmap_dfr(.l = list(test_conditions$size, test_conditions$pre_size, test_conditions$post_size, test_conditions$selection), .f=bc.sampling)

#converts numeric size values into factor labels
#(this needs to be done because the pre-size and post-size differ between extinction methods)
outcomes_bc$t0_size <- as.factor(outcomes_bc$pre_size)
levels(outcomes_bc$t0_size) <- c("Small", "Moderate", "Perfect")
outcomes_bc$t0_size <- factor(outcomes_bc$t0_size,
                                          levels = c("Perfect", "Moderate", "Small"))

outcomes_bc$t2_size <- as.factor(outcomes_bc$post_size)
levels(outcomes_bc$t2_size) <- c("Small", "Moderate", "Perfect")

#adds a column with the method label
outcomes_bc$method <- "BC"

#runs for three-timer taxa
outcomes_tt <- pmap_dfr(.l = list(test_conditions$size, test_conditions$pre_size, test_conditions$post_size, test_conditions$selection), .f=tt.sampling)

#converts numeric size values into factor labels
outcomes_tt$t0_size <- as.factor(outcomes_tt$pre_size)
levels(outcomes_tt$t0_size) <- c("Small", "Moderate", "Perfect")
outcomes_tt$t0_size <- factor(outcomes_tt$t0_size,
                               levels = c("Perfect", "Moderate", "Small"))

outcomes_tt$t2_size <- as.factor(outcomes_tt$post_size)
levels(outcomes_tt$t2_size) <- c("Small", "Moderate", "Perfect")

#adds a column with the method label
outcomes_tt$method <- "TT"

#runs for range-through extinction
outcomes_rt <- pmap_dfr(.l = list(test_conditions$size, test_conditions$pre_size, test_conditions$post_size, test_conditions$selection), .f=rt.sampling)

#converts numeric size values into factor labels
outcomes_rt$t0_size <- as.factor(outcomes_rt$pre_size)
levels(outcomes_rt$t0_size) <- c("Small", "Moderate", "Perfect")
outcomes_rt$t0_size <- factor(outcomes_rt$t0_size,
                               levels = c("Perfect", "Moderate", "Small"))


outcomes_rt$t2_size <- as.factor(outcomes_rt$post_size)
levels(outcomes_rt$t2_size) <- c("Small", "Moderate", "Perfect")

#adds a column with the method label
outcomes_rt$method <- "RT"

#combines three extinction methods into a single data frame
outcomes_combined <- rbind(outcomes_bc, outcomes_tt, outcomes_rt)

#converts numeric size values into factor labels
outcomes_combined$t1_size <- as.factor(outcomes_combined$size)
levels(outcomes_combined$t1_size) <- c("Small", "Moderate")

#plots results
ggplot(outcomes_combined, aes(selectivity.genus_cts, extinction_bias,
                              color=method)) +
  stat_smooth(se=F, size=2) +
  geom_vline(xintercept = 0, color="gray", size=1, linetype=2) +
  geom_hline(yintercept = 0, color="gray", size=1, linetype=2) +
  facet_grid(t0_size ~ t2_size + t1_size, labeller = label_both) +
  xlab("Abundance selectivity") +
  ylab("Extinction proportion bias") +
  scale_color_manual(values = c("#807dba", "#f1b6da", "#f16913"),
                     name = "Method") +
  theme_bw()


