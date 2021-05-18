#Bethany Allen   22nd September 2020
#Code to plot outputs of the evolutionary rate simulation

#setwd("#####")

library(tidyverse)
library(pspearman)


###Reread results###
results <- read_csv("data/Sim_res_overall_main.csv")
abundances <- read_csv("data/Sim_abund_overall_main.csv")
differences <- read_csv("data/Sim_diffs_overall_main.csv")
sampling <- read_csv("data/Sim_samp_overall_main.csv")
extremes <- read_csv("data/Sim_extremes_overall_main.csv")
gradients <- read_csv("data/Sim_grads_overall_main.csv")
shifts <- read_csv("data/Sim_shifts_overall_main.csv")


###Summarise results###

#Examine distribution of global "true" origination and extinction rates
results_g <- filter(results, bin_no == "global") %>% filter(sampling == "100")
results_g$raw_origination_rate <- as.numeric(results_g$raw_origination_rate)
results_g$raw_extinction_rate <- as.numeric(results_g$raw_extinction_rate)

ggplot(results_g, aes(raw_origination_rate)) +
  geom_histogram(binwidth = 0.05, colour = "black", fill = "lightblue") + theme_classic()
ggplot(results_g, aes(raw_extinction_rate)) +
  geom_histogram(binwidth = 0.05, colour = "black", fill = "salmon") + theme_classic()


#Examine distribution of within-bin "true" origination and extinction rates
results_b <- filter(results, bin_no != "global") %>% filter(sampling == "100")
results_b$raw_origination_rate <- as.numeric(results_b$raw_origination_rate)
results_b$raw_extinction_rate <- as.numeric(results_b$raw_extinction_rate)

ggplot(results_b, aes(raw_origination_rate)) +
  geom_histogram(binwidth = 0.05, colour = "black", fill = "lightblue") + theme_classic()
ggplot(results_b, aes(raw_extinction_rate)) +
  geom_histogram(binwidth = 0.05, colour = "black", fill = "salmon") + theme_classic()


#Plot species abundance distributions
abundances_g <- filter(abundances, bin_no == "global")
abundances_g <- abundances_g[,4:ncol(abundances_g)]
abundances_g[] <- lapply(abundances_g, as.numeric)
mean_abundances_g <- colMeans(abundances_g)
max_abundances_g <- apply(abundances_g, 2, FUN = max)
min_abundances_g <- apply(abundances_g, 2, FUN = min)
abundances_g_sum <- cbind(c(1:ncol(abundances_g)), mean_abundances_g, max_abundances_g, min_abundances_g)
abundances_g_sum <- as.data.frame(abundances_g_sum)
ggplot(abundances_g_sum, aes(x = V1, y = mean_abundances_g, ymax = max_abundances_g, ymin = min_abundances_g)) +
  geom_ribbon(fill = "grey") + geom_line(size = 2, colour = "black") +
  labs(x = "Species identity", y = "Number of counts") +
  theme_classic()

abundances_b <- filter(abundances, bin_no != "global")
abundances_b <- abundances_b[,4:ncol(abundances_b)]
abundances_b[] <- lapply(abundances_b, as.numeric)
mean_abundances_b <- colMeans(abundances_b)
max_abundances_b <- apply(abundances_b, 2, FUN = max)
min_abundances_b <- apply(abundances_b, 2, FUN = min)
abundances_b_sum <- cbind(c(1:ncol(abundances_b)), mean_abundances_b, max_abundances_b, min_abundances_b)
abundances_b_sum <- as.data.frame(abundances_b_sum)
ggplot(abundances_b_sum, aes(x = V1, y = mean_abundances_b, ymax = max_abundances_b, ymin = min_abundances_b)) +
  geom_ribbon(fill = "grey") + geom_line(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 500)) +
  labs(x = "Species identity", y = "Number of counts") +
  theme_classic()


#Plot difference between "true" and measured rates at different sampling levels
sampled_g <- filter(differences, bin_size == "global")
sampled_g$occs <- as.numeric(as.character(sampled_g$occs))
sampled_g$difference <- as.numeric(as.character(sampled_g$difference))
sampled_g$sampling <- as.factor(sampled_g$sampling)

sampled_b <- filter(differences, bin_size == "lat_band")
sampled_b$occs <- as.numeric(as.character(sampled_b$occs))
sampled_b$difference <- as.numeric(as.character(sampled_b$difference))
sampled_b$sampling <- as.factor(sampled_b$sampling)

ggplot(sampled_g, aes(x = sampling, y = difference, fill = rate)) + geom_hline(aes(yintercept = 0)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  #scale_y_continuous(limits = c(-1, 1)) +
  theme_classic()
ggplot(sampled_b, aes(x = sampling, y = difference, fill = rate)) + geom_hline(aes(yintercept = 0)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  #scale_y_continuous(limits = c(-1, 1)) +
  theme_classic()


#Plot difference between "true" and measured rates post-sampling by sampling size or completeness
sampled_b <- filter(sampled_b, sampling == 0) %>%
  mutate(size_group = cut_width(occs, width = 100, boundary = 0))

sample_sizes <- filter(results, bin_no != "global") %>%
  select(iteration_no, bin_no, sampling, occs_t1, occs_t2) %>%
  pivot_wider(names_from = sampling, values_from = c(occs_t1, occs_t2))
sample_sizes[] <- lapply(sample_sizes, as.numeric)
sample_sizes <- mutate(sample_sizes, t1_compl = (occs_t1_0/occs_t1_100)) %>%
  mutate(t2_compl = (occs_t2_0/occs_t2_100))

t1_sampling <- rep(sample_sizes$t1_compl, each = 3)
t2_sampling <- rep(sample_sizes$t2_compl, each = 3)
sampling_vector <- c(rbind(t2_sampling, t1_sampling))
sampled_b <- mutate(sampled_b, sampling_prop = sampling_vector)
sampled_b <- filter(sampled_b, sampling_prop != "Inf")
sampled_b <- mutate(sampled_b, sampling_group = cut_width(sampling_prop, width = 0.2, boundary = 0))

ggplot(sampled_b, aes(x = size_group, y = difference, counts_cut_width, fill = rate)) + geom_hline(aes(yintercept = 0)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  #scale_y_continuous(limits = c(-1, 1)) +
  theme_classic()

ggplot(sampled_b, aes(x = sampling_group, y = difference, counts_cut_width, fill = rate)) + geom_hline(aes(yintercept = 0)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  #scale_y_continuous(limits = c(-1, 1)) +
  theme_classic()


#Plot the relationship between sample size and rate differences
sampled_g_o <- filter(sampled_g, rate == "origination")
sampled_g_e <- filter(sampled_g, rate == "extinction")

sampled_b_o <- filter(sampled_b, rate == "origination")
sampled_b_e <- filter(sampled_b, rate == "extinction")

ggplot(sampled_g_o, aes(x = occs, y = difference)) + geom_hline(aes(yintercept = 0)) +
  geom_point(colour = "lightblue") + facet_wrap(~method + sampling) + theme_classic()
ggplot(sampled_g_e, aes(x = occs, y = difference)) + geom_hline(aes(yintercept = 0)) +
  geom_point(colour = "salmon") + facet_wrap(~method + sampling) + theme_classic()

ggplot(sampled_b_o, aes(x = occs, y = difference)) + geom_hline(aes(yintercept = 0)) +
  geom_point(colour = "lightblue") + facet_wrap(~method + sampling) + theme_classic()
ggplot(sampled_b_e, aes(x = occs, y = difference)) + geom_hline(aes(yintercept = 0)) +
  geom_point(colour = "salmon") + facet_wrap(~method + sampling) + theme_classic()


#Plot number of iterations with matching maximum and minimum rates
min_same <- extremes %>% group_by(sampling, rate, method) %>% count(min_match) %>% filter(min_match == T)
max_same <- extremes %>% group_by(sampling, rate, method) %>% count(max_match) %>% filter(max_match == T)

ggplot(min_same, aes(x = rate, y = n, fill = rate)) + geom_hline(aes(yintercept = iterations)) +
  geom_col() + facet_wrap(~method + sampling) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()

ggplot(max_same, aes(x = rate, y = n, fill = rate)) + geom_hline(aes(yintercept = iterations)) +
  geom_col() + facet_wrap(~method + sampling) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()


#Plot Pearson's p-values comparing linear correlation of latitude bins
gradients$t.p_value <- as.numeric(as.character(gradients$t.p_value))
gradients$cor <- as.numeric(as.character(gradients$cor))

ggplot(gradients, aes(x = sampling, y = t.p_value, fill = rate)) + geom_hline(aes(yintercept = 0.05)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()
ggplot(gradients, aes(x = sampling, y = cor, fill = rate)) + geom_hline(aes(yintercept = 0)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()


#Plot Spearman's p-values comparing rank order of latitude bins
gradients$S.p_value <- as.numeric(as.character(gradients$S.p_value))
gradients$rho <- as.numeric(as.character(gradients$rho))

ggplot(gradients, aes(x = sampling, y = S.p_value, fill = rate)) + geom_hline(aes(yintercept = 0.05)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()
ggplot(gradients, aes(x = sampling, y = rho, fill = rate)) + geom_hline(aes(yintercept = 0)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()


#Plot number of bins which overestimate rates at different sampling levels
shifts <- filter(shifts, sampling != "100")
shifts$bins_over <- as.numeric(as.character(shifts$bins_over))

ggplot(shifts, aes(x = bins_over, fill = rate)) +
  geom_bar() + facet_wrap(~method + rate) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()


#Plot sampling level estimated by the 3t method
sampling <- pivot_longer(sampling, cols = c("estimated_o_t1", "estimated_e_t2"), names_to = "rate")
sampling$value <- as.numeric(as.character(sampling$value))
sampling$sampled <- as.character(sampling$sampled/100)
ggplot(sampling, aes(x = sampled, y = value, fill = rate)) + geom_boxplot() + 
  scale_fill_manual(values = c("salmon", "lightblue")) + theme_classic()
