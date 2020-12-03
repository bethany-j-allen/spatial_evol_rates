#Bethany Allen   22nd September 2020
#Code to build simulated datasets of species richness within latitude bands for time interval t1
#  (compared to prior interval t0 and following interval t2)

#setwd("#####")

###Designate input parameters###

#Number of iterations
iterations <- 10000

#Number of latitudinal bins, likely 6 (i.e. 30deg bins)
nbins <- 6

#Limits of occurrences in each latitudinal bin
#Eyeballing the empirical data, opted to use random number up to 1000 for each bin
occs_range <- c(1:1000)       #Initial (t0)
add_occs_range <- c(1:500)   #t1 and t2

#Limits for number of species in initial and additional global species pools
#Ranges from roughly 100 to 800 for each clade in the empirical data
sp_range <- c(100:800)

#Range of extinction percentages to sample from for t1 and t2
ext_range <- c(1:100)

#Designate sampling levels
sample_pc <- c(1, 0.75, 0.5, 0.25)

#Create data frame to store results
results <- data.frame()
differences <- data.frame()


for (x in 1:iterations){

  #Task 1: Build t0
  
  #Designate number of species in the initial pool
  sp <- base::sample(sp_range, 1)

  #Designate number of occurrences in each latitudinal bin
  occs <- base::sample(occs_range, size = nbins, replace = T)

  #Assign each occurrence a species at random
  t0 <- list()          #Document IDs of all occurrences
  t0rich <- list()      #Document unique occurrences
  t0rich_count <- c()   #Create vector of species richness for each latitude band

  #For each latitude band
  for (a in 1:nbins){
    #Designate a species ID for each occurrence in the latitude band from the pool
    species_ids <- base::sample(sp, size = occs[a], replace = T)
    #Add bin to global list
    t0[[a]] <- species_ids
    t0rich[[a]] <- unique(species_ids)
    t0rich_count <- append(t0rich_count, length(unique(species_ids)))
  }

  #Task 2: Facilitate origination and extinction to create t1 and t2

  ###t0 -> t1###
  #Designate number of occurrences to survive (those that don't are local extinctions)
  surv_prop1 <- (base::sample(ext_range, size = nbins, replace = T)/100) #Generate a proportion
  surv_occs1 <- round((t0rich_count * surv_prop1), 0) #Convert proportion to whole occurrences

  #Designate new species pool for origination (those present in t0 plus preset additional number)
  #[This allows 'migration' while preventing taxa occurring in t0 and t2 but not t1, as would be
  #  expected under perfect sampling]
  #Designate number of species in the initial pool
  add_sp1 <- base::sample(sp_range, 1)
  new_sp1 <- append(unique(unlist(t0)), c((sp + 1):(sp + add_sp1)))

  #Desginate the number of new occurrences in each latitude bin
  new_occs1 <- base::sample(add_occs_range, size = nbins, replace = T)

  #Implement extinction and origination on t0
  t1 <- list()          #Document IDs of all occurrences
  t1rich <- list()      #Document unique occurrences
  t1rich_count <- c()   #Create vector of species richness for each latitude band

  for (b in 1:nbins){
    #Pull out one latitude bin
    focal_bin1 <- t0[[(b)]]
    #Remove 'local' extinctions by sampling survivors at random from t0
    focal_bin1 <- base::sample(focal_bin1, size = surv_occs1[b], replace = F)
    #Add origination
    new_occ_ids1 <- base::sample(new_sp1, size = new_occs1, replace = T)
    focal_bin1 <- append(focal_bin1, new_occ_ids1)
    #Add bin to global list
    t1[[b]] <- focal_bin1
    t1rich[[b]] <- unique(focal_bin1)
    t1rich_count <- append(t1rich_count, length(unique(focal_bin1)))
  }

  ###t1 -> t2###
  #Designate number of occurrences to survive (those that don't are local extinctions)
  surv_prop2 <- (base::sample(ext_range, size = nbins, replace = T)/100) #Generate a proportion
  surv_occs2 <- round((t1rich_count * surv_prop2), digits = 0)      #Convert proportion to whole occurrences

  #Designate new species pool for origination (those present in t1 plus preset additional number)
  add_sp2 <- base::sample(sp_range, 1)
  new_sp2 <- append(unique(unlist(t1)), c((sp + add_sp1 + 1):(sp + add_sp1 + add_sp2)))

  #Desginate the number of new occurrences in each latitude bin (used half of that in t0)
  new_occs2 <- base::sample(add_occs_range, size = nbins, replace = T)

  #Implement extinction and origination on t0
  t2 <- list()          #Document IDs of all occurrences
  t2rich <- list()      #Document unique occurrences
  t2rich_count <- c()   #Create vector of species richness for each latitude band

  for (d in 1:nbins){
    #Pull out one latitude bin
    focal_bin2 <- t1[[(d)]]
    #Remove 'local' extinctions by sampling survivors at random from t0
    focal_bin2 <- base::sample(focal_bin2, size = surv_occs2[d], replace = F)
    #Add origination
    new_occ_ids2 <- base::sample(new_sp2, size = new_occs2, replace = T)
    focal_bin2 <- append(focal_bin2, new_occ_ids2)
    #Add bin to global list
    t2[[d]] <- focal_bin2
    t2rich[[d]] <- unique(focal_bin2)
    t2rich_count <- append(t2rich_count, length(unique(focal_bin2)))
  }

  #See LDGs across the three time intervals
  #par(mfrow = c(3,1))
  #plot(t2rich_count, type = "o")
  #plot(t1rich_count, type = "o")
  #plot(t0rich_count, type = "o")

  #Task 3: Calculate origination and extinction rates globally for t1
  measured_rates <- data.frame()
  measured_diffs <- data.frame()

  #Produce global lists of unique occurrences (i.e. richness)
  t0_global <- unique(unlist(t0))
  t1_global <- unique(unlist(t1))
  t2_global <- unique(unlist(t2))

  #Raw (these counts include singletons)
  global_orig <- length(setdiff(t1_global, t0_global))   #Present in t1 but not in t0
  global_ext <- length(setdiff(t1_global, t2_global))    #Present in t1 but not in t2
  global_orig_p <- global_orig / length(t1_global)
  global_ext_p <- global_ext / length(t1_global)

  #Boundary crosser
  global_originations <- length(setdiff(intersect(t1_global, t2_global), t0_global))
                                                       #Present in t1 and t2 but not in t0
  global_extinctions <- length(setdiff(intersect(t0_global, t1_global), t2_global))
                                                       #Present in t0 and t1 but not in t2
  global_through <- length(intersect(t0_global, t2_global))
                                                       #This = three-timers in perfect sampling
  global_bc_orig <- log((global_through + global_originations)/global_through)
                                                       #Per-capita BC origination rate
  global_bc_ext <- log((global_through + global_extinctions)/global_through)
                                                       #Per-capita BC extinction rate
  global_bc_orig_diff <- global_orig_p - global_bc_orig
                                                       #Difference between raw and BC origination
  global_bc_ext_diff <- global_ext_p - global_bc_ext
                                                       #Difference between raw and BC extinction

  #Three-timer (as there is perfect sampling / no part-timers, the sampling probability is 1)
  global_2t_1 <- length(intersect(t0_global, t1_global)) #Present in t0 and t1 irrespective of t2
  global_2t_2 <- length(intersect(t1_global, t2_global)) #Present in t1 and t2 irrespective of t0
  global_3t_orig <- log(global_2t_2/global_through)      #3t origination rate
  global_3t_ext <- log(global_2t_1/global_through)       #3t extinction rate
  global_3t_orig_diff <- global_orig_p - global_3t_orig  #Difference between raw and 3t origination
  global_3t_ext_diff <- global_ext_p - global_3t_ext     #Difference between raw and 3t extinction

  #Add global rates to data frame
  global_rates <- c(x, "global", length(t1_global), global_orig, global_ext,
                    (round(c(global_orig_p, global_ext_p, global_bc_orig, global_bc_ext,
                             global_3t_orig, global_3t_ext), 3)))
  measured_rates <- rbind(measured_rates, global_rates)
  measured_diffs <- rbind(measured_diffs, c(x, "global", "origination", "boundary-crosser", round(global_bc_orig_diff, 3)))
  measured_diffs <- rbind(measured_diffs, c(x, "global", "extinction", "boundary-crosser", round(global_bc_ext_diff, 3)))
  measured_diffs <- rbind(measured_diffs, c(x, "global", "origination", "three-timer", round(global_3t_orig_diff, 3)))
  measured_diffs <- rbind(measured_diffs, c(x, "global", "extinction", "three-timer", round(global_3t_ext_diff, 3)))

  #Task 4: Calculate origination and extinction rates for individual latitude bands

  for (e in 1:nbins){
    #Pull out one latitude band
    focal_bin_t1 <- t1rich[[(e)]]
  
    #Calculate per-band origination and extinction rates
    #Raw (these counts include singletons)
    bin_orig <- length(setdiff(focal_bin_t1, t0_global))   #Present in t1 band but nowhere in t0
    bin_ext <- length(setdiff(focal_bin_t1, t2_global))    #Present in t1 band but nowhere in t2
    bin_orig_prop <- bin_orig / length(focal_bin_t1)
    bin_ext_prop <- bin_ext / length(focal_bin_t1)
  
    #Boundary crosser
    bin_originations <- length(setdiff(intersect(focal_bin_t1, t2_global), t0_global))
                                                         #Present in t1 and t2 but not in t0
    bin_extinctions <- length(setdiff(intersect(t0_global, focal_bin_t1), t2_global))
                                                         #Present in t0 and t1 but not in t2
    bin_through <- length(intersect(intersect(t0_global, focal_bin_t1), t2_global))
                                                         #Present in band for t1, globally for t0 & t2
    bin_bc_orig <- log((bin_through + bin_originations)/bin_through)
                                                         #Per-capita BC origination rate
    bin_bc_ext <- log((bin_through + bin_extinctions)/bin_through)
                                                         #Per-capita BC extinction rate
    bin_bc_orig_diff <- bin_orig_prop - bin_bc_orig
                                                         #Difference between raw and BC origination
    bin_bc_ext_diff <- bin_ext_prop - bin_bc_ext
                                                         #Difference between raw and BC extinction
  
    #Three-timer (as there is perfect sampling / no part-timers, the sampling probability is 1)
    bin_2t_1 <- length(intersect(t0_global, focal_bin_t1)) #Present in t0 and t0 irrespective of t2
    bin_2t_2 <- length(intersect(focal_bin_t1, t2_global)) #Present in t1 and t2 irrespective of t0
    bin_3t_orig <- log(bin_2t_2/bin_through)      #3t origination rate
    bin_3t_ext <- log(bin_2t_1/bin_through)       #3t extinction rate
    bin_3t_orig_diff <- bin_orig_prop - bin_3t_orig  #Difference between raw and 3t origination
    bin_3t_ext_diff <- bin_ext_prop - bin_3t_ext     #Difference between raw and 3t extinction
  
    #Save in a vector
    rates_vector <- c(x, e, length(focal_bin_t1), bin_orig, bin_ext, (round(c(bin_orig_prop,
                          bin_ext_prop, bin_bc_orig, bin_bc_ext, bin_3t_orig, bin_3t_ext), 3)))
    measured_rates <- rbind(measured_rates, rates_vector)
    
    measured_diffs <- rbind(measured_diffs, c(x, "lat_band", "origination", "boundary-crosser", round(bin_bc_orig_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, "lat_band", "extinction", "boundary-crosser", round(bin_bc_ext_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, "lat_band", "origination", "three-timer", round(bin_3t_orig_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, "lat_band", "extinction", "three-timer", round(bin_3t_ext_diff, 3)))
  }
  
  #Label columns in rates data frame
  colnames(measured_rates) <- c("iteration_no", "bin_no", "t1_n", "raw_origination", "raw_extinction",
                         "raw_origination_rate", "raw_extinction_rate", "BC_origination_pc",
                         "BC_extinction_pc", "tt_origination_rate", "tt_extinction_rate")
  colnames(measured_diffs) <- c("iteration_no", "bin_size", "rate", "method", "difference")
  results <- rbind(measured_rates, results)
  differences <- rbind(measured_diffs, differences)
}

#Summarise results
library(tidyverse)

differences$difference <- as.numeric(as.character(differences$difference))

ggplot(differences, aes(x = bin_size, y = difference, fill = rate)) + geom_boxplot() +
  facet_wrap(~method) + theme_classic()

