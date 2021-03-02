#Bethany Allen   22nd September 2020
#Code to build simulated datasets of species richness within latitude bands
#Currently uses four time bins (t0, t1, t2, t3) and calculates rates across t1-t2 boundary

#setwd("#####")

library(tidyverse)
library(pspearman)

###Designate input parameters###

#Number of iterations
iterations <- 1000

#Number of latitudinal bins (e.g. 6 -> 30deg bins)
nbins <- 6

#Do occurrences vary against fixed sampling rates, or are occurrences fixed while sampling varies?
to_vary <- "occurrences"

#Set occurrence numbers
#If occurrences vary, this is the upper limit of occurrences in each latitudinal bin
#If sampling varies, this is the fixed number of occurrences used
occs_range <- 1000       #Initial (t0)
add_occs_range <- 500   #t1 and t2

#If occurrences vary, sampling levels are 25%, 50% and 75% of occurrences
#If sampling varies, it ranges between 0% and 100% of occurrences for each latitude bin

#Limits for number of species in initial and additional global species pools
#Ranges from roughly 100 to 800 for each clade in the empirical data
sp_range <- c(500:2000)        #Initial (t0)
add_sp_range <- c(100:1500)    #t1 and t2

#Range of survival (not extinction) percentages to sample from for t1 and t2
ext_range <- seq(from = 0, to = 50, by = 0.01)

#Create data frames to store results
results <- data.frame(); differences <- data.frame(); sampling <- data.frame()
extremes <- data.frame(); gradients <- data.frame(); shifts <- data.frame()

#Add a progress bar to show simulation completion
pb <- txtProgressBar(min = 0, max = iterations, initial = 0, style = 3)


for (x in 1:iterations){

  #Task 1: Build t0
  
  #Designate number of species in the initial pool
  sp <- base::sample(sp_range, 1)

  #Designate number of occurrences in each latitudinal bin
  if (to_vary == "occurrences") {occs <- base::sample(c(0:occs_range), size = nbins, replace = T)} else
    if (to_vary == "sampling") {occs <- rep(occs_range, nbins)}

  #Create lists to document IDs of all occurrences, and samples of 75%, 50% and 25% of occurrences
  t0_100 <- list()
  if (to_vary == "occurrences") {t0_75 <- list(); t0_50 <- list(); t0_25 <- list()} else
    if (to_vary == "sampling") {t0_0 <- list()}

  #For each latitude band
  for (a in 1:nbins){
    #Designate a species ID for each occurrence in the latitude band from the pool
    species_ids <- base::sample(sp, size = occs[a], replace = T)
    #Add bin to global list
    t0_100[[a]] <- species_ids
    if (to_vary == "occurrences") {t0_75[[a]] <- species_ids[1:round((0.75 * length(species_ids)), 0)];
      t0_50[[a]] <- species_ids[1:round((0.5 * length(species_ids)), 0)];
      t0_25[[a]] <- species_ids[1:round((0.25 * length(species_ids)), 0)]} else
      if (to_vary == "sampling") {t0_0[[a]] <- species_ids[1:base::sample(length(species_ids), 1)]}
  }

  #Task 2: Facilitate origination and extinction to create t1, t2 and t3

  ###t0 -> t1###
  #Designate number of occurrences to survive (those that don't are local extinctions)
  surv_prop1 <- (base::sample(ext_range, size = nbins, replace = T)/100) #Generate a proportion
  surv_occs1 <- round((lengths(t0_100) * surv_prop1), 0) #Convert proportion to whole occurrences

  #Designate new species pool for origination (those present in t0 plus preset additional number)
  #[This allows 'migration' while preventing taxa occurring in t0 and t2 but not t1, as would be
  #  expected under perfect sampling]
  #Designate number of species in the initial pool
  add_sp1 <- base::sample(add_sp_range, 1)
  new_sp1 <- append(unique(unlist(t0_100)), c((sp + 1):(sp + add_sp1)))

  #Designate the number of new occurrences in each latitude bin
  if (to_vary == "occurrences") {new_occs1 <- base::sample(c(0:add_occs_range), size = nbins, replace = T)} else
    if (to_vary == "sampling") {new_occs1 <- rep(add_occs_range, nbins)}

  #Implement extinction and origination on t0
  #Create lists to document IDs of all occurrences, and samples of 75%, 50% and 25% of occurrences
  t1_100 <- list()
  if (to_vary == "occurrences") {t1_75 <- list(); t1_50 <- list(); t1_25 <- list()} else
    if (to_vary == "sampling") {t1_0 <- list()}

  for (b in 1:nbins){
    #Pull out one latitude bin
    focal_bin1 <- t0_100[[(b)]]
    #Remove 'local' extinctions by sampling survivors at random from t0
    focal_bin1 <- base::sample(focal_bin1, size = surv_occs1[b], replace = F)
    #Add origination
    new_occ_ids1 <- base::sample(new_sp1, size = new_occs1[b], replace = T)
    focal_bin1 <- append(focal_bin1, new_occ_ids1)
    #Add bin to global list
    t1_100[[b]] <- focal_bin1
    if (to_vary == "occurrences") {t1_75[[b]] <- focal_bin1[1:round((0.75 * length(focal_bin1)), 0)];
    t1_50[[b]] <- focal_bin1[1:round((0.5 * length(focal_bin1)), 0)];
    t1_25[[b]] <- focal_bin1[1:round((0.25 * length(focal_bin1)), 0)]} else
      if (to_vary == "sampling") {t1_0[[b]] <- focal_bin1[1:base::sample(length(focal_bin1), 1)]}
  }

  ###t1 -> t2###
  #Designate number of occurrences to survive (those that don't are local extinctions)
  surv_prop2 <- (base::sample(ext_range, size = nbins, replace = T)/100) #Generate a proportion
  surv_occs2 <- round((lengths(t1_100) * surv_prop2), digits = 0)      #Convert proportion to whole occurrences

  #Designate new species pool for origination (those present in t1 plus preset additional number)
  add_sp2 <- base::sample(add_sp_range, 1)
  new_sp2 <- append(unique(unlist(t1_100)), c((sp + add_sp1 + 1):(sp + add_sp1 + add_sp2)))

  #Designate the number of new occurrences in each latitude bin
  if (to_vary == "occurrences") {new_occs2 <- base::sample(c(0:add_occs_range), size = nbins, replace = T)} else
    if (to_vary == "sampling") {new_occs2 <- rep(add_occs_range, nbins)}

  #Implement extinction and origination on t1
  #Create lists to document IDs of all occurrences, and samples of 75%, 50% and 25% of occurrences
  t2_100 <- list()
  if (to_vary == "occurrences") {t2_75 <- list(); t2_50 <- list(); t2_25 <- list()} else
    if (to_vary == "sampling") {t2_0 <- list()}

  for (d in 1:nbins){
    #Pull out one latitude bin
    focal_bin2 <- t1_100[[(d)]]
    #Remove 'local' extinctions by sampling survivors at random from t1
    focal_bin2 <- base::sample(focal_bin2, size = surv_occs2[d], replace = F)
    #Add origination
    new_occ_ids2 <- base::sample(new_sp2, size = new_occs2[d], replace = T)
    focal_bin2 <- append(focal_bin2, new_occ_ids2)
    #Add bin to global list
    t2_100[[d]] <- focal_bin2
    if (to_vary == "occurrences") {t2_75[[d]] <- focal_bin2[1:round((0.75 * length(focal_bin2)), 0)];
    t2_50[[d]] <- focal_bin2[1:round((0.5 * length(focal_bin2)), 0)];
    t2_25[[d]] <- focal_bin2[1:round((0.25 * length(focal_bin2)), 0)]} else
      if (to_vary == "sampling") {t2_0[[d]] <- focal_bin2[1:base::sample(length(focal_bin2), 1)]}
  }
  
  ###t2 -> t3###
  #Designate number of occurrences to survive (those that don't are local extinctions)
  surv_prop3 <- (base::sample(ext_range, size = nbins, replace = T)/100) #Generate a proportion
  surv_occs3 <- round((lengths(t2_100) * surv_prop3), digits = 0) #Convert proportion to whole occurrences
  
  #Designate new species pool for origination (those present in t1 plus preset additional number)
  add_sp3 <- base::sample(add_sp_range, 1)
  new_sp3 <- append(unique(unlist(t2_100)), c((sp + add_sp1 + add_sp2 + 1):(sp + add_sp1 + add_sp2 + add_sp3)))
  
  #Designate the number of new occurrences in each latitude bin
  if (to_vary == "occurrences") {new_occs3 <- base::sample(c(0:add_occs_range), size = nbins, replace = T)} else
    if (to_vary == "sampling") {new_occs3 <- rep(add_occs_range, nbins)}
  
  #Implement extinction and origination on t1
  #Create lists to document IDs of all occurrences, and samples of 75%, 50% and 25% of occurrences
  t3_100 <- list()
  if (to_vary == "occurrences") {t3_75 <- list(); t3_50 <- list(); t3_25 <- list()} else
    if (to_vary == "sampling") {t3_0 <- list()}
  
  for (e in 1:nbins){
    #Pull out one latitude bin
    focal_bin3 <- t2_100[[(e)]]
    #Remove 'local' extinctions by sampling survivors at random from t2
    focal_bin3 <- base::sample(focal_bin3, size = surv_occs3[e], replace = F)
    #Add origination
    new_occ_ids3 <- base::sample(new_sp3, size = new_occs3[e], replace = T)
    focal_bin3 <- append(focal_bin3, new_occ_ids3)
    #Add bin to global list
    t3_100[[e]] <- focal_bin3
    if (to_vary == "occurrences") {t3_75[[e]] <- focal_bin3[1:round((0.75 * length(focal_bin3)), 0)];
    t3_50[[e]] <- focal_bin3[1:round((0.5 * length(focal_bin3)), 0)];
    t3_25[[e]] <- focal_bin3[1:round((0.25 * length(focal_bin3)), 0)]} else
      if (to_vary == "sampling") {t3_0[[e]] <- focal_bin3[1:base::sample(length(focal_bin3), 1)]}
  }
  
  #Task 3: Calculate origination and extinction rates, using each sampling level, at global
  #  scale and for individual latitude bands
  
  #Create data frames to store results in
  measured_rates <- data.frame(); measured_diffs <- data.frame(); sampling_3t_est <- data.frame()
  
  #Designate sampling levels, starting with 100%
  if (to_vary == "occurrences") {sample_pc <- c(100, 75, 50, 25)} else
    if (to_vary == "sampling") {sample_pc <- c(100, 0)}
  
  for (f in 1:length(sample_pc)){

    #Produce global lists of unique occurrences (i.e. richness) for the given sampling level
    t0_global <- unique(unlist(eval(parse(text = paste0("t0_", sample_pc[f])))))
    t1_global <- unique(unlist(eval(parse(text = paste0("t1_", sample_pc[f])))))
    t2_global <- unique(unlist(eval(parse(text = paste0("t2_", sample_pc[f])))))
    t3_global <- unique(unlist(eval(parse(text = paste0("t3_", sample_pc[f])))))
    
    #Produce global lists of occurrence numbers for the given sampling level
    t0_occs <- lengths(eval(parse(text = paste0("t0_", sample_pc[f]))))
    t1_occs <- lengths(eval(parse(text = paste0("t1_", sample_pc[f]))))
    t2_occs <- lengths(eval(parse(text = paste0("t2_", sample_pc[f]))))
    t3_occs <- lengths(eval(parse(text = paste0("t3_", sample_pc[f]))))
  
    #Raw (these counts include singletons)
    global_orig <- length(setdiff(t2_global, t1_global))   #Present in t2 but not in t1
    global_ext <- length(setdiff(t1_global, t2_global))    #Present in t1 but not in t2
    global_orig_p <- global_orig / length(t2_global)
    global_ext_p <- global_ext / length(t1_global)
  
    #Store 100% sampled values as the benchmark for comparison
    if (f == 1){global_true_orig <- global_orig_p; global_true_ext <- global_ext_p}
    
    #Compare sampled raw counts to true value
    global_orig_diff <- global_orig_p - global_true_orig
    global_ext_diff <- global_ext_p - global_true_ext

    #Boundary crosser
    global_originations <- length(setdiff(intersect(t2_global, t3_global), t1_global))
                                                       #Present in t2 and t3 but not in t1
    global_extinctions <- length(setdiff(intersect(t0_global, t1_global), t2_global))
                                                       #Present in t0 and t1 but not in t2
    global_through_o <- length(intersect(t1_global, t3_global))
    global_through_e <- length(intersect(t0_global, t2_global))
                                                       #This = three-timers in perfect sampling
    global_bc_orig <- global_originations/(global_through_o + global_originations)
                                                       #BC origination proportion
    global_bc_ext <- global_extinctions/(global_through_e + global_extinctions)
                                                       #BC extinction proportion
    
    #Compare sampled BC proportions to 100% value
    global_bc_orig_diff <- global_bc_orig - global_true_orig
    global_bc_ext_diff <- global_bc_ext - global_true_ext

    #Three-timer
    global_2t_o <- length(intersect(t2_global, t3_global)) #Present in t2 and t3 irrespective of t1
    global_2t_e <- length(intersect(t0_global, t1_global)) #Present in t0 and t1 irrespective of t2
    
    global_3t_1 <- length(intersect(intersect(t0_global, t1_global), t2_global)) #Present in t0-2
    global_3t_2 <- length(intersect(intersect(t1_global, t2_global), t3_global)) #Present in t1-3

    global_pt_1 <- length(setdiff(intersect(t0_global, t2_global), t1_global)) #Ghost ranges for t1
    global_pt_2 <- length(setdiff(intersect(t1_global, t3_global), t2_global)) #Ghost ranges for t2
    
    global_sampling_o <- global_3t_1/(global_3t_1 + global_pt_1) #Sampling completeness est. for t1
    global_sampling_e <- global_3t_2/(global_3t_2 + global_pt_2) #Sampling completeness est. for t2
    
    global_3t_orig <- 1 - (global_3t_2/(global_sampling_o*global_2t_o)) #3t origination proportion
    global_3t_ext <- 1 - (global_3t_1/(global_sampling_e*global_2t_e))  #3t extinction proportion
    
    #Compare sampled 3T proportions to 100% value
    global_3t_orig_diff <- global_3t_orig - global_true_orig
    global_3t_ext_diff <- global_3t_ext - global_true_ext
    
    #Add global rates to data frame
    global_rates <- c(x, "global", sample_pc[f], sum(t1_occs), sum(t2_occs), length(t1_global), length(t2_global),
                      global_orig, global_ext, (round(c(global_orig_p, global_ext_p, global_bc_orig,
                                                        global_bc_ext, global_3t_orig, global_3t_ext), 3)))
    measured_rates <- rbind(measured_rates, global_rates)
    sampling_3t_est <- rbind(sampling_3t_est, c(x, sample_pc[f], global_sampling_o, global_sampling_e))
    measured_diffs <- rbind(measured_diffs, c(x, sample_pc[f], sum(t2_occs), length(t2_global), "global", "origination", "raw", round(global_orig_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, sample_pc[f], sum(t1_occs), length(t1_global), "global", "extinction", "raw", round(global_ext_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, sample_pc[f], sum(t2_occs), length(t2_global), "global", "origination", "boundary-crosser", round(global_bc_orig_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, sample_pc[f], sum(t1_occs), length(t1_global), "global", "extinction", "boundary-crosser", round(global_bc_ext_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, sample_pc[f], sum(t2_occs), length(t2_global), "global", "origination", "three-timer", round(global_3t_orig_diff, 3)))
    measured_diffs <- rbind(measured_diffs, c(x, sample_pc[f], sum(t1_occs), length(t1_global), "global", "extinction", "three-timer", round(global_3t_ext_diff, 3)))
  }

  #Add column names to sampling estimate data frame
  colnames(sampling_3t_est) <- c("iteration_no", "sampled", "estimated_o_t1", "estimated_e_t2")
  
  #Calculate rates for each latitude bin
  for (e in 1:nbins){
    for (g in 1:length(sample_pc)){
      
      #Pull out unique occurrences from one latitude band
      focal_bin_t1 <- unique(eval(parse(text = paste0("t1_", sample_pc[g], "[[(", e, ")]]"))))
      focal_bin_t2 <- unique(eval(parse(text = paste0("t2_", sample_pc[g], "[[(", e, ")]]"))))
      t0_global <- unique(unlist(eval(parse(text = paste0("t0_", sample_pc[g])))))
      t1_global <- unique(unlist(eval(parse(text = paste0("t1_", sample_pc[g])))))
      t2_global <- unique(unlist(eval(parse(text = paste0("t2_", sample_pc[g])))))
      t3_global <- unique(unlist(eval(parse(text = paste0("t3_", sample_pc[g])))))
      
      focal_bin_t1_occs <- length(eval(parse(text = paste0("t1_", sample_pc[g], "[[(", e, ")]]"))))
      focal_bin_t2_occs <- length(eval(parse(text = paste0("t2_", sample_pc[g], "[[(", e, ")]]"))))
  
      #Calculate per-band origination and extinction rates
      #Raw (these counts include singletons)
      bin_orig <- length(setdiff(focal_bin_t2, t1_global))   #Present in t2 band but nowhere in t1
      bin_ext <- length(setdiff(focal_bin_t1, t2_global))    #Present in t1 band but nowhere in t2
      bin_orig_prop <- bin_orig / length(focal_bin_t2)
      bin_ext_prop <- bin_ext / length(focal_bin_t1)
    
      #Store 100% sampled values as the benchmark for comparison
      if (g == 1){bin_true_orig <- bin_orig_prop; bin_true_ext <- bin_ext_prop}
      
      #Compare sampled raw counts to true value
      bin_orig_diff <- bin_orig_prop - bin_true_orig
      bin_ext_diff <- bin_ext_prop - bin_true_ext
  
      #Boundary crosser
      bin_originations <- length(setdiff(intersect(focal_bin_t2, t3_global), t1_global))
                                                         #Present in t2 and t3 but not in t1
      bin_extinctions <- length(setdiff(intersect(t0_global, focal_bin_t1), t2_global))
                                                         #Present in t0 and t1 but not in t2
      bin_through_o <- length(intersect(intersect(t1_global, focal_bin_t2), t3_global))
                                                         #Present in band for t2, globally for t1 & t3
      bin_through_e <- length(intersect(intersect(t0_global, focal_bin_t1), t2_global))
                                                         #Present in band for t1, globally for t0 & t2
      bin_bc_orig <- bin_originations/(bin_through_o + bin_originations)
                                                         #BC origination proportion
      bin_bc_ext <- bin_extinctions/(bin_through_e + bin_extinctions)
                                                         #BC extinction proportion
      
      #Compare sampled BC rates to 100% value
      bin_bc_orig_diff <- bin_bc_orig - bin_true_orig
      bin_bc_ext_diff <- bin_bc_ext - bin_true_ext
  
      #Three-timer - uses global sampling probability
      bin_2t_o <- length(intersect(focal_bin_t2, t3_global)) #Present in t2 and t3 irrespective of t1
      bin_2t_e <- length(intersect(t0_global, focal_bin_t1)) #Present in t0 and t1 irrespective of t2
      
      bin_3t_1 <- length(intersect(intersect(t0_global, focal_bin_t1), t2_global)) #Present in t0-2
      bin_3t_2 <- length(intersect(intersect(t1_global, focal_bin_t2), t3_global)) #Present in t1-3
      
      global_est_sampling_o <- sampling_3t_est$estimated_o_t1[g]    #Pull the global sampling estimate
      global_est_sampling_e <- sampling_3t_est$estimated_e_t2[g]    #Pull the global sampling estimate
      
      bin_3t_orig <- 1 - (bin_3t_2/(global_est_sampling_o*bin_2t_o))      #3t origination proportion
      bin_3t_ext <- 1 - (bin_3t_1/(global_est_sampling_e*bin_2t_e))       #3t extinction proportion
      
      #Compare sampled 3T rates to 100% value
      bin_3t_orig_diff <- bin_3t_orig - bin_true_orig
      bin_3t_ext_diff <- bin_3t_ext - bin_true_ext
  
      #Save in a vector
      rates_vector <- c(x, e, sample_pc[g], focal_bin_t1_occs, focal_bin_t2_occs,
                        length(focal_bin_t1), length(focal_bin_t2), bin_orig, bin_ext,
                        (round(c(bin_orig_prop, bin_ext_prop, bin_bc_orig, bin_bc_ext, bin_3t_orig, bin_3t_ext), 3)))
      measured_rates <- rbind(measured_rates, rates_vector)
      measured_diffs <- rbind(measured_diffs, c(x, sample_pc[g], focal_bin_t2_occs, length(focal_bin_t2), "lat_band", "origination", "raw", round(bin_orig_diff, 3)))
      measured_diffs <- rbind(measured_diffs, c(x, sample_pc[g], focal_bin_t1_occs, length(focal_bin_t1), "lat_band", "extinction", "raw", round(bin_ext_diff, 3)))
      measured_diffs <- rbind(measured_diffs, c(x, sample_pc[g], focal_bin_t2_occs, length(focal_bin_t2), "lat_band", "origination", "boundary-crosser", round(bin_bc_orig_diff, 3)))
      measured_diffs <- rbind(measured_diffs, c(x, sample_pc[g], focal_bin_t1_occs, length(focal_bin_t1), "lat_band", "extinction", "boundary-crosser", round(bin_bc_ext_diff, 3)))
      measured_diffs <- rbind(measured_diffs, c(x, sample_pc[g], focal_bin_t2_occs, length(focal_bin_t2), "lat_band", "origination", "three-timer", round(bin_3t_orig_diff, 3)))
      measured_diffs <- rbind(measured_diffs, c(x, sample_pc[g], focal_bin_t1_occs, length(focal_bin_t1), "lat_band", "extinction", "three-timer", round(bin_3t_ext_diff, 3)))
    }
  }
  
  #Label columns in rates data frame
  colnames(measured_rates) <- c("iteration_no", "bin_no", "sampling", "occs_t1", "occs_t2", "richness_t1",
                                "richness_t2", "raw_origination", "raw_extinction", "raw_origination_rate",
                                "raw_extinction_rate", "BC_origination_pc", "BC_extinction_pc",
                                "tt_origination_rate", "tt_extinction_rate")
  colnames(measured_diffs) <- c("iteration_no", "sampling", "occs", "richness", "bin_size", "rate", "method", "difference")
  results <- rbind(results, measured_rates)
  differences <- rbind(differences, measured_diffs)
  sampling <- rbind(sampling, sampling_3t_est)
  
  #Task 4: Compare gradient of rates across latitude bands in this iteration
  
  for (h in 1:length(sample_pc)){
    
    #Filter rates to one sampling level
    bins_to_rank <- filter(measured_rates, sampling == sample_pc[h]) %>% filter(bin_no != "global")
    
    #Store 100% sampled values as the benchmark for comparison
    if (h == 1) {true_bin_orig <- as.numeric(bins_to_rank$raw_origination_rate);
                  true_bin_ext <- as.numeric(bins_to_rank$raw_extinction_rate)}
    
    #Evaluate whether max and min bins are the same
    #(can tolerate ties in the true values, but not in the calculated rates)
    extremes <- rbind(extremes, c(bins_to_rank[1,1], sample_pc[h], "origination", "raw",
                                    which.min(as.numeric(bins_to_rank$raw_origination_rate)) %in% which(true_bin_orig == min(true_bin_orig)),
                                    which.max(as.numeric(bins_to_rank$raw_origination_rate)) %in% which(true_bin_orig == max(true_bin_orig))))
    extremes <- rbind(extremes, c(bins_to_rank[1,1], sample_pc[h], "extinction", "raw",
                                    which.min(as.numeric(bins_to_rank$raw_extinction_rate)) %in% which(true_bin_ext == min(true_bin_ext)),
                                    which.max(as.numeric(bins_to_rank$raw_extinction_rate)) %in% which(true_bin_ext == max(true_bin_ext))))
    extremes <- rbind(extremes, c(bins_to_rank[1,1], sample_pc[h], "origination", "boundary-crosser",
                                    which.min(as.numeric(bins_to_rank$BC_origination_pc)) %in% which(true_bin_orig == min(true_bin_orig)),
                                    which.max(as.numeric(bins_to_rank$BC_origination_pc)) %in% which(true_bin_orig == max(true_bin_orig))))
    extremes <- rbind(extremes, c(bins_to_rank[1,1], sample_pc[h], "extinction", "boundary-crosser",
                                    which.min(as.numeric(bins_to_rank$BC_extinction_pc)) %in% which(true_bin_ext == min(true_bin_ext)),
                                    which.max(as.numeric(bins_to_rank$BC_extinction_pc)) %in% which(true_bin_ext == max(true_bin_ext))))
    extremes <- rbind(extremes, c(bins_to_rank[1,1], sample_pc[h], "origination", "three-timer",
                                    which.min(as.numeric(bins_to_rank$tt_origination_rate)) %in% which(true_bin_orig == min(true_bin_orig)),
                                    which.max(as.numeric(bins_to_rank$tt_origination_rate)) %in% which(true_bin_orig == max(true_bin_orig))))
    extremes <- rbind(extremes, c(bins_to_rank[1,1], sample_pc[h], "extinction", "three-timer",
                                    which.min(as.numeric(bins_to_rank$tt_extinction_rate)) %in% which(true_bin_ext == min(true_bin_ext)),
                                    which.max(as.numeric(bins_to_rank$tt_extinction_rate)) %in% which(true_bin_ext == max(true_bin_ext))))
    
    #Compare rates from different methods to the true values using Pearson's correlation coefficient
    raw_orig_cor <- cor.test(true_bin_orig, as.numeric(bins_to_rank$raw_origination_rate), method = "pearson")
    raw_ext_cor <- cor.test(true_bin_ext, as.numeric(bins_to_rank$raw_extinction_rate), method = "pearson")
    bc_orig_cor <- cor.test(true_bin_orig, as.numeric(bins_to_rank$BC_origination_pc), method = "pearson")
    bc_ext_cor <- cor.test(true_bin_ext, as.numeric(bins_to_rank$BC_extinction_pc), method = "pearson")
    tt_orig_cor <- cor.test(true_bin_orig, as.numeric(bins_to_rank$tt_origination_rate), method = "pearson")
    tt_ext_cor <- cor.test(true_bin_ext, as.numeric(bins_to_rank$tt_extinction_rate), method = "pearson")
    
    #Compare rates from different methods to the true values using Spearman's rank
    raw_orig_spear <- spearman.test(true_bin_orig, as.numeric(bins_to_rank$raw_origination_rate))
    raw_ext_spear <- spearman.test(true_bin_ext, as.numeric(bins_to_rank$raw_extinction_rate))
    bc_orig_spear <- spearman.test(true_bin_orig, as.numeric(bins_to_rank$BC_origination_pc))
    bc_ext_spear <- spearman.test(true_bin_ext, as.numeric(bins_to_rank$BC_extinction_pc))
    tt_orig_spear <- spearman.test(true_bin_orig, as.numeric(bins_to_rank$tt_origination_rate))
    tt_ext_spear <- spearman.test(true_bin_ext, as.numeric(bins_to_rank$tt_extinction_rate))
    
    #Add test outputs to a list of rate gradients
    gradients <- rbind(gradients, c(bins_to_rank[1,1], sample_pc[h], "origination", "raw",
                                    raw_orig_cor$statistic, raw_orig_cor$p.value, raw_orig_cor$estimate,
                                    raw_orig_spear$statistic, raw_orig_spear$p.value, raw_orig_spear$estimate))
    gradients <- rbind(gradients, c(bins_to_rank[1,1], sample_pc[h], "extinction", "raw",
                                    raw_ext_cor$statistic, raw_ext_cor$p.value, raw_ext_cor$estimate,
                                    raw_ext_spear$statistic, raw_ext_spear$p.value, raw_ext_spear$estimate))
    gradients <- rbind(gradients, c(bins_to_rank[1,1], sample_pc[h], "origination", "boundary-crosser",
                                    bc_orig_cor$statistic, bc_orig_cor$p.value, bc_orig_cor$estimate,
                                    bc_orig_spear$statistic, bc_orig_spear$p.value, bc_orig_spear$estimate))
    gradients <- rbind(gradients, c(bins_to_rank[1,1], sample_pc[h], "extinction", "boundary-crosser",
                                    bc_ext_cor$statistic, bc_ext_cor$p.value, bc_ext_cor$estimate,
                                    bc_ext_spear$statistic, bc_ext_spear$p.value, bc_ext_spear$estimate))
    gradients <- rbind(gradients, c(bins_to_rank[1,1], sample_pc[h], "origination", "three-timer",
                                    tt_orig_cor$statistic, tt_orig_cor$p.value, tt_orig_cor$estimate,
                                    tt_orig_spear$statistic, tt_orig_spear$p.value, tt_orig_spear$estimate))
    gradients <- rbind(gradients, c(bins_to_rank[1,1], sample_pc[h], "extinction", "three-timer",
                                    tt_ext_cor$statistic, tt_ext_cor$p.value, tt_ext_cor$estimate,
                                    tt_ext_spear$statistic, tt_ext_spear$p.value, tt_ext_spear$estimate))
    
    #Filter differences to one sampling level, then split into different methods
    bin_shifts <- filter(measured_diffs, sampling == sample_pc[h]) %>% filter(bin_size != "global")
    bin_shifts_or <- filter(bin_shifts, rate == "origination") %>% filter(method == "raw")
    bin_shifts_er <- filter(bin_shifts, rate == "extinction") %>% filter(method == "raw")
    bin_shifts_obc <- filter(bin_shifts, rate == "origination") %>% filter(method == "boundary-crosser")
    bin_shifts_ebc <- filter(bin_shifts, rate == "extinction") %>% filter(method == "boundary-crosser")
    bin_shifts_ott <- filter(bin_shifts, rate == "origination") %>% filter(method == "three-timer")
    bin_shifts_ett <- filter(bin_shifts, rate == "extinction") %>% filter(method == "three-timer")
    
    #Count the number of bins where the discrepancy between the true and calculated rates is more than 0
    #Currently compares using same comparison as "differences" i.e. to 100% values for same method
    shifts <- rbind(shifts, c(bin_shifts[1,1], bin_shifts[1,2], "origination", "raw", sum(as.numeric(bin_shifts_or$difference > 0))))
    shifts <- rbind(shifts, c(bin_shifts[1,1], bin_shifts[1,2], "extinction", "raw", sum(as.numeric(bin_shifts_er$difference > 0))))
    shifts <- rbind(shifts, c(bin_shifts[1,1], bin_shifts[1,2], "origination", "boundary-crosser", sum(as.numeric(bin_shifts_obc$difference > 0))))
    shifts <- rbind(shifts, c(bin_shifts[1,1], bin_shifts[1,2], "extinction", "boundary-crosser", sum(as.numeric(bin_shifts_ebc$difference > 0))))
    shifts <- rbind(shifts, c(bin_shifts[1,1], bin_shifts[1,2], "origination", "three-timer", sum(as.numeric(bin_shifts_ott$difference > 0))))
    shifts <- rbind(shifts, c(bin_shifts[1,1], bin_shifts[1,2], "extinction", "three-timer", sum(as.numeric(bin_shifts_ett$difference > 0))))
  }
  setTxtProgressBar(pb, x)
}

#Add column names to gradients and shifts tables
colnames(extremes) <- c("iteration_no", "sampling", "rate", "method", "min_match", "max_match")
colnames(gradients) <- c("iteration_no", "sampling", "rate", "method", "t", "t.p_value", "cor", "S", "S.p_value", "rho")
colnames(shifts) <- c("iteration_no", "sampling", "rate", "method", "bins_over")


###Save results###
write.csv(results, "data/Sim_res_overall.csv", row.names = F)
write.csv(differences, "data/Sim_diffs_overall.csv", row.names = F)
write.csv(sampling, "data/Sim_samp_overall.csv", row.names = F)
write.csv(gradients, "data/Sim_grads_overall.csv", row.names = F)
write.csv(shifts, "data/Sim_shifts_overall.csv", row.names = F)


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


#Plot difference between "true" and measured rates at different sampling levels
sampled_g <- filter(differences, sampling != "100") %>% filter(bin_size == "global")
sampled_g$occs <- as.numeric(as.character(sampled_g$occs))
sampled_g$difference <- as.numeric(as.character(sampled_g$difference))

sampled_b <- filter(differences, sampling != "100") %>% filter(bin_size == "lat_band")
sampled_b$occs <- as.numeric(as.character(sampled_b$occs))
sampled_b$difference <- as.numeric(as.character(sampled_b$difference))

ggplot(sampled_g, aes(x = sampling, y = difference, fill = rate)) + geom_hline(aes(yintercept = 0)) +
  geom_boxplot() + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  #scale_y_continuous(limits = c(-1, 1)) +
  theme_classic()
ggplot(sampled_b, aes(x = sampling, y = difference, fill = rate)) + geom_hline(aes(yintercept = 0)) +
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
shifts$bins_over <- as.numeric(as.character(shifts$bins_over))

ggplot(shifts, aes(x = sampling, y = bins_over, fill = rate)) +
  geom_violin(scale = "width") + facet_wrap(~method) + scale_fill_manual(values = c("salmon", "lightblue")) +
  theme_classic()


#Plot sampling level estimated by the 3t method
sampling$estimated <- as.numeric(as.character(sampling$estimated))
sampling$sampled <- as.character(sampling$sampled/100)
ggplot(sampling, aes(x = sampled, y = estimated)) + geom_boxplot(aes(group = sampled)) + theme_classic()
