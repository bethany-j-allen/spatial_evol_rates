#Bethany Allen   19th November 2020 (updated 5th December 2020, 10th March 2021)
#Code to calculate raw, boundary-crosser and three-timer evolutionary rates by latitude bin

#setwd("#####")

#Load packages
library(tidyverse)

#Create a vector giving the chronological order of stages
stages <- c("Artinskian", "Kungurian", "Roadian", "Wordian", "Capitanian", "Wuchiapingian",
            "Changhsingian", "Induan", "Olenekian", "Anisian", "Ladinian", "Carnian", "Norian")

#Designate latitude bins
bins <- seq(from = -90, to = 90, by = 30)
labels <- seq(from = -75, to = 75, by = 30)

#Read in dataset
fossils <- read_csv("data/PT_marine_inverts_clean.csv")
glimpse(fossils)

###Designate the taxonomic resolution ("species" or "genera")###
tax_res <- "species"

#Apply filters
if (tax_res == "species"){fossils <- fossils %>% filter(accepted_rank == "species") %>%
  filter(!str_detect(identified_name, " cf")) %>% filter(!str_detect(identified_name, " aff")) %>%
  filter(!str_detect(identified_name, '"')) %>% filter(!str_detect(identified_name, " \\?")) %>%
  filter(!str_detect(identified_name, "ex gr."))
} else if (tax_res == "genera"){fossils <- filter(fossils, !is.na(genus))}


#Generate list of taxa sampled in each latitude band per stage for both clades
brachs_list <- list()
bivs_list <- list()

for (i in 1:length(stages)) {
  #Filter stage, then separate the clades
  one_stage <- filter(fossils, stage_bin == stages[i])
  stage_brachs <- filter(one_stage, phylum == "Brachiopoda")
  stage_bivs <- filter(one_stage, phylum == "Mollusca")
  #Filter latitude bins, find unique taxa, and put the vector in the list
  for (j in 1:length(labels)){
    brachs_one_bin <- filter(stage_brachs, paleolat_code == labels[j])
    if (tax_res == "species"){brachs_in_bin <- sort(unique(brachs_one_bin$accepted_name))
    } else if (tax_res == "genera"){brachs_in_bin <- sort(unique(brachs_one_bin$genus))}
    brachs_list[[j]] <- brachs_in_bin
    bivs_one_bin <- filter(stage_bivs, paleolat_code == labels[j])
    if (tax_res == "species"){bivs_in_bin <- sort(unique(bivs_one_bin$accepted_name))
    } else if (tax_res == "genera"){bivs_in_bin <- sort(unique(bivs_one_bin$genus))}
    bivs_list[[j]] <- bivs_in_bin
  }
  #Save list labelled with its stage name
  names(brachs_list) <- labels
  names(bivs_list) <- labels
  eval(parse(text = paste0(stages[i], "_brachs <- brachs_list")))
  eval(parse(text = paste0(stages[i], "_bivs <- bivs_list")))
}


###Choose stage###
#Can't choose first or last stage in the vector
stage <- "Wordian"

#Find stage in stages vector, and use this to label the t0, t1 and t2 lists
stage_ID <- match(stage, stages)

t0_brachs <- eval(parse(text = paste0(stages[(stage_ID - 1)], "_brachs")))
t1_brachs <- eval(parse(text = paste0(stages[(stage_ID)], "_brachs")))
t2_brachs <- eval(parse(text = paste0(stages[(stage_ID + 1)], "_brachs")))

t0_bivs <- eval(parse(text = paste0(stages[(stage_ID - 1)], "_bivs")))
t1_bivs <- eval(parse(text = paste0(stages[(stage_ID)], "_bivs")))
t2_bivs <- eval(parse(text = paste0(stages[(stage_ID + 1)], "_bivs")))
  
evol_rates <- data.frame()

#Produce global lists of unique occurrences (i.e. richness)
t0_g_brachs <- unique(unlist(t0_brachs))
t1_g_brachs <- unique(unlist(t1_brachs))
t2_g_brachs <- unique(unlist(t2_brachs))

t0_g_bivs <- unique(unlist(t0_bivs))
t1_g_bivs <- unique(unlist(t1_bivs))
t2_g_bivs <- unique(unlist(t2_bivs))
  
#Raw (these counts include singletons)
br_g_orig <- length(setdiff(t1_g_brachs, t0_g_brachs))   #Present in t1 but not in t0
br_g_ext <- length(setdiff(t1_g_brachs, t2_g_brachs))    #Present in t1 but not in t2
br_g_orig_p <- br_g_orig / length(t1_g_brachs)
br_g_ext_p <- br_g_ext / length(t1_g_brachs)

bi_g_orig <- length(setdiff(t1_g_bivs, t0_g_bivs))   #Present in t1 but not in t0
bi_g_ext <- length(setdiff(t1_g_bivs, t2_g_bivs))    #Present in t1 but not in t2
bi_g_orig_p <- bi_g_orig / length(t1_g_bivs)
bi_g_ext_p <- bi_g_ext / length(t1_g_bivs)
  
#Boundary crosser
br_g_originations <- length(setdiff(intersect(t1_g_brachs, t2_g_brachs), t0_g_brachs)) #Present in t1 and t2 but not in t0
br_g_extinctions <- length(setdiff(intersect(t0_g_brachs, t1_g_brachs), t2_g_brachs))  #Present in t0 and t1 but not in t2
br_g_through <- length(intersect(t0_g_brachs, t2_g_brachs))                         
br_g_bc_orig <- log((br_g_through + br_g_originations)/br_g_through)       #Per-capita BC origination rate
br_g_bc_ext <- log((br_g_through + br_g_extinctions)/br_g_through)         #Per-capita BC extinction rate

bi_g_originations <- length(setdiff(intersect(t1_g_bivs, t2_g_bivs), t0_g_bivs)) #Present in t1 and t2 but not in t0
bi_g_extinctions <- length(setdiff(intersect(t0_g_bivs, t1_g_bivs), t2_g_bivs))  #Present in t0 and t1 but not in t2
bi_g_through <- length(intersect(t0_g_bivs, t2_g_bivs))                         
bi_g_bc_orig <- log((bi_g_through + bi_g_originations)/bi_g_through)       #Per-capita BC origination rate
bi_g_bc_ext <- log((bi_g_through + bi_g_extinctions)/bi_g_through)         #Per-capita BC extinction rate
  
#Three-timer
  #As sampling is being fixed through time, the sampling rate here is calculated from t1
    #global_2t_1 <- length(intersect(t0_global, t1_global)) #Present in t0 and t1 irrespective of t2
    #global_2t_2 <- length(intersect(t1_global, t2_global)) #Present in t1 and t2 irrespective of t0
    #global_3t <- length(intersect(intersect(t0_global, t1_global), t2_global)) #Present in all t
    #global_pt <- length(setdiff(intersect(t0_global, t2_global), t1_global)) #t1 ghost ranges
    #global_t1_sampling <- global_3t/(global_3t + global_pt)
    #global_3t_orig <- log(global_2t_2/global_3t) + log(global_t1_sampling) #3t origination rate
    #global_3t_ext <- log(global_2t_1/global_3t) + log(global_t1_sampling)  #3t extinction rate
    #global_3t_orig_diff <- global_orig_p - global_3t_orig  #Difference between raw and 3t origination
    #global_3t_ext_diff <- global_ext_p - global_3t_ext     #Difference between raw and 3t extinction
  
#Add global rates to data frames
evol_rates <- rbind(evol_rates, c("Brachiopods", "global", length(t1_g_brachs), br_g_orig, "origination", "raw", round(br_g_orig_p, 3)))
evol_rates <- rbind(evol_rates, c("Brachiopods", "global", length(t1_g_brachs), br_g_ext, "extinction", "raw", round(br_g_ext_p, 3)))
evol_rates <- rbind(evol_rates, c("Brachiopods", "global", length(t1_g_brachs), br_g_orig, "origination", "boundary-crosser", round(br_g_bc_orig, 3)))
evol_rates <- rbind(evol_rates, c("Brachiopods", "global", length(t1_g_brachs), br_g_ext, "extinction", "boundary-crosser", round(br_g_bc_ext, 3)))
evol_rates <- rbind(evol_rates, c("Bivalves", "global", length(t1_g_bivs), bi_g_orig, "origination", "raw", round(bi_g_orig_p, 3)))
evol_rates <- rbind(evol_rates, c("Bivalves", "global", length(t1_g_bivs), bi_g_ext, "extinction", "raw", round(bi_g_ext_p, 3)))
evol_rates <- rbind(evol_rates, c("Bivalves", "global", length(t1_g_bivs), bi_g_orig, "origination", "boundary-crosser", round(bi_g_bc_orig, 3)))
evol_rates <- rbind(evol_rates, c("Bivalves", "global", length(t1_g_bivs), bi_g_ext, "extinction", "boundary-crosser", round(bi_g_bc_ext, 3)))

#Calculate origination and extinction rates for individual latitude bands
  
for (e in 1:length(labels)){
    #Pull out one latitude band
    focal_bin_t1_br <- t1_brachs[[e]]
    focal_bin_t1_bi <- t1_bivs[[e]]
    
    #Raw (these counts include singletons)
    br_b_orig <- length(setdiff(focal_bin_t1_br, t0_g_brachs))   #Present in t1 band but nowhere in t0
    br_b_ext <- length(setdiff(focal_bin_t1_br, t2_g_brachs))    #Present in t1 band but nowhere in t2
    br_b_orig_prop <- br_b_orig / length(focal_bin_t1_br)
    br_b_ext_prop <- br_b_ext / length(focal_bin_t1_br)
    
    bi_b_orig <- length(setdiff(focal_bin_t1_bi, t0_g_bivs))   #Present in t1 band but nowhere in t0
    bi_b_ext <- length(setdiff(focal_bin_t1_bi, t2_g_bivs))    #Present in t1 band but nowhere in t2
    bi_b_orig_prop <- bi_b_orig / length(focal_bin_t1_bi)
    bi_b_ext_prop <- bi_b_ext / length(focal_bin_t1_bi)
    
    #Boundary crosser
    br_b_originations <- length(setdiff(intersect(focal_bin_t1_br, t2_g_brachs), t0_g_brachs))  #Present in t1 and t2 but not in t0
    br_b_extinctions <- length(setdiff(intersect(t0_g_brachs, focal_bin_t1_br), t2_g_brachs))   #Present in t0 and t1 but not in t2
    br_b_through <- length(intersect(intersect(t0_g_brachs, focal_bin_t1_br), t2_g_brachs))     #Present in band for t1, globally for t0 & t2
    br_b_bc_orig <- log((br_b_through + br_b_originations)/br_b_through)                        #Per-capita BC origination rate
    br_b_bc_ext <- log((br_b_through + br_b_extinctions)/br_b_through)                          #Per-capita BC extinction rate
    
    bi_b_originations <- length(setdiff(intersect(focal_bin_t1_bi, t2_g_bivs), t0_g_bivs))  #Present in t1 and t2 but not in t0
    bi_b_extinctions <- length(setdiff(intersect(t0_g_bivs, focal_bin_t1_bi), t2_g_bivs))   #Present in t0 and t1 but not in t2
    bi_b_through <- length(intersect(intersect(t0_g_bivs, focal_bin_t1_bi), t2_g_bivs))     #Present in band for t1, globally for t0 & t2
    bi_b_bc_orig <- log((bi_b_through + bi_b_originations)/bi_b_through)                        #Per-capita BC origination rate
    bi_b_bc_ext <- log((bi_b_through + bi_b_extinctions)/bi_b_through)                          #Per-capita BC extinction rate
    
    #Three-timer - uses global sampling probability
      #bin_2t_1 <- length(intersect(t0_global, focal_bin_t1)) #Present in t0 and t1 irrespective of t2
      #bin_2t_2 <- length(intersect(focal_bin_t1, t2_global)) #Present in t1 and t2 irrespective of t0
      #bin_3t <- length(intersect(intersect(t0_global, focal_bin_t1), t2_global)) #Present in all t
      #bin_3t_orig <- log(bin_2t_2/bin_3t) + log(global_t1_sampling)      #3t origination rate
      #bin_3t_ext <- log(bin_2t_1/bin_3t) + log(global_t1_sampling)       #3t extinction rate
      #bin_3t_orig_diff <- bin_orig_prop - bin_3t_orig  #Difference between raw and 3t origination
      #bin_3t_ext_diff <- bin_ext_prop - bin_3t_ext     #Difference between raw and 3t extinction
    
    #Save in a vector
    evol_rates <- rbind(evol_rates, c("Brachiopods", labels[e], length(focal_bin_t1_br), br_b_orig, "origination", "raw", round(br_b_orig_prop, 3)))
    evol_rates <- rbind(evol_rates, c("Brachiopods", labels[e], length(focal_bin_t1_br), br_b_ext, "extinction", "raw", round(br_b_ext_prop, 3)))
    evol_rates <- rbind(evol_rates, c("Brachiopods", labels[e], length(focal_bin_t1_br), br_b_orig, "origination", "boundary-crosser", round(br_b_bc_orig, 3)))
    evol_rates <- rbind(evol_rates, c("Brachiopods", labels[e], length(focal_bin_t1_br), br_b_ext, "extinction", "boundary-crosser", round(br_b_bc_ext, 3)))
    evol_rates <- rbind(evol_rates, c("Bivalves", labels[e], length(focal_bin_t1_bi), bi_b_orig, "origination", "raw", round(bi_b_orig_prop, 3)))
    evol_rates <- rbind(evol_rates, c("Bivalves", labels[e], length(focal_bin_t1_bi), bi_b_ext, "extinction", "raw", round(bi_b_ext_prop, 3)))
    evol_rates <- rbind(evol_rates, c("Bivalves", labels[e], length(focal_bin_t1_bi), bi_b_orig, "origination", "boundary-crosser", round(bi_b_bc_orig, 3)))
    evol_rates <- rbind(evol_rates, c("Bivalves", labels[e], length(focal_bin_t1_bi), bi_b_ext, "extinction", "boundary-crosser", round(bi_b_bc_ext, 3)))
}
  
#Label columns in rates data frames
colnames(evol_rates) <- c("clade", "bin", "bin_size", "raw_n", "rate", "method", "value")

#Amend table for plot
evol_rates_b <- filter(evol_rates, bin != "global")
evol_rates_b <- filter(evol_rates_b, method == "raw")
evol_rates_b$bin <- as.numeric(as.character(evol_rates_b$bin))
evol_rates_b$value <- as.numeric(as.character(evol_rates_b$value))

#Plot -> geom_line(size = 1.5m aes(linetype = method))
ggplot(evol_rates_b, aes(x = bin, y = value, colour = rate)) +
  geom_line(size = 1.5) + geom_point(size = 2) + facet_wrap(~clade) +
  labs(x = "Palaeolatitude", y = "Rate") + 
  coord_flip() + scale_x_continuous(limits = c(-80, 80)) + scale_y_continuous(limits = c(0, 1)) +
  #scale_colour_manual(values = c("blue", "limegreen")) +
  geom_vline(aes(xintercept = 0), colour = "black", size = 0.7) +
  theme_classic()
