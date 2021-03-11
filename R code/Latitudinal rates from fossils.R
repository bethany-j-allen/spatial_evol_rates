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


#Generate list of taxa sampled in each latitude band per stage for all clades
brachiopod_list <- list()
bivalve_list <- list()
ammonoid_list <- list()
gastropod_list <- list()

for (i in 1:length(stages)) {
  #Filter stage, then separate the clades
  one_stage <- filter(fossils, stage_bin == stages[i])
  stage_brachs <- filter(one_stage, phylum == "Brachiopoda")
  stage_bivs <- filter(one_stage, class == "Bivalvia")
  stage_amms <- filter(one_stage, class == "Cephalopoda")
  stage_gast <- filter(one_stage, class == "Gastropoda")
  
  #Filter latitude bins, find unique taxa, and put the vector in the list
  for (j in 1:length(labels)){
    
    brachs_one_bin <- filter(stage_brachs, paleolat_code == labels[j])
    if (tax_res == "species"){brachs_in_bin <- sort(unique(brachs_one_bin$accepted_name))
    } else if (tax_res == "genera"){brachs_in_bin <- sort(unique(brachs_one_bin$genus))}
    brachiopod_list[[j]] <- brachs_in_bin
    
    bivs_one_bin <- filter(stage_bivs, paleolat_code == labels[j])
    if (tax_res == "species"){bivs_in_bin <- sort(unique(bivs_one_bin$accepted_name))
    } else if (tax_res == "genera"){bivs_in_bin <- sort(unique(bivs_one_bin$genus))}
    bivalve_list[[j]] <- bivs_in_bin
    
    amms_one_bin <- filter(stage_amms, paleolat_code == labels[j])
    if (tax_res == "species"){amms_in_bin <- sort(unique(amms_one_bin$accepted_name))
    } else if (tax_res == "genera"){amms_in_bin <- sort(unique(amms_one_bin$genus))}
    ammonoid_list[[j]] <- amms_in_bin
    
    gast_one_bin <- filter(stage_gast, paleolat_code == labels[j])
    if (tax_res == "species"){gast_in_bin <- sort(unique(gast_one_bin$accepted_name))
    } else if (tax_res == "genera"){gast_in_bin <- sort(unique(gast_one_bin$genus))}
    gastropod_list[[j]] <- gast_in_bin
  }
  
  #Save list labelled with its stage name
  names(brachiopod_list) <- labels
  names(bivalve_list) <- labels
  names(ammonoid_list) <- labels
  names(gastropod_list) <- labels
  
  eval(parse(text = paste0(stages[i], "_brachiopods <- brachiopod_list")))
  eval(parse(text = paste0(stages[i], "_bivalves <- bivalve_list")))
  eval(parse(text = paste0(stages[i], "_ammonoids <- ammonoid_list")))
  eval(parse(text = paste0(stages[i], "_gastropods <- gastropod_list")))
}


###Choose stage###
#Can't choose first two or last two in vector
stage <- "Wordian"

#Find stage in stages vector, and use this to label the t0, t1 and t2 lists
stage_ID <- match(stage, stages)

#Create vector of clades
clades <- c("brachiopods", "bivalves", "ammonoids", "gastropods")

#Create data frame to store outputs
evol_rates <- data.frame()

for (x in 1:length(clades)){
  
  #Find relevant time bins
  t0_occs <- eval(parse(text = paste0(stages[(stage_ID - 2)], "_", clades[x])))
  t1_occs <- eval(parse(text = paste0(stages[(stage_ID - 1)], "_", clades[x])))
  t2_occs <- eval(parse(text = paste0(stages[(stage_ID)], "_", clades[x])))
  t3_occs <- eval(parse(text = paste0(stages[(stage_ID + 1)], "_", clades[x])))
  t4_occs <- eval(parse(text = paste0(stages[(stage_ID + 2)], "_", clades[x])))

  #Produce global lists of unique occurrences (i.e. richness)
  t0_g_occs <- unique(unlist(t0_occs))
  t1_g_occs <- unique(unlist(t1_occs))
  t2_g_occs <- unique(unlist(t2_occs))
  t3_g_occs <- unique(unlist(t3_occs))
  t4_g_occs <- unique(unlist(t4_occs))

  #Calculate proportions
  #Raw (these counts include singletons)
  global_orig <- length(setdiff(t2_g_occs, t1_g_occs))   #Present in t2 but not in t1
  global_ext <- length(setdiff(t2_g_occs, t3_g_occs))    #Present in t2 but not in t3
  global_orig_p <- global_orig / length(t2_g_occs)
  global_ext_p <- global_ext / length(t2_g_occs)

  #Boundary crosser
  global_originations <- length(setdiff(intersect(t2_g_occs, t3_g_occs), t1_g_occs)) #Present in t2 and t3 but not in t1
  global_extinctions <- length(setdiff(intersect(t1_g_occs, t2_g_occs), t3_g_occs))  #Present in t1 and t2 but not in t3
  global_through <- length(intersect(t1_g_occs, t3_g_occs))
  global_bc_orig <- global_originations/(global_through + global_originations) #BC origination proportion
  global_bc_ext <- global_extinctions/(global_through + global_extinctions)    #BC extinction proportion
  
  #Three-timer
  global_2t_o <- length(intersect(t2_g_occs, t3_g_occs)) #Present in t2 and t3 irrespective of t1
  global_2t_e <- length(intersect(t1_g_occs, t2_g_occs)) #Present in t1 and t2 irrespective of t3
  
  global_3t_1 <- length(intersect(intersect(t0_g_occs, t1_g_occs), t2_g_occs)) #Present in t0-2
  global_3t_2 <- length(intersect(intersect(t1_g_occs, t2_g_occs), t3_g_occs)) #Present in t1-3
  global_3t_3 <- length(intersect(intersect(t2_g_occs, t3_g_occs), t4_g_occs)) #Present in t2-4
  
  global_pt_1 <- length(setdiff(intersect(t0_g_occs, t2_g_occs), t1_g_occs)) #Ghost ranges for t1
  global_pt_2 <- length(setdiff(intersect(t2_g_occs, t4_g_occs), t3_g_occs)) #Ghost ranges for t3
  
  global_sampling_o <- global_3t_1/(global_3t_1 + global_pt_1) #Sampling completeness est. for t1
  global_sampling_e <- global_3t_3/(global_3t_3 + global_pt_2) #Sampling completeness est. for t2
  
  global_3t_orig <- 1 - (global_3t_2/(global_sampling_o*global_2t_o)) #3t origination proportion
  global_3t_ext <- 1 - (global_3t_2/(global_sampling_e*global_2t_e))  #3t extinction proportion
  
  #Add global rates to data frames
  evol_rates <- rbind(evol_rates, c(clades[x], "global", length(t2_g_occs), "origination", "raw", round(global_orig_p, 3)))
  evol_rates <- rbind(evol_rates, c(clades[x], "global", length(t2_g_occs), "extinction", "raw", round(global_ext_p, 3)))
  evol_rates <- rbind(evol_rates, c(clades[x], "global", length(t2_g_occs), "origination", "boundary-crosser", round(global_bc_orig, 3)))
  evol_rates <- rbind(evol_rates, c(clades[x], "global", length(t2_g_occs), "extinction", "boundary-crosser", round(global_bc_ext, 3)))
  evol_rates <- rbind(evol_rates, c(clades[x], "global", length(t2_g_occs), "origination", "three-timer", round(global_3t_orig, 3)))
  evol_rates <- rbind(evol_rates, c(clades[x], "global", length(t2_g_occs), "extinction", "three-timer", round(global_3t_ext, 3)))
  
  
  #Calculate origination and extinction rates for individual latitude bands
  for (y in 1:length(labels)){
    
    #Pull out one latitude band
    focal_bin_t2 <- t2_occs[[y]]

    #Raw (these counts include singletons)
    bin_orig <- length(setdiff(focal_bin_t2, t1_g_occs))   #Present in t2 but not in t1
    bin_ext <- length(setdiff(focal_bin_t2, t3_g_occs))    #Present in t2 but not in t3
    bin_orig_p <- bin_orig / length(focal_bin_t2)
    bin_ext_p <- bin_ext / length(focal_bin_t2)

    #Boundary crosser
    bin_originations <- length(setdiff(intersect(focal_bin_t2, t3_g_occs), t1_g_occs)) #Present in t2 and t3 but not in t1
    bin_extinctions <- length(setdiff(intersect(t1_g_occs, focal_bin_t2), t3_g_occs))  #Present in t1 and t2 but not in t3
    bin_through <- length(intersect(intersect(t1_g_occs, focal_bin_t2), t3_g_occs))    #Present in band for t1, globally for t0 & t2
    bin_bc_orig <- bin_originations/(bin_through + bin_originations)                   #BC origination proportion
    bin_bc_ext <- bin_extinctions/(bin_through + bin_extinctions)                      #BC extinction proportion
    
    #Three-timer
    bin_2t_o <- length(intersect(focal_bin_t2, t3_g_occs)) #Present in t2 and t3 irrespective of t1
    bin_2t_e <- length(intersect(t1_g_occs, focal_bin_t2)) #Present in t1 and t2 irrespective of t3

    bin_3t <- length(intersect(intersect(t1_g_occs, focal_bin_t2), t3_g_occs)) #Present in t1-3
    
    #Uses global sampling probabilities
    
    bin_3t_orig <- 1 - (bin_3t/(global_sampling_o*bin_2t_o)) #3t origination proportion
    bin_3t_ext <- 1 - (bin_3t/(global_sampling_e*bin_2t_e))  #3t extinction proportion
    
    #Save in a vector
    evol_rates <- rbind(evol_rates, c(clades[x], labels[y], length(focal_bin_t2), "origination", "raw", round(bin_orig_p, 3)))
    evol_rates <- rbind(evol_rates, c(clades[x], labels[y], length(focal_bin_t2), "extinction", "raw", round(bin_ext_p, 3)))
    evol_rates <- rbind(evol_rates, c(clades[x], labels[y], length(focal_bin_t2), "origination", "boundary-crosser", round(bin_bc_orig, 3)))
    evol_rates <- rbind(evol_rates, c(clades[x], labels[y], length(focal_bin_t2), "extinction", "boundary-crosser", round(bin_bc_ext, 3)))
    evol_rates <- rbind(evol_rates, c(clades[x], labels[y], length(focal_bin_t2), "origination", "three-timer", round(bin_3t_orig, 3)))
    evol_rates <- rbind(evol_rates, c(clades[x], labels[y], length(focal_bin_t2), "extinction", "three-timer", round(bin_3t_ext, 3)))
    }
}
  
#Label columns in rates data frames
colnames(evol_rates) <- c("clade", "bin", "richness", "rate", "method", "value")

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
