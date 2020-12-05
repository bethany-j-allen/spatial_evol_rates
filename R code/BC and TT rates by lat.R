#Bethany Allen   19th November 2020 (updated 5th December 2020)
#Code to run boundary crosser and gap-filler methods of calculating origination and extinction rates

#setwd("#####")

#Load packages
library(tidyverse)

#Create a vector giving the chronological order of stages
stages <- c("Roadian", "Wordian", "Capitanian", "Wuchiapingian", "Changhsingian", "Induan", "Olenekian",
            "Anisian", "Ladinian", "Carnian")

#Designate latitude bins
bins <- seq(from = -90, to = 90, by = 30)
labels <- seq(from = -75, to = 75, by = 30)

#Read in dataset
fossils <- read_csv("data/PT_brach_biv_clean.csv")
glimpse(fossils)

###Designate the taxonomic resolution ("species" or "genera")###
tax_res <- "species"

#Apply filters
if (tax_res == "species"){fossils <- filter(fossils, accepted_rank == "species")
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



rates <- data.frame()

for (i in 3:(length(stages) - 2)) {
  #Find raw rates
  global_orig <- length(setdiff(taxon_lists[[i]], taxon_lists[[i-1]]))
  global_ext <- length(setdiff(taxon_lists[[i]], taxon_lists[[i+1]]))
  global_orig_p <- global_orig / length(taxon_lists[[i]])
  global_ext_p <- global_ext / length(taxon_lists[[i]])
  
  #Find forward boundary crossers (not previously known)
  bc_forward <- length(setdiff(intersect(taxon_lists[[i]], taxon_lists[[i+1]]), taxon_lists[[i-1]]))
  
  #Find backward boundary crossers (not subsequently known)
  bc_backward <- length(setdiff(intersect(taxon_lists[[i-1]], taxon_lists[[i]]), taxon_lists[[i+1]]))
  
  #Find range-throughs
  bc_range <- length(intersect(taxon_lists[[(i-1)]], taxon_lists[[i+1]]))
  
  #Calculate per-capita BC origination rate
  bc_orig <- log((bc_range + bc_forward)/bc_range)
  
  #Calculate per-capita BC extinction rate
  bc_ext <- log((bc_range + bc_backward)/bc_range)

  
  #Find upper two timers, those that continue through the upper boundary of the bin
  orig_twotimers <- length(intersect(taxon_lists[[i]], taxon_lists[[(i+1)]]))
  
  #Find lower two timers, those that are known from the previous bin and continue into the focal one
  ext_twotimers <- length(intersect(taxon_lists[[(i-1)]], taxon_lists[[i]]))
  
  #Find three timers, those known from the focal bin, the one before and the one after
  threetimers <- length(intersect((intersect(taxon_lists[[(i-1)]], taxon_lists[[i]])), taxon_lists[[(i+1)]]))
  
  #Find part timers, those known from the bins before and after, but not the focal one
  parttimers <- length(setdiff((intersect(taxon_lists[[(i-1)]], taxon_lists[[(i+1)]])), taxon_lists[[i]]))
  
  #Find backward gap fillers, known from second previous and forward bins but not one previous
  orig_gapfillers <- length(setdiff((intersect(taxon_lists[[(i-2)]], taxon_lists[[(i+1)]])), taxon_lists[[i-1]]))
  
  #Find forward gap fillers, known from previous and second forward bins but not one forward
  ext_gapfillers <- length(setdiff((intersect(taxon_lists[[(i-1)]], taxon_lists[[(i+2)]])), taxon_lists[[i+1]]))
  
  #Find origination rate (lambda) and origination proportion
  orig_rate <- log((orig_twotimers + parttimers)/(threetimers + parttimers + orig_gapfillers))
  orig_prop <- 1 - ((threetimers + parttimers + orig_gapfillers)/(orig_twotimers + parttimers))
  
  #Find extinction rate (miu) and extinction proportion
  ext_rate <- log((ext_twotimers + parttimers)/(threetimers + parttimers + ext_gapfillers))
  ext_prop <- 1 - ((threetimers + parttimers + ext_gapfillers)/(ext_twotimers + parttimers))
  
  rates <- rbind(rates, c(global_orig_p, global_ext_p, bc_orig, bc_ext, orig_rate, orig_prop, ext_rate, ext_prop))
}

colnames(rates) <- c("raw_orig_prop", "raw_ext_prop", "bc_orig_rate", "bc_ext_rate", "gf_orig_rate", "gf_orig_prop", "gf_ext_rate", "gf_ext_prop")
rownames(rates) <- stages[3:(length(stages) - 2)]

View(rates)
