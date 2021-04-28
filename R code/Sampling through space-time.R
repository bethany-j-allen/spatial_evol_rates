#Bethany Allen   28th April 2021
#Code to plot sampling through space and time

#setwd("#####")

#Load packages
library(tidyverse)

#Create a vector giving the chronological order of stages
stages <- c("Roadian", "Wordian", "Capitanian", "Wuchiapingian", "Changhsingian", "Induan", "Olenekian",
            "Anisian", "Ladinian")

#Define latitude bins
bins <- seq(from = -90, to = 90, by = 30); labels <- seq(from = -75, to = 75, by = 30)

#Read in data
fossils <- read_csv("data/PT_marine_inverts_clean.csv")

#Indicate phyla
clade <- c("Brachiopoda", "Bivalvia", "Cephalopoda", "Gastropoda")

#Create empty data frame
occ_counts <- data.frame()


for (p in 1:length(clade)){
  if(clade[p] == "Brachiopoda"){one_clade <- filter(fossils, phylum == clade[p])} else {
    one_clade <- filter(fossils, class == clade[p])}
  
  for (k in 1:length(stages)){
  
    #Occurrence counts
    occs_stage <- one_clade %>% filter(stage_bin == stages[k]) %>% count(paleolat_code)
    for(m in 1:length(labels)) {
      if((labels[m] %in% occs_stage$paleolat_code) == FALSE) occs_stage <- rbind(occs_stage, c(labels[m], 0))
    }
    occs_stage$stage <- stages[k]
    occs_stage$clade <- clade[p]
    occ_counts <- rbind(occ_counts, occs_stage)
  }
}

group_by(occ_counts, clade) %>% count(n >= 200)

occ_counts2 <- pivot_wider(occ_counts, names_from = stage, values_from = n)

#Plot occurrences
occ_counts_stage <- filter(occ_counts, stage == "Roadian")

ggplot(occ_counts_stage, aes(x = paleolat_code, y = n, group = clade, colour = clade)) +
  geom_line(size = 2) + geom_point(size = 2) + labs(x = "Palaeolatitude", y = "Occurrence counts") +
  coord_flip() + scale_x_continuous(limits = c(-90, 90)) +
  geom_vline(aes(xintercept = 0), colour = "black") + 
  scale_colour_manual(values = c("grey", "black")) + theme_classic()
