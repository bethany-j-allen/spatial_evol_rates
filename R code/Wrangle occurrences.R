#Bethany Allen   30th June 2020 (updated 29th July 2020, 5th Dec 2020)
#Code to clean the PBDB marine invertebrate data

#setwd("#####")

#Load packages
library(tidyverse)

#Create a vector giving the chronological order of stages
stages <- c("Artinskian", "Kungurian", "Roadian", "Wordian", "Capitanian", "Wuchiapingian",
            "Changhsingian", "Induan", "Olenekian", "Anisian", "Ladinian", "Carnian", "Norian")

#Create a vector giving the chronological order of substages
substages <- c("Artinskian", "Kungurian", "Roadian", "Wordian", "Capitanian", "Wuchiapingian",
               "Changhsingian", "Griesbachian", "Dienerian", "Smithian", "Spathian", "Aegean",
               "Bithynian", "Pelsonian", "Illyrian", "Fassanian", "Longobardian", "Julian",
               "Tuvalian", "Lacian", "Alaunian", "Sevatian")


#Read in dataset
fossils <- read_csv("data/PT_marine_inverts.csv")
glimpse(fossils)

#Add filters to remove lacustrine occurrences
#fluvial  <- c("fluvial-lacustrine indet.", "fluvial indet.", "lacustrine - large",
#              "lacustrine delta front", "lacustrine indet.", "pond", "terrestrial indet.")
#fossils <- filter(fossils, !environment %in% fluvial)

#Add filters to remove uncertain IDs
#fossils <- fossils %>% filter(!str_detect(identified_name, " cf"))
#fossils <- fossils %>% filter(!str_detect(identified_name, " aff"))
#fossils <- fossils %>% filter(!str_detect(identified_name, '"'))
#fossils <- fossils %>% filter(!str_detect(identified_name, " \\?"))
#fossils <- fossils %>% filter(!str_detect(identified_name, "ex gr."))

###Bin occurrences by stage and substage###
#Create columns for stage and substage designation
fossils$stage_bin <- NA
fossils$substage_bin <- NA

#For each occurrence
for (i in 1:nrow(fossils)){
  #If occurrence is dated to a single stage, allocate it to that bin
  if (fossils$early_interval[i] %in% stages & is.na(fossils$late_interval[i])){
    fossils$stage_bin[i] <- fossils$early_interval[i]}
  #If occurrence is dated to a single substage, allocate it to that bin
  if (fossils$early_interval[i] %in% substages & is.na(fossils$late_interval[i])){
    fossils$substage_bin[i] <- fossils$early_interval[i]}
  #If occurrence is dated to Griesbachian/Dienerian or both, it is Induan
  if (fossils$early_interval[i] %in% substages[8:9] & is.na(fossils$late_interval[i])){
    fossils$stage_bin[i] <- "Induan"}
  if (fossils$early_interval[i] == substages[8] & !is.na(fossils$late_interval[i])){
    if(fossils$late_interval[i] == substages[9]){fossils$stage_bin[i] <- "Induan"}}
  #If occurrence is dated to Smithian/Spathian or both, it is Olenekian
  if (fossils$early_interval[i] %in% substages[10:11] & is.na(fossils$late_interval[i])){
    fossils$stage_bin[i] <- "Olenekian"}
  if (fossils$early_interval[i] == substages[10] & !is.na(fossils$late_interval[i])){
    if(fossils$late_interval[i] == substages[11]){fossils$stage_bin[i] <- "Olenekian"}}
  #If occurrence is dated to Aegean/Bithynian/Pelsonian/Illyrian or a combination, it is Anisian
  if (fossils$early_interval[i] %in% substages[12:15] & is.na(fossils$late_interval[i])){
    fossils$stage_bin[i] <- "Anisian"}
  if (fossils$early_interval[i] %in% substages[12:15] & !is.na(fossils$late_interval[i])){
    if(fossils$late_interval[i] %in% substages[13:15]){fossils$stage_bin[i] <- "Anisian"}}
  #If occurrence is dated to Fassanian/Longobardian or both, it is Ladinian
  if (fossils$early_interval[i] %in% substages[16:17] & is.na(fossils$late_interval[i])){
    fossils$stage_bin[i] <- "Ladinian"}
  if (fossils$early_interval[i] == substages[16] & !is.na(fossils$late_interval[i])){
    if(fossils$late_interval[i] == substages[17]){fossils$stage_bin[i] <- "Ladinian"}}
  #If occurrence is dated to Julian/Tuvalian or both, it is Carnian
  if (fossils$early_interval[i] %in% substages[18:19] & is.na(fossils$late_interval[i])){
    fossils$stage_bin[i] <- "Carnian"}
  if (fossils$early_interval[i] == substages[18] & !is.na(fossils$late_interval[i])){
    if(fossils$late_interval[i] == substages[19]){fossils$stage_bin[i] <- "Carnian"}}
  #If occurrence is dated to Lacian/Alaunian/Sevatian or a combination, it is Norian
  if (fossils$early_interval[i] %in% substages[20:22] & is.na(fossils$late_interval[i])){
    fossils$stage_bin[i] <- "Norian"}
  if (fossils$early_interval[i] %in% substages[20:22] & !is.na(fossils$late_interval[i])){
    if(fossils$late_interval[i] %in% substages[21:22]){fossils$stage_bin[i] <- "Norian"}}
}

#Remove occurrences undated at stage resolution
fossils <- filter(fossils, !is.na(stage_bin))


###Retain uncatalogued species###
#If an occurrence is to species level but the species hasn't been entered into the database, convert
# its accepted name/rank to the species rather than the genus
for (j in 1:nrow(fossils)){
  if(!is.na(fossils$difference[j]))
    (if (fossils$difference[j] == "species not entered"){
      fossils$accepted_name[j] <- fossils$identified_name[j]
      fossils$accepted_rank[j] <- "species"})
}


###Remove synonymy repeats (combinations of the same collection no. AND accepted name)###
fossils <- distinct(fossils, accepted_name, collection_no, .keep_all = T)


###Rotate palaeo-occurrences if desired - can just use "paleolat" column below for PBDB rotations###


###Allocate occurrences to a latitude bin###
#Create latitude bins
bins <- seq(from = -90, to = 90, by = 30)
labels <- seq(from = -75, to = 75, by = 30)

#Add latitude band as an extra variable to the original dataset
fossils <- mutate(fossils, paleolat_code = cut(rot_lat, breaks = bins, labels = labels))
fossils <- filter(fossils, !is.na(paleolat_code))


###Save cleaned dataset###
write.csv(fossils, "data/PT_marine_inverts_clean.csv", row.names = F)
