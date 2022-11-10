### Bethany Allen 03.11.2022 ###
### BES Palaeo SIG R workshop ###
### Example: evolutionary simulations ###

#Imagine we want to simulate species occurrences across six latitudinal bins

#Determine how many species we want
species <- base::sample(x = 5, size = 1)
species

#How would this work for one bin?
#Determine how many occurrences
occ_count <- base::sample(x = 10, size = 1)
#Sample species from our list of species identities
occurrences <- base::sample(x = species, size = occ_count, replace = TRUE)
occurrences

#Now we can do this for all six of our latitude bins
#Set the number of bins
nbins <- 6
#Determine how many occurrences, sampling a number for each of the bins
occ_count <- base::sample(x = 10, size = nbins, replace = TRUE)
occ_count
#We can now write a for loop that samples species for each of our latitude bins
#Create an empty list to store our occurrences in
world <- list()
#Iterate through each bin, drawing species from the species pool
for (i in 1:nbins){
  occurrences <- base::sample(x = species, size = occ_count[i], replace = TRUE)
  world[[i]] <- occurrences
}
#This is the starting set-up of our simulation!
world

#We can determine simple summary statistics using length() and unique()
#How many occurrences are there in the 1st bin?
length(world[[1]])
#What is the species diversity in the 1st bin?
length(unique(world[[1]]))
#How many occurrences are there across all six bins, or the whole world?
length(unlist(world))
#What is our global species diversity?
length(unique(unlist(world)))

#What would our world look like after sampling (in the fossil record)?
#Take the occurrences in one bin
occurrences <- world[[1]]
#Determine how many occurrences to fossilise, between 1 and all of them
sample_size <- base::sample(x = length(occurrences), size = 1)
#Sample the occurrences to fossilise
fossils <- base::sample(x = occurrences, size = sample_size, replace = FALSE)
#Compare the full and sampled vectors
occurrences
fossils

#What if we wanted the probability of selecting our samples to be non-uniform?
species <- base::sample(x = 5, size = 1, prob = c(1, 1, 1, 1, 2))
species
#There are also functions in the stats package that allow you to sample from
# various distributions (e.g. rnorm(), rlnorm(), see ?Distributions)

#We now want to subject our initial "world" to extinction and speciation
#We will do this bin-by-bin
#First, describe the range of possible survival percentages (with non-survivors
# becoming extinct)
survival_val <- seq(from = 0, to = 50, by = 0.01)
survival_val
#Sample survival percentages from this vector for each latitude bin
survival_prop <- base::sample(x = survival_val, size = nbins,
                              replace = TRUE) / 100
#Convert the survival decimal to a whole number of occurrences
survival_counts <- round((lengths(world) * survival_prop), 0)
survival_counts

#Now, speciation: we will start by increasing the size of our species pool
species <- species + base::sample(x = 5, size = 1)
species
#Determine how many new occurrences each latitude bin will gain
new_occs <- base::sample(x = 10, size = nbins, replace = TRUE)
new_occs

#Finally, we will write a loop which actually carries out the extinction and
# origination on our initial latitude bins, and saves the outcomes
#Create an empty list to store our new occurrences in
new_world <- list()
#Iterate through each bin, selecting occurrences to survive and adding new ones
for (j in 1:nbins){
  one_bin <- world[[j]]
  one_bin <- base::sample(x = one_bin, size = survival_counts[j],
                          replace = FALSE)
  new_occ_id <- base::sample(x = species, size = new_occs[j], replace = TRUE)
  one_bin <- append(one_bin, new_occ_id)
  new_world[[j]] <- one_bin
}
#And that's it! You can compare the first time slice to the second:
world
new_world
