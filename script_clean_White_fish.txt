#this is a quick script to clean the names in file "Whiteetal_fish.csv"

library(ape)
library(dplyr)

tree <- read.tree("fish.tre")

#set strings as factors to false to avoid problems later

traits <- read.csv("Whiteetal_fish.csv", stringsAsFactors=FALSE)

#get rid of unnecessary columns

traits <- select(traits, species, mass, temp, MR, mean.mass, q10smr)

#what names are in the trait database that aren't in the tree?

examine <- setdiff(traits$species, tree$tip.label)

#ok, it's down to a reasonable number. send the wrong names out to excel to check

write.csv(examine, "look_here.csv", row.names=FALSE)

#manually cleaned. some species aren't in the phylogeny. THIS INCLUDES ELASMOBRANCHS,
#HAGFISH, LUNGFISH, LAMPREYS. Cutting these. Would need to add these into the phylogeny
#manually if we want them! 

use <- read.csv("White_Fish_clean1.csv", stringsAsFactors=FALSE)

#write a stupid for loop because i don't have a comprehensive file to run a merge command
#over, and because there are multiple instances in the input file, so can't just replace
#easily (could combine with all unique names and do a better job of this)

for(i in 1:length(traits$species))
{
	if(traits$species[i] %in% use$original)
	{
		original <- traits$species[i]
		replacement <- use$cleaned[use$original == original]
		traits$species[i] <- replacement
	}
	else
	{
		traits$species[i] <- traits$species[i]
	}
}

#now check to see if there are any additional problems

examine <- setdiff(traits$species, tree$tip.label)

#ok everything matches. need to cut the ones that are labeled "CUT" and CUT_ETC

traits <- traits[traits$species != "CUT", ]
traits <- traits[traits$species != "CUT_ELASMO", ]
traits <- traits[traits$species != "CUT_LAMPREY", ]
traits <- traits[traits$species != "CUT_HAGFISH", ]
traits <- traits[traits$species != "CUT_LUNGFISH", ]

write.csv(traits, "White_Fishes.csv", row.names=FALSE)





