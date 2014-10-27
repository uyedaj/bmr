#this is a quick script to clean the names in file "Whiteetal_amph.csv"

library(ape)
library(dplyr)

tree <- read.tree("squa.tre")

#set strings as factors to false to avoid problems later

traits <- read.csv("Whiteetal_squa.csv", stringsAsFactors=FALSE)

names(traits) <- c("group","species","mass","temp","MR","mean.mass","q10smr","ref","obs")

#get rid of unnecessary columns

traits <- select(traits, species, mass, temp, MR, mean.mass, q10smr)

#what names are in the trait database that aren't in the tree?

examine <- setdiff(traits$species, tree$tip.label)

#ok, it's down to a reasonable number. send the wrong names out to excel to check

write.csv(examine, "look_here.csv", row.names=FALSE)

#manually cleaned. some species aren't in the phylogeny. note some useful info in the file
#we use here, we can add some of these species in later following these notes if we want

use <- read.csv("White_squa_clean1.csv", stringsAsFactors=FALSE)

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

#ok everything matches. need to cut the ones that are labeled "CUT"

traits <- traits[traits$species != "CUT", ]

write.csv(traits, "White_Squa.csv", row.names=FALSE)





