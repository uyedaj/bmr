#this script combines our cleaned and prepped White et al BMR~mass databases into a
#single CSV

library(ape)
library(dplyr)

#load the datasets

fish <- read.csv("White_Fishes.csv")
amph <- read.csv("White_Amph.csv")
bird <- read.csv("Whiteetal_bird.csv")
mamm <- read.csv("White_Mamm.csv")
squa <- read.csv("White_Squa.csv")

#load the trees

fishTree <- read.tree("fish.tre")
amphTree <- read.tree("amph.tre")
birdTree <- read.tree("bird.tre")
mammTree <- read.tree("mamm.tre")
squaTree <- read.tree("squa.tre")

#confirm all names match. ugly function with no checks in it for different formatting 
#issues. that that it expects a list of lists, where first set of lists are the data
#frames (with a column "species") and second are the trees.

checkTaxa <- function(name.tree.list)
{
	output <- list()
	for(i in 1:length(name.tree.list[[1]]))
	{
		output[[i]] <- setdiff(name.tree.list[[1]][[i]]$species, name.tree.list[[2]][[i]]$tip.label)
	}
	names(output) <- names(name.tree.list[[1]])
	return(output)
}

nameList <- list(fish, amph, bird, mamm, squa)
names(nameList) <- c("fish","amph","bird","mamm","squa")
treeList <- list(fishTree, amphTree, birdTree, mammTree, squaTree)
nameTreeList <- list(nameList, treeList)

checked <- checkTaxa(nameTreeList)

#ok all checks out. notice that the bird database is in a different format. get rid of the
#unnecessary columns to facilitate binding

bird <- select(bird, species, mass, temp, MR, mean.mass, q10smr)

#now just rbind them all, then add a column describing what group it is

traits <- rbind(fish, amph, bird, mamm, squa)

traits$group <- c(rep("fish", length(fish$species)), rep("amph", length(amph$species)), 
	rep("bird", length(bird$species)), rep("mamm", length(mamm$species)),
	rep("squa", length(squa$species)))

traits$color <- c(rep("#1B9E77", length(fish$species)), rep("#D95F02", length(amph$species)), 
	rep("#7570B3", length(bird$species)), rep("#E7298A", length(mamm$species)),
	rep("#66A61E", length(squa$species)))

#take out any where mean mass is na. Once you do this, there are no entries where q10smr 
#is na either

final <- traits[!is.na(traits$mean.mass),]

#write.csv(final, "final_white_etal.csv", row.names=FALSE)

plot(log(final$q10smr)~log(final$mean.mass), col=final$color, xlab="Ln(body mass)", 
	ylab="Ln(basal metabolic rate)", pch=20)

legend("topleft", c("Fishes", "Amphibians", "Birds", "Mammals", "Squamates"), pch=20, 
	col = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E"))



