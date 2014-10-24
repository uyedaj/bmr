#this is a quick script to clean the names in the appendix from McKechnie & Wolf

library(ape)

tree <- read.tree("bird.tre")

#set strings as factors to false to avoid problems later

traits <- read.csv("McKechnie_Wolf.csv", stringsAsFactors=FALSE)

#what names are in the trait database that aren't in the tree?

examine <- setdiff(traits$species, tree$tip.label)

#send the wrong names out to excel to check

write.csv(examine, "look_here.csv", row.names=FALSE)

#manually cleaned. bring the cleaned back in, make sure to set no factors

use <- read.csv("McKechnie_Wolf_clean1.csv", stringsAsFactors=FALSE)

#subset to instances where the bad name equals the original name, and replace with cleaned

traits$species[traits$species %in% use$original] <- use$cleaned

#now check to see if there are any additional problems

examine <- setdiff(traits$species, tree$tip.label)

write.csv(examine, "look_here.csv", row.names=FALSE)

#bring in the twice cleaned data

use <- read.csv("McKechnie_Wolf_clean2.csv", stringsAsFactors=FALSE)

traits$species[traits$species %in% use$original] <- use$cleaned

#ok everything matches. save this file out

write.csv(traits, "McKechnie_Wolf.csv", row.names=FALSE)