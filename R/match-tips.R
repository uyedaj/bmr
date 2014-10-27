library(geiger)

dir.create("output/data")

remove_dups <- function(x){
    x <- x[complete.cases(x),]
    x[,"species"] <- as.character(x[,"species"])
    x <- x[!duplicated(x$species),]
    rownames(x) <- x$species
    x
}

dat.mamm <- read.csv("datasets/White_Mamm.csv")
dat.mamm <- remove_dups(dat.mamm)

phy.mamm <- read.tree("datasets/mamm.tre")
td.mamm  <- treedata(phy.mamm, dat.mamm)
saveRDS(td.mamm, "output/data/mammals.rds")

dat.squa <- read.csv("datasets/White_Squa.csv")
dat.squa <- remove_dups(dat.squa)

phy.squa <- read.tree("datasets/squa.tre")
td.squa  <- treedata(phy.squa, dat.squa)
saveRDS(td.squa, "output/data/squamates.rds")

dat.amph <- read.csv("datasets/White_Amph.csv")
dat.amph <- remove_dups(dat.amph)

phy.amph <- read.tree("datasets/amph.tre")
td.amph  <- treedata(phy.amph, dat.amph)
saveRDS(td.squa, "output/data/amphibians.rds")

dat.fish <- read.csv("datasets/White_Fishes.csv")
dat.fish <- remove_dups(dat.fish)

phy.fish <- read.tree("datasets/fish.tre")
td.fish  <- treedata(phy.fish, dat.fish)
saveRDS(td.squa, "output/data/fish.rds")

dat.bird <- read.csv("datasets/Whiteetal_bird.csv")
dat.bird <- remove_dups(dat.bird)

phy.bird <- read.tree("datasets/bird.tre")
td.bird  <- treedata(phy.bird, dat.bird)
saveRDS(td.bird, "output/data/birds.rds")





