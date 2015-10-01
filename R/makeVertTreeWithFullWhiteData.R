## Stitch together a full tree from all taxa
require(phytools)
require(geiger)
require(aRbor)
setwd("~/repos/bmr/R/")
#amphib <- readRDS("../output/data/amphibians.rds")
#birds <- readRDS("../output/data/mckenchie.rds")
#squam <- readRDS("../output/data/squamates.rds")
#mamm <- readRDS("../output/data/newmammals.rds")
#fish <- readRDS("../output/data/fish.rds")
allwhite <- read.csv("../datasets/all_white_etal.csv")
k <- 8.617332478E-5

tmp <- filter_(allwhite, paste("!is.na(mass)", sep=""), paste("!is.na(temp)", sep=""), paste("!is.na(MR)", sep=""))
lnBMR <- log(tmp$MR)
lnMass <- log(tmp$mass)
temp <- 1/k*1/(tmp$temp+273.15)
mod1 <- lm(lnBMR~lnMass+temp+tmp$species)
Ei <- mod1$coef[3]
uspecies <- unique(allwhite$species)
ndat <- lapply(1:length(uspecies), function(x){tmp <- filter_(allwhite, paste("species=='", uspecies[x], "'", sep=""),
                                                      paste("!is.na(mass)", sep=""), paste("!is.na(temp)", sep=""), paste("!is.na(MR)", sep=""));
                                       data.frame("species"=uspecies[x], "lnMass"=mean(log(tmp$mass)), 
                                                 "lnBMR"= mean(log(tmp$MR))+Ei/k*(mean(1/(tmp$temp+273.15))-1/T2))})
ndat <- do.call(rbind, ndat)
rownames(ndat) <- NULL

tds <- list(amphib=amphib, birds=birds, squam=squam, mamm=mamm, fish=fish)
tds <- lapply(tds, function(x) make.treedata(x$phy, as.data.frame(x$dat))) 


ntips <- sum(sapply(trees, function(x) length(x$phy$tip.label)))
tip.labels <- c("fish", "amphib", "squam", "birds", "mamm")
edge <- matrix(c(9, 4,
  9, 3,
  8, 5,
  8, 9,
  7, 8,
  7, 2,
  6, 7,
  6, 1), byrow=TRUE, ncol=2)
edge.length <- c(274.9, 274.9, 324.5, 324.5-274.9, 382.9-324.5, 382.9, 454.6-382.9 , 454.6)
Nnode <- 4
ordertree <- list(edge=edge, Nnode=Nnode, tip.label=tip.labels, edge.length=edge.length)
class(ordertree) <- 'phylo'
ordertree <- reorder(ordertree, "postorder")
plot(ordertree)

otax <- data.frame("Class"= ordertree$tip.label, "Superclass"=c("Actinopterygii", rep("Tetrapoda",4)))
rownames(otax) <- ordertree$tip.label
classtree <- nodelabel.phylo(ordertree, otax, ncores=1)

trees <- lapply(tds, function(x) x$phy)
trees <- lapply(trees, multi2di)
class(trees) <- "multiPhylo"
plot(classtree)
abline(v=sapply(trees, function(x) max(nodeHeights(x))),lty=2)
res <- glomogram.phylo(classtree, trees)

tdsall <- lapply(tds, function(x) select(x, mean.mass, q10smr))
tdsall <- lapply(tdsall, function(x) mutate(x, lnMass=log(mean.mass), lnBMR=log(q10smr)))
dats <- lapply(tdsall, function(x) x$dat)
rn <- unname(unlist(sapply(tdsall, function(x) rownames(x$dat))))
resdat <- do.call(rbind, dats)
rownames(resdat) <- rn
td <- make.treedata(res, resdat)
plot(td$dat$lnMass, td$dat$lnBMR)

saveRDS(td, "../output/data/tetrapods.rds")

td <- reorder(td, "postorder")
tmp <- identifyBranches(td$phy, 30, plot.simmap=TRUE)


glomogram.phylo(classtree, trees[4])

plot(res)


