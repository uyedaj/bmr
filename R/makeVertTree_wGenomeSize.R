## Stitch together a full tree from all taxa
require(phytools)
require(geiger)
require(aRbor)
setwd("~/repos/bmr/R/")
amphib <- readRDS("../output/data/amphibians.rds")
birds <- readRDS("../output/data/birds.rds")
squam <- readRDS("../output/data/squamates.rds")
mamm <- readRDS("../output/data/newmammals.rds")
fish <- readRDS("../output/data/fish.rds")

corrected <- read.csv("../datasets/BMR_corrected.csv")
to.drop <- corrected$taxon[which(corrected$action=="cut")]
to.replace <- list(taxa=corrected$taxon[which(corrected$action=="replace")], nos=which(corrected$action=="replace"))

tds <- list(amphib=amphib, birds=birds, squam=squam, mamm=mamm, fish=fish)
tds <- lapply(tds, function(x) make.treedata(x$phy, as.data.frame(x$dat), name_column=0))
trees <- lapply(tds, function(x) x$phy)


drop.subspecies <- function(species){
  tmp <- lapply(species, function(x) strsplit(x, "_")[[1]])
  sp <- sapply(tmp, function(x) paste(x[1], x[2], sep="_"))
  return(sp)
}

genomesize <- read.csv("../datasets/all_genome_size.csv", fileEncoding="latin1")
genomesize$Species <- sapply(genomesize$Species, function(x) gsub(" ", "_", x))
genomesize$C.value <- as.character(genomesize$C.value)
genomesize$C.value[grep("-", genomesize$C.value)] <- sapply(genomesize$C.value[grep("-", genomesize$C.value)], function(x) mean(as.numeric(unlist(strsplit(as.character(x), "-")))))
cvalues <- tapply(as.numeric(genomesize$C.value), drop.subspecies(genomesize$Species), mean, na.rm=TRUE)
cvalues <- cbind(cvalues)
colnames(cvalues) <- "cvalue"
dats <- lapply(tds, function(x) cbind(x$dat, cvalue=cvalues[match(x$phy$tip.label, rownames(cvalues)),]))

tds <- lapply(1:length(tds), function(x) make.treedata(tds[[x]]$phy, dats[[x]], name_column=0))

#tdgs <- lapply(1:length(tds), function(x) make.treedata(tds[[x]]$phy, cvalues))
#sapply(tdgs, function(x) length(x$phy$tip.label))/sapply(tds, function(x) length(x$phy$tip.label))
#sum(sapply(tdgs, function(x) length(x$phy$tip.label)))/sum(sapply(tds, function(x) length(x$phy$tip.label)))

ntips <- sum(sapply(trees, function(x) length(x$tip.label)))
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
names(tds)<- names(trees) <- c("amphib", "birds", "squam", "mamm", "fish")
res <- glomogram.phylo(classtree, trees)

tdsall <- lapply(tds, function(x) select(x, mean.mass, q10smr, cvalue))
tdsall <- lapply(tdsall, function(x) mutate(x, lnMass=log(as.numeric(mean.mass)), lnBMR=log(as.numeric(q10smr)), lnGenSize=log(as.numeric(cvalue))))
dats <- lapply(tdsall, function(x) x$dat)
rn <- unname(unlist(sapply(tdsall, function(x) rownames(x$dat))))
resdat <- do.call(rbind, lapply(dats, function(x) x[,4:6]))
rownames(resdat) <- rn
td <- make.treedata(res, resdat)
td <- treeply(td, drop.tip, as.character(to.drop))
td$dat[["lnMass"]][match(to.replace$taxa, td$phy$tip.label)] <- log(corrected$newMass..g.[to.replace$nos])
td$dat[["lnBMR"]][match(to.replace$taxa, td$phy$tip.label)] <- log(corrected$newRate..standardized.to.38.[to.replace$nos])

plot(td$dat$lnMass, td$dat$lnBMR)
td <- reorder(td, "postorder")
saveRDS(td, "../output/data/tetrapods_gs.rds")

#tmp <- identifyBranches(td$phy, 30, plot.simmap=TRUE)

#reps <- sapply(tds, function(x) nrow(x$dat))
#cladeid <- factor(unlist(sapply(1:5, function(x) rep(names(reps[x]), reps[x]))))
#thermy <- factor(ifelse(cladeid %in% c("mamm", "birds"), "endo", "ecto"))
#tmp2 <- lm(td$dat$lnBMR ~ td$dat$lnMass+I(td$dat$lnMass^2)+td$dat$lnGenSize+thermy)
#summary(tmp2)
#plot(tmp2$resid)

#glomogram.phylo(classtree, trees[4])

#plot(res)


