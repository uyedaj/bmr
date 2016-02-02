## Making Kolokotrones dataset
## Load packages. If you haven't installed bayou from github, uncomment out that line. 
require(devtools)
#install_github("uyedaj/bayou", ref="dev")
require(ape)
require(aRbor)
require(mvnfast)
require(Matrix)
require(bayou)
source("./betaBayouFunctions.R")
#load_all("~/repos/bayou/bayou_1.0")

setwd("~/repos/bmr/R/")
td <- readRDS("../output/data/tetrapods_ei.rds")
td_gs <- readRDS("../output/data/tetrapods_gs.rds")
lnGenSize <- setNames(td_gs$dat$lnGenSize, td_gs$phy$tip.label)
lnGenSize <- lnGenSize[!is.na(lnGenSize)]

tree <- td$phy
dat <- td$dat

koko <-  readxl::read_excel("../datasets/koko.xls", skip=6)
pan <- read.csv("../datasets/PanTHERIA_1-0_WR05_Aug2008.csv")
pan$genspec <- gsub(" ", "_", pan$MSW05_Binomial)
pan <- select(pan, genspec, X5.1_AdultBodyMass_g)
colnames(pan)[2] <- "lnMass.pan"
pan[,2] <- log(as.numeric(as.character(pan[,2])))

koko$genspec <- gsub(" ", "_" , koko[["Genus Species"]])
colnames(koko) <- gsub(" ", ".", colnames(koko)); colnames(koko) <- gsub(")", "", colnames(koko), fixed=TRUE); colnames(koko) <- gsub("(", "", colnames(koko), fixed=TRUE); colnames(koko) <- gsub("/", "", colnames(koko), fixed=TRUE)
koko <- koko[,!duplicated(colnames(koko))]
koko <- mutate(koko, lnBMR.k = log(BMR.kJh), lnMass.k = log(Mass.g))
merged <- left_join(koko, data.frame(genspec=tree$tip.label, dat), by="genspec")
merged <- left_join(merged, pan)
plot(merged$lnBMR.k, merged$lnBMR)

fulltree <- read.tree("../datasets/tetrapods.tre")
#nH <- phytools::nodeHeights(fulltree)
#saveRDS(nH, file='../output/data/nodeHeights.rds')
nH <- readRDS(file="../output/data/nodeHeights.rds")
external <- which(fulltree$edge[,2] <= length(fulltree$tip.label))
fulltree$edge.length[external] <- fulltree$edge.length[external] + (max(nH[,2]) - nH[external,2])

phynd <- phyndr::phyndr_genus(fulltree, koko$genspec)
saveRDS(phynd, file="../output/data/kolokoPhynd.rds")
phynd <- readRDS("../output/data/kolokoPhynd.rds")
tdkoko <- make.treedata(phyndr::phyndr_sample(phynd), merged)
tdkoko <- select(tdkoko, lnMass.k, lnBMR.k, lnMass, lnBMR, lnMass.pan, Temperature.C)
tdkoko <- filter(tdkoko, !is.na(lnMass.k), !is.na(lnBMR.k))
saveRDS(tdkoko, file = "../output/data/kolokotrones.rds")
plot(tdkoko$dat$lnMass.k, tdkoko$dat$lnBMR.k)
