## Stitch together a full tree from all taxa
require(phytools)
require(geiger)
require(aRbor)
setwd("~/repos/bmr/R/")
allwhite <- read.csv("../datasets/all_white_etal.csv")
attributes(allwhite$iT) <- NULL

uSpecies <- unique(allwhite$species)
nSpecies <- tapply(allwhite$species, allwhite$species, length)
spp <- names(nSpecies)[which(nSpecies > 3)]

lmm <- function(species){
  tmp <- filter_(allwhite, paste("species=='",species,"'", sep=""), paste("!is.na(mass)", sep=""), paste("!is.na(temp)", sep=""), paste("!is.na(MR)", sep=""))
  if(nrow(tmp) < 3){return("error")}
  lnBMR <- log(tmp$MR)
  lnMass <- log(tmp$mass)
  temp <- 1/(293.15+tmp$temp)*1/k
  if(length(unique(lnMass)) ==1 | length(unique(temp))==1){
    return("error")
  } else {
    lmfit <- lm(lnBMR ~ lnMass + temp)
  }
  return(lmfit)
}

tmp <- filter_(allwhite, paste("!is.na(mass)", sep=""), paste("!is.na(temp)", sep=""), paste("!is.na(MR)", sep=""))
lnBMR <- log(tmp$MR)
lnMass <- log(tmp$mass)
temp <- 1/k*1/(tmp$temp+273.15)
mod1 <- lm(lnBMR~lnMass+temp+tmp$species)


lmfits <- lapply(spp, lmm)
names(lmfits) <- spp
lmfits <- lmfits[!(lmfits=="error")]
sumfits <- lapply(lmfits, summary)

tempCoef <- lapply(sumfits, function(x) x$coef[3, 1:2])
tempCoefs <- do.call(rbind, tempCoef)
o <- order(tempCoefs[,2])
tempCoefs <- tempCoefs[o,]

massCoef <- lapply(sumfits, function(x) x$coef[2,1:2])
massCoefs <- do.call(rbind, massCoef)
massCoefs <- massCoefs[o,]
coefs <- cbind(tempCoefs, massCoefs)
colnames(coefs) <- c("temp", "temp.se", "mass", "mass.se")

par(mfrow=c(1,2))
plot(massCoefs[,1], pch=21, bg=1, ylim=c(0, 2))
lapply(1:nrow(massCoefs), function(x) lines(c(x,x), c(massCoefs[x,1]-2*massCoefs[x,2],massCoefs[x,1]+2*massCoefs[x,2] ), lty=2))
abline(h=c(0.66666667, 0.75))

plot(tempCoefs[,1], pch=21, bg=1, ylim=c(-2, 0.5))
lapply(1:nrow(tempCoefs), function(x) lines(c(x,x), c(tempCoefs[x,1]-2*tempCoefs[x,2],tempCoefs[x,1]+2*tempCoefs[x,2] ), lty=2))
abline(h=-0.65)

tdall <- readRDS("../output/data/tetrapods.rds")
vertree <- tdall$phy
tdtmp <- make.treedata(vertree, coefs)
phenogram(tdtmp$phy, setNames(tdtmp$dat[[1]], tdtmp$phy$tip.label), fsize=0.25)
phenogram(tdtmp$phy, setNames(tdtmp$dat[[3]], tdtmp$phy$tip.label), fsize=0.25)

fctemp <- fitContinuous(tdtmp$phy, setNames(tdtmp$dat[[1]], tdtmp$phy$tip.label), SE =  setNames(tdtmp$dat[[2]], tdtmp$phy$tip.label), model="OU")
fctemp
log(2)/fctemp$opt$alpha

fctemp <- fitContinuous(tdtmp$phy, setNames(tdtmp$dat[[3]], tdtmp$phy$tip.label), SE =  setNames(tdtmp$dat[[4]], tdtmp$phy$tip.label), model="OU")
fctemp
log(2)/fctemp$opt$alpha


mean(tempCoefs[,1])




