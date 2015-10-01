## This script runs the most complicated fixed shift model that we are going to run

## Load in command args:
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
print(args)

## Uncomment out the following lines to run locally:
rm(list=ls(all=TRUE))
setwd("~/repos/bmr/R/")
args <- list("tetrapods_raw", 10000, "delete_me", "rj.shift_lnMass.fixed_endo.impute_TempQ10")

## Future switches....These do nothing at the moment:
#pred_vars <- c("lnmass", "lnMass2", "lnGS") #Options: lnMass, lnMass2, lnGS, endo
#taxa_shifts <- c("Tetrapods", "Amniota","Varanus_exanthematicus", "Pseudoeurycea_brunnata", "Sceloporus_variabilis", "Sceloporus_occidentalis", "Uta_stansburiana", "Masticophis_flagellum", "Typhlogobius_californiensis", "Sebastolobus_altivelis", "Lacertidae", "Boidae", "Hyla_arborea", "Aves", "Mammals", "Lichanura_trivirgata", "Bunopus_tuberculatus", "Fishes", "Scaphiophidae", "Euphlyctis", "Labeo", "Salmonidae", "Cricetidae", "Ambystoma_mexicanum")
  #Options: c(Tetrapods, Amniota, Varanus_exanthematicus, Pseudoeurycea_brunnata, Sceloporus_variabilis, Sceloporus_occidentalis, Uta_stansburiana, Masticophis_flagellum, Typhlogobius_californiensis, Sebastolobus_altivelis, Lacertidae, Boidae, Hyla_arborea, Aves, Mammals, Lichanura_trivirgata, Bunopus_tuberculatus, Fishes, Scaphiophidae, Euphlyctis, Labeo, Salmonidae, Cricetidae, Ambystoma_mexicanum
#missing_data <- "impute" #Options: impute, exclude
modelCode <- args[[4]]
## Other example modelCodes: "rj.shift_lnMass.fixed_endo.missing_drop", "rj.shift_lnMass_lnMass2.fixed_endo_lnGS.missing_impute"; "fixed.shift_lnMass.fixed_endo_lnMass2.missing_drop"


## Flags: If a line is only used for imputation (genome size), it will be flagged with $impute

## Load packages. If you haven't installed bayou from github, uncomment out that line. 
require(devtools)
#install_github("uyedaj/bayou", ref="dev")
require(ape)
require(aRbor)
require(mvnfast)
require(Matrix)
require(bayou)
source("./betaBayouFunctions.R")
load_all("~/repos/bayou/bayou_1.0")

## Match and prepare the datasets
td <- readRDS(paste("../output/data/", args[[1]], ".rds", sep=""))
tree <- td$phy
dat <- td$dat
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches
rownames(dat) <- attributes(td)$tip.label
td <- make.treedata(tree, dat) # Rematch data and tree
td <- reorder(td, "postorder") # Reorder tree
td_gs <- filter(td, !is.na(lnGenSize)) # Produce a 2nd dataset of complete genome size data
tree <- td$phy
dat <- td$dat
rownames(dat) <- attributes(td)$tip.label
## BM likelihood function for genome size
gs.lik <- bm.lik(td_gs$phy, setNames(td_gs$dat[[3]], td_gs$phy$tip.label), SE=0, model="BM")

lnBMR <- setNames(dat[['lnBMR']], tree$tip.label)
lnMass <- setNames(dat[['lnMass']], tree$tip.label)
pred <- cbind(setNames(dat[['lnMass']], tree$tip.label), setNames(dat[['lnMass']]^2, tree$tip.label),setNames(dat[['lnGenSize']], tree$tip.label))
colnames(pred) <- c("lnMass", "lnMass2", "lnGS")

## Create a bayou cache object
cache <- bayou:::.prepare.ou.univariate(tree, setNames(dat$lnBMR, tree$tip.label), pred = pred)

## Calculate values that will be useful for quickly imputing
missing <- which(is.na(cache$pred[,3])) #$impute
pv <- getPreValues(cache) #$impute

## We're going to define a custom likelihood function; this is nearly identical to bayou.lik (bayou's standard
## OU function), except we use pars$beta1 to calculate the residuals after accounting for lnMass.
## 1. Define likelihood function
custom.lik <- liks[[modelCode]]

## List of predetermined shifts with high support
if(substr(modelCode,1,2)=="rj"){
  sb <- NULL
  k <- rpois(1, lambda=30)
} else {
  sb <- c("Tetrapods"=1707, "Amniota"=1705,"Varanus_exanthematicus"=563, "Pseudoeurycea_brunnata"=277, "Sceloporus_variabilis"=515, "Sceloporus_occidentalis"=509, "Uta_stansburiana"=518, "Masticophis_flagellum"=443, "Typhlogobius_californiensis"=54, "Sebastolobus_altivelis"=85, "Lacertidae"=614, "Boidae"=501, "Hyla_arborea"=159, "Aves"=841, "Mammals"=1703, "Lichanura_trivirgata"=495, "Bunopus_tuberculatus"=655, "Salamandridae"=378, "Fishes"=115, "Scaphiophidae"=261, "Euphlyctis"=227, "Labeo"=7, "Bolitoglossinae"=294, "Salmonidae"=48, "Plethodontidae"=344, "Cricetidae"=1222, "Ambystoma_mexicanum"=368)
  k <- NULL
}


## Define starting parameters 
startpar <- startpars[[modelCode]](sb, k)

## Get optimized starting values and trait evolutionary models for genome size. 
if(identical(as.numeric(grep("impute", modelCode)),1)){
  require(phylolm)#$impute 
  tdgs <- filter(td, !is.na(lnGenSize))#$impute 
  fits <- lapply(c("BM", "OUrandomRoot", "OUfixedRoot", "EB"), function(x) phylolm(lnGenSize~1, data=tdgs$dat, phy=tdgs$phy, model=x))#$impute 
  aics <- sapply(fits, function(x) x$aic)#$impute 
  bestfit <- fits[[which(aics == min(aics))]]#$impute 
  phenogram(tdgs$phy, setNames(tdgs$dat$lnGenSize, tdgs$phy$tip.label), fsize=0.5)#$impute 

  ## Set the starting imputation parameters at the ML estimate.
  startpar$pred.root <- unname(bestfit$coeff)#$impute 
  startpar$pred.sig2 <- unname(bestfit$sigma2)#$impute 
  startpar <- .imputePredBM(cache, startpar, d=1, NULL, ct=NULL, prevalues=pv)$pars#$impute 
}

## This is a function to monitor the output of bayou for our custom model
## Optional:
BetaBMR.monitor <- monitors[[modelCode]] 

## We're going to define a custom model with variable slopes (beta1)
model.BetaBMR <- models[[modelCode]](startpar, BetaBMR.monitor, custom.lik)


## Now we define the prior:
prior <- priors[[modelCode]]

tr <- pars2simmap(startpar, cache$phy)
plotSimmap(tr$tree, colors=tr$col, fsize=0.5)

## Test to make sure shit works
prior(startpar)
custom.lik(startpar, cache, cache$dat)$loglik

mymcmc <- bayou.makeMCMC(tree, lnBMR, pred=pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)

mymcmc$run(args[[2]])


chain <- mymcmc$load()
chain <- set.burnin(chain, 0.3)
out <- summary(chain)
print(out)

require(foreach)
require(doParallel)
registerDoParallel(cores=6)
Bk <- seq(0, 1, length.out=6)
ss <- mymcmc$steppingstone(args[[2]], chain, Bk, burnin=0.3, plot=TRUE)
plot(ss)

saveRDS(chain, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".chain.rds",sep=""))
saveRDS(mymcmc, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".mcmc.rds",sep=""))
saveRDS(ss, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".ss.rds",sep=""))

#ss <- readRDS(file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".ss.rds",sep=""))

#plotSimmap.mcmc(tree, chain, burnin=0.3)


#chain <- mymcmc$load()
#out <- list(tree=mymcmc$tree, dat=mymcmc$dat, outname="newmammals_quadr1", model="custom", model.pars=mymcmc$model.pars, dir="~/repos/bmr/output/runs/newmammals/")
#chain <- load.bayou(out,save.Rdata = TRUE, file="../output/runs/newmammals/newmammals_quad.rds")


#pdf("../output/figures/QuadraticMammals.pdf")
#plotSimmap.mcmc(tree, chain, fsize=0.25, burnin=0.3)
#postburn <- round(0.3*length(chain$gen), 0):length(chain$gen)
#plot(density(sapply(chain$beta2[postburn], function(x) x[1])))#

#par(mfrow=c(2,1))
#plot(c(0, length(chain$gen)), c(0.5, 1), type="n", xlab="Sample", ylab="beta1")
#dum <- lapply(1:length(chain$gen), function(x) points(rep(x, length(chain$beta1[[x]])), chain$beta1[[x]], pch=".", col=1:length(chain$beta2[[x]])))
#plot(c(0, length(chain$gen)), c(-0.03, 0.03), type="n", xlab="Sample", ylab="beta2")
#dum <- lapply(1:length(chain$gen), function(x) points(rep(x, length(chain$beta2[[x]])), chain$beta2[[x]], pch=".", col=1:length(chain$beta2[[x]])))

#par(mfrow=c(1,1))
#samp <- round(seq(postburn[1], postburn[length(postburn)], length.out=1000),0)

#{burnin=0.3
#LP <- Lposterior(chain, tree, burnin=burnin)
#focalSB <- c(0, which(LP[,1] > cutoff))
#PostSB <- lapply(focalSB, function (y) getBranchPosterior(y, chain, burnin=burnin, thin=thin))
#postBetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[1]]))
#postThetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[2]]))
#postBetas <- data.frame(do.call(rbind, postBetas))
#postBetas$var <- rep("Beta", nrow(postBetas))
#postThetas <- data.frame(do.call(rbind, postThetas))
#postThetas$var <- rep("Theta", nrow(postThetas))
#allpar <- rbind(postBetas, postThetas)
#allpar$x <- factor(allpar$x)
#focalpars <- list(k=length(focalSB[-1]), 
#                  ntheta=length(focalSB), 
#                  sb=focalSB[-1], 
#                  loc=rep(0, length(focalSB[-1])), 
#                  t2=2:length(focalSB))#

#tr <- pars2simmap(focalpars, tree)
#tipreg <- bayou:::.tipregime(focalpars, tree)

#plot(lnMass, lnBMR, type="n")
#beta1s <- chain$beta1[samp]
#beta2s <- chain$beta2[samp]
#thetas <- chain$theta[samp]
#lapply(1:length(samp), function(j) sapply(1:length(beta1s[[j]]), function(y) curve(thetas[[j]][y]+beta1s[[j]][y]*x + beta2s[[j]][y]*x^2, add=TRUE, col=makeTransparent(y, 5))))
#points(lnMass, lnBMR, pch=21, bg=makeTransparent(tipreg,alpha=50), col=makeTransparent(tipreg,alpha=50))#

#}
#dev.off()

