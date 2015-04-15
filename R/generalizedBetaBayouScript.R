## This script runs the most complicated fixed shift model that we are going to run

## Load in command args:
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
print(args)

## Uncomment out the following lines to run locally:
rm(list=ls(all=TRUE))
setwd("~/repos/bmr/R/")
args <- list("tetrapods_gs", 10000, "delete_me")

## Future switches....These do nothing at the moment:
pred_vars <- c("lnmass", "lnMass2", "lnGS") #Options: lnMass, lnMass2, lnGS, endo
taxa_shifts <- c("Tetrapods", "Amniota","Varanus_exanthematicus", "Pseudoeurycea_brunnata", "Sceloporus_variabilis", "Sceloporus_occidentalis", "Uta_stansburiana", "Masticophis_flagellum", "Typhlogobius_californiensis", "Sebastolobus_altivelis", "Lacertidae", "Boidae", "Hyla_arborea", "Aves", "Mammals", "Lichanura_trivirgata", "Bunopus_tuberculatus", "Fishes", "Scaphiophidae", "Euphlyctis", "Labeo", "Salmonidae", "Cricetidae", "Ambystoma_mexicanum")
  #Options: c(Tetrapods, Amniota, Varanus_exanthematicus, Pseudoeurycea_brunnata, Sceloporus_variabilis, Sceloporus_occidentalis, Uta_stansburiana, Masticophis_flagellum, Typhlogobius_californiensis, Sebastolobus_altivelis, Lacertidae, Boidae, Hyla_arborea, Aves, Mammals, Lichanura_trivirgata, Bunopus_tuberculatus, Fishes, Scaphiophidae, Euphlyctis, Labeo, Salmonidae, Cricetidae, Ambystoma_mexicanum
missing_data <- "impute" #Options: impute, exclude

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
td <- make.treedata(tree, dat) # Rematch data and tree
td <- reorder(td, "postorder") # Reorder tree
td_gs <- filter(td, !is.na(lnGenSize)) # Produce a 2nd dataset of complete genome size data
tree <- td$phy
dat <- td$dat
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
custom.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  pred[is.na(pred[,3]),3] <- pars$missing.pred #$impute
  betaID <- getTipMap(pars, cache)
  ## Specify the model here
  X = X - pars$beta1[betaID]*pred[,1] - pars$beta2[betaID]*pred[,2] - pars$beta3*pred[,3]
  cache$dat <- X
  ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
  dpars <- pars
  dpars$theta[dpars$t2[which(dpars$sb %in% c(841, 1703))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(841, 1703))]]+dpars$endo
  ### The part below mostly does not change
  X.c <- bayou:::C_weightmatrix(cache, dpars)$resid
  transf.phy <- bayou:::C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
  transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
  comp <- bayou:::C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
  if(pars$alpha==0){
    inv.yVy <- comp$PP
    detV <- comp$logd
  } else {
    inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
    detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
  }
  llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
  llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

## List of predetermined shifts with high support
sb <- c("Tetrapods"=1707, "Amniota"=1705,"Varanus_exanthematicus"=563, "Pseudoeurycea_brunnata"=277, "Sceloporus_variabilis"=515, "Sceloporus_occidentalis"=509, "Uta_stansburiana"=518, "Masticophis_flagellum"=443, "Typhlogobius_californiensis"=54, "Sebastolobus_altivelis"=85, "Lacertidae"=614, "Boidae"=501, "Hyla_arborea"=159, "Aves"=841, "Mammals"=1703, "Lichanura_trivirgata"=495, "Bunopus_tuberculatus"=655, "Salamandridae"=378, "Fishes"=115, "Scaphiophidae"=261, "Euphlyctis"=227, "Labeo"=7, "Bolitoglossinae"=294, "Salmonidae"=48, "Plethodontidae"=344, "Cricetidae"=1222, "Ambystoma_mexicanum"=368)
k <- length(sb)

## Define starting parameters 
startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(k+1, 0.7, 0.1), beta2=rnorm(k+1, 0, 0.01), beta3=rnorm(1, 0, 0.05), endo=2, k=k, ntheta=k+1, theta=rnorm(k+1, 0, 1), sb=sb, loc=rep(0, k), t2=2:(k+1))

## Get optimized starting values and trait evolutionary models for genome size. 
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

## This is a function to monitor the output of bayou for our custom model
## Optional:
BetaBMR.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1", "endo", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$endo, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.BetaBMR <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                   beta1=".vectorMultiplier", beta2=".vectorSlidingWindow",
                                   beta3=".slidingWindowProposal", endo=".slidingWindowProposal", 
                                   theta=".adjustTheta", slide=".slide", 
                                   pred.sig2=".multiplierProposal" , pred.root=".slidingWindowProposal", #$impute 
                                   missing.pred=".imputePredBM" #$impute
                                   ),
                      control.weights = list(alpha=4, sig2=2, beta1=10, 
                                             beta2=8, beta3=4, endo=3, 
                                             theta=10, slide=2, k=0,
                                             pred.sig2=1, pred.root=1, missing.pred=3 #$impute
                                             ),
                      D = list(alpha=1, sig2= 0.75, beta1=0.75, beta2=0.05, beta3=0.05, endo=0.25, theta=2, slide=1, 
                               pred.sig2=1, pred.root=1, missing.pred=1 #$impute
                               ),
                      parorder = names(startpar),
                      rjpars = c("theta"),
                      shiftpars = c("sb", "loc", "t2"),
                      monitor.fn = BetaBMR.monitor,
                      lik.fn = custom.lik)


## Now we define the prior:
prior <- make.prior(tree, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta1="dnorm",
                               dbeta2="dnorm", dbeta3="dnorm", dendo="dnorm", 
                               dsb="fixed", dk="fixed", dtheta="dnorm",
                               dpred.sig2="dhalfcauchy", dpred.root="dnorm" #$impute
                               ), 
                    param=list(dalpha=list(scale=1), dsig2=list(scale=1), dbeta1=list(mean=0.7, sd=0.1), 
                               dbeta2=list(mean=0, sd=0.05),dbeta3=list(mean=0, sd=0.05), 
                               dendo=list(mean=0, sd=4),  dk="fixed", dsb="fixed", 
                               dtheta=list(mean=0, sd=4),
                               dpred.sig2=list(scale=1), dpred.root=list(mean=1, sd=1) #$impute
                               ), 
                    fixed=list(sb=startpar$sb, k=startpar$k)
                    )

tr <- pars2simmap(startpar, cache$phy)
plotSimmap(tr$tree, colors=tr$col, fsize=0.5)

## Test to make sure shit works
prior(startpar)
custom.lik(startpar, cache, cache$dat)$loglik
custom.lik(.imputePredBM(cache, startpar, 1, move=NULL)$pars, cache, cache$dat)$loglik #$impute


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

