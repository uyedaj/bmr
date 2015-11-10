## This script runs the most complicated fixed shift model that we are going to run

## Load in command args:
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
print(args)

## Uncomment out the following lines to run locally:
#rm(list=ls(all=TRUE))
#setwd("~/repos/bmr/R/")
#args <- list("kolokotrones", 1000000, "_koko2")

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
#load_all("~/repos/bayou/bayou_1.0")

## Match and prepare the datasets
td <- readRDS(paste("../output/data/", args[[1]], ".rds", sep=""))
tree <- td$phy
dat <- td$dat
rownames(dat) <- attributes(td)$tip.label
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches
td <- make.treedata(tree, dat) # Rematch data and tree
td <- reorder(td, "postorder") # Reorder tree
td <- filter(td, !is.na(lnBMR.k), !is.na(lnMass.k))
#td_gs <- filter(td, !is.na(lnGenSize)) # Produce a 2nd dataset of complete genome size data
tree <- td$phy
dat <- td$dat
## BM likelihood function for genome size
#gs.lik <- bm.lik(td_gs$phy, setNames(td_gs$dat[[3]], td_gs$phy$tip.label), SE=0, model="BM")

lnBMR <- setNames(dat[['lnBMR.k']], tree$tip.label)
lnMass <- setNames(dat[['lnMass.k']], tree$tip.label)
.pred <- cbind(setNames(dat[['lnMass.k']], tree$tip.label), setNames(dat[['lnMass.k']]^2, tree$tip.label))#,setNames(dat[['lnGenSize']], tree$tip.label))
colnames(.pred) <- c("lnMass", "lnMass2")#, "lnGS")

## Create a bayou cache object
cache <- bayou:::.prepare.ou.univariate(tree, setNames(dat$lnBMR.k, tree$tip.label), pred = .pred)
pred <- cache$pred
# Birds branch 847; mammals branch 1709
#endotherms <- lapply(cache$desc$tips[cache$phy$edge[c(847, 1709), 2]], function(x) cache$phy$tip.label[x])
#pred <- data.frame(pred, endo=as.numeric(cache$phy$tip.label %in% unlist(endotherms)))
#cache$pred <- pred

## Calculate values that will be useful for quickly imputing
#missing <- which(is.na(cache$pred[,3])) #$impute
#pv <- getPreValues(cache) #$impute

## We're going to define a custom likelihood function; this is nearly identical to bayou.lik (bayou's standard
## OU function), except we use pars$beta1 to calculate the residuals after accounting for lnMass.
## 1. Define likelihood function
custom.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  #pred[is.na(pred[,3]),3] <- pars$missing.pred #$impute
  betaID <- getTipMap(pars, cache)
  ## Specify the model here
  X = X - pars$beta1[betaID]*pred[,1] - pars$beta2*pred[,2]#- pars$endo*pred[,3] #- pars$beta3*pred[,3]
  cache$dat <- X
  ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
  #dpars <- pars
  #dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]]+dpars$endo
  ### The part below mostly does not change
  X.c <- bayou:::C_weightmatrix(cache, pars)$resid
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
  #llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

## List of predetermined shifts with high support
## Since we have the endothermy parameter, there is an identifiability issue with birds and mammals; keeping only one of these shifts
## As an allowable shift ensures that this parameter is the difference between birds and mammals, and the other parameter is common to both.
bmax <- prob <- rep(1, nrow(cache$edge))
#bmax[1709] <- 0  
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0
sb <- NULL
k <- rpois(1, lambda=0.025*sum(tree$edge.length != .Machine$double.eps))
startpars_gen <- function(sb, k) {
  sb <- sample((1:length(cache$bdesc))[-(which(bmax==0))], k, replace=FALSE, prob = sapply(cache$bdesc, length)[-(which(bmax==0))])
  startpar <- list(alpha=5, sig2=5, beta1=rnorm(k+1, 0.7, 0.05), beta2=rnorm(1,0,0.01), k=k, ntheta=k+1, theta=rnorm(k+1, -2.5, 1), sb=sb, loc=rep(0, k), t2=2:(k+1))
  return(startpar)
}


## Define starting parameters 
startpar <- startpars_gen(sb, k)

## Get optimized starting values and trait evolutionary models for genome size. 
#require(phylolm)#$impute 
#tdgs <- filter(td, !is.na(lnGenSize))#$impute 
#fits <- lapply(c("BM", "OUrandomRoot", "OUfixedRoot", "EB"), function(x) phylolm(lnGenSize~1, data=tdgs$dat, phy=tdgs$phy, model=x))#$impute 
#aics <- sapply(fits, function(x) x$aic)#$impute 
#bestfit <- fits[[which(aics == min(aics))]]#$impute 
#phenogram(tdgs$phy, setNames(tdgs$dat$lnGenSize, tdgs$phy$tip.label), fsize=0.5)#$impute 

## Set the starting imputation parameters at the ML estimate.
#startpar$pred.root <- unname(bestfit$coeff)#$impute 
#startpar$pred.sig2 <- unname(bestfit$sigma2)#$impute 
#startpar <- .imputePredBM(cache, startpar, d=1, NULL, ct=NULL, prevalues=pv)$pars#$impute 

## This is a function to monitor the output of bayou for our custom model
## Optional:
BetaBMR.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1", "beta2", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$beta2 , pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.BetaBMR <- {list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                   beta1=".slidingWindowProposal", 
                                   beta2=".slidingWindowProposal", 
                                   theta=".adjustTheta", slide=".slide",
                                   k=".splitmergebd"
                                   #,pred.sig2=".multiplierProposal" , pred.root=".slidingWindowProposal", #$impute 
                                   #missing.pred=".imputePredBM" #$impute
                                   ),
                      control.weights = list(alpha=5, sig2=3, beta1=20, 
                                             beta2=8, #beta3=4, 
                                             #endo=2, 
                                             theta=20, slide=3, k=8#,
                                             #pred.sig2=1, pred.root=1, missing.pred=3 #$impute
                                             ),
                      D = list(alpha=0.5, sig2= 0.5, beta1=0.025, 
                               beta2=0.002, #beta3=0.05, 
                               #endo=0.5, 
                               k=c(4,0.5), theta=2, slide=1 
                               #,pred.sig2=1, pred.root=1, missing.pred=1 #$impute
                               ),
                      parorder = names(startpar),
                      rjpars = c("theta", "beta1"),
                      shiftpars = c("sb", "loc", "t2"),
                      monitor.fn = BetaBMR.monitor,
                      lik.fn = custom.lik)}

model.prime <- {list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                   beta1=".vectorMultiplier", 
                                   beta2=".slidingWindowProposal", 
                                   theta=".adjustTheta", #slide=".slide",
                                   k=".splitmergebd"
                                   #,pred.sig2=".multiplierProposal" , pred.root=".slidingWindowProposal", #$impute 
                                   #missing.pred=".imputePredBM" #$impute
),
control.weights = list(alpha=6, sig2=4, beta1=20, 
                       beta2=8, #beta3=4, 
                       #endo=2, 
                       theta=20, k=0#,
                       #pred.sig2=1, pred.root=1, missing.pred=3 #$impute
),
D = list(alpha=0.5, sig2= 0.5, beta1=0.75, 
         beta2=0.005, #beta3=0.05, 
         #endo=0.5, 
         k=c(3,0.5), theta=2, slide=1 
         #,pred.sig2=1, pred.root=1, missing.pred=1 #$impute
),
parorder = names(startpar),
rjpars = c("theta", "beta1"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = BetaBMR.monitor,
lik.fn = custom.lik)}


## Now we define the prior:

prior <- make.prior(tree, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta1="dnorm",
                               dbeta2="dnorm", #dbeta3="dnorm", 
                               #dendo="dnorm", 
                               dsb="dsb", dk="cdpois", dtheta="dnorm"
                               #,dpred.sig2="dhalfcauchy", dpred.root="dnorm" #$impute
                               ), 
                    param=list(dalpha=list(scale=1), dsig2=list(scale=1), dbeta1=list(mean=0.7, sd=0.3), 
                               dbeta2=list(mean=0, sd=0.01), #dbeta3=list(mean=0, sd=0.05), 
                               #dendo=list(mean=5, sd=1),  
                               dk=list(lambda=26.15, kmax=53), dsb=list(bmax=bmax, prob=1), 
                               dtheta=list(mean=-3, sd=1.5)
                               #,dpred.sig2=list(scale=1), dpred.root=list(mean=1, sd=1) #$impute
                               ), 
                    #fixed=list(sb=startpar$sb, k=startpar$k)
                    )

tr <- pars2simmap(startpar, cache$phy)
plotRegimes(tr$tree, cex=0.5)

## Test to make sure shit works
prior(startpar)
custom.lik(startpar, cache, cache$dat)$loglik
primerMcmc <- bayou.makeMCMC(tree, lnBMR, pred=pred, SE=0, model=model.prime, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]], "primer",sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)
primerMcmc$run(10000)
primerChain <- primerMcmc$load()
startpar <- pull.pars(length(primerChain$gen), primerChain, model.prime)
#custom.lik(.imputePredBM(cache, startpar, 1, move=NULL)$pars, cache, cache$dat)$loglik #$impute

mymcmc <- bayou.makeMCMC(tree, lnBMR, pred=pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)
mymcmc$run(args[[2]])

chain <- mymcmc$load()
chain <- set.burnin(chain, 0.3)
out <- summary(chain)
print(out)

#require(foreach)
#require(doParallel)
#registerDoParallel(cores=6)
#Bk <- seq(0, 1, length.out=6)
#ss <- mymcmc$steppingstone(args[[2]], chain, Bk, burnin=0.3, plot=TRUE)
#plot(ss)

saveRDS(chain, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".chain.rds",sep=""))
saveRDS(mymcmc, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".mcmc.rds",sep=""))
#saveRDS(ss, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".ss.rds",sep=""))

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

sumstats <- summary(chain)
cutoff <- 0.4
sumpars <- list(sb = which(sumstats$branch.posteriors$pp > cutoff))
sumpars$k <- length(sumpars$sb)
sumpars$ntheta <- length(sumpars$sb)+1
sumpars$loc <- rep(0, sumpars$k)
sumpars$t2 <- 2:sumpars$ntheta
tr <- pars2simmap(sumpars, tree)
pdf("../output/figures/koko1.pdf", height=50, width=8)
plotRegimes(tr$tree, cex=0.5)
dev.off()

summarizeDerivedState <- function(branch, chain){
  if(branch==0){
    Th <- sapply(chain$theta, function(x) x[1])
    B1 <- sapply(chain$beta1, function(x) x[1])
    #B2 <- sapply(chain$beta2, function(x) x[1])
    #B3 <- unlist(chain$beta3)
  } else {
    SB <- unlist(chain$sb)
    gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
    ind <- which(SB==branch)
    gen <- gen[ind]
    T2 <- unlist(chain$t2)[ind]
    B1 <- sapply(1:length(T2), function(x) chain$beta1[[gen[x]]][T2[x]])
    Th <- sapply(1:length(T2), function(x) chain$theta[[gen[x]]][T2[x]])
    #B2 <- unlist(chain$beta2)[gen]
    #B3 <- unlist(chain$beta3)[gen]
  }
  medians = list(theta=median(Th), beta1=median(B1))#, beta2=median(B2))
  densities = list(theta=density(Th), beta1=density(B1))#, beta2=density(B2))
  return(list(medians=medians, densities=densities))
}

sb <- sumpars$sb
cladesummaries <- lapply(c(0, sb), function(x) summarizeDerivedState(x, chain))
regressions <- t(sapply(cladesummaries, function(x) unlist(x$medians)))
sumpars$endo <- median(chain$endo[floor(0.3*length(chain$endo)):length(chain$endo)])
rownames(regressions) = c("root", sumpars$sb[sb])
#cache <- bayou:::.prepare.ou.univariate(tree, cache$dat, SE=0, pred)
tipregs <- bayou:::.tipregime(sumpars, tree)
descendents <- lapply(1:(length(sb)+1), function(x) names(tipregs[tipregs==x])) 
#descendents <- lapply(sb, function(x) na.exclude(cache$tip.label[cache$edge[c(cache$bdesc[[x]], x), 2]] ))
nodesc <- which(sapply(descendents, length)==0)


pdf(paste("../output/figures/regressions_", args[[1]],"_",args[[3]],".pdf"), height=8, width=12)
par(mfrow=c(2,2), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
c1 <- "pink"
palx <- rainbow
#beta3 <- median(chains$beta3)
#missing.pred <- do.call(rbind, chainsmp$missing.pred)
#mps <- apply(missing.pred, 2, median)
#impPred <- pred
#impPred[is.na(pred[,3]),3] <- mps
dat <- cache$dat
if(length(nodesc)>0){
  if(length(nodesc)==1 & nodesc==1){
    regressions <- regressions[-nodesc, ]
    descendents <- descendents[-nodesc]
    cladesummaries <- cladesummaries[-nodesc]
  }
  else {
  regressions <- regressions[-nodesc, ]
  descendents <- descendents[-nodesc]
  cladesummaries <- cladesummaries[-nodesc]
  sumpars$sb <- sumpars$sb[-(nodesc-1)]
  sumpars$t2 <- sumpars$t2[-(nodesc-1)]
  sumpars$loc <- sumpars$loc[-(nodesc-1)]
  }
}
for(i in (1:nrow(regressions))){
  plotBayoupars(sumpars, tree, col=setNames(c(palx(nrow(regressions))[i], rep("gray80", nrow(regressions)-1)), c(i, (1:nrow(regressions))[-i])), cex=0.2)
  plot(pred[,1], dat, xlab="lnMass", ylab="lnBMR", pch=21, bg=makeTransparent("gray20", 100), col =makeTransparent("gray20", 100) )
  include <- which(names(dat) %in% descendents[[i]])
  text(pred[include, 1], dat[include], labels=names(dat[include]), col="white", cex=0.4, pos = 2)
  points(pred[include,1], dat[include], pch=21, bg=palx(nrow(regressions))[i])
  print(descendents[[i]])
  #expected <- regressions[i,1]+regressions[i,2]*pred[include,1]+sumpars$endo*pred[include,3]#+beta3*impPred[include,3]
  #o <- order(pred[include,1])
  #lines(pred[include,1][o], expected[o], col=palx(nrow(regressions))[i], lwd=2)
  #plot(cladesummaries[[i]]$densities$beta2, col=palx(nrow(regressions))[i], xlim=c(-0.05, 0.05), lwd=3, main="Beta2")
  #abline(v=0, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(cladesummaries[[i]]$densities$beta1, col=palx(nrow(regressions))[i], xlim=c(0.5, 1.25), lwd=3, main="Beta1")
  abline(v=0.75, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(cladesummaries[[i]]$densities$theta, col=palx(nrow(regressions))[i], xlim=c(-6,3), lwd=3, main="Theta")
  
}
dev.off()
