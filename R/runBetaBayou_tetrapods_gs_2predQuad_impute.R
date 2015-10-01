## Load in command args:
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
print(args)
rm(list=ls(all=TRUE))
setwd("~/repos/bmr/R/")
args <- list("tetrapods_gs", 10000, "impute_r1")

## Load tree and data into R
require(ape)
require(aRbor)
require(devtools)
require(mvnfast)
#load_all("~/repos/bayou/bayou_1.0")
#install_github("uyedaj/bayou", ref="dev")
load_all("~/repos/bayou/bayou_1.0")

require(bayou)

td <- readRDS(paste("../output/data/", args[[1]], ".rds", sep=""))
tree <- td$phy
dat <- td$dat
tree <- multi2di(tree)
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps
td <- make.treedata(tree, dat)
td <- reorder(td, "postorder")
td_gs <- filter(td, !is.na(lnGenSize))
#td <- mutate(td, lnMass = log(mean.mass), lnBMR = log(q10smr))
#td <- filter(td, !is.na(lnMass), !is.na(lnBMR), !is.na(lnGenSize))
tree <- td$phy
dat <- td$dat
## BM likelihood function for genome size
gs.lik <- bm.lik(td_gs$phy, setNames(td_gs$dat[[3]], td_gs$phy$tip.label), SE=0, model="BM")


#par(mfrow=c(1,2))
#plot(tree, show.tip.label=FALSE)
#plot(dat$lnMass, dat$lnBMR)
lnBMR <- setNames(dat[['lnBMR']], tree$tip.label)
lnMass <- setNames(dat[['lnMass']], tree$tip.label)
pred <- cbind(setNames(dat[['lnMass']], tree$tip.label),setNames(dat[['lnMass']]^2, tree$tip.label),setNames(dat[['lnGenSize']], tree$tip.label))

cache <- bayou:::.prepare.ou.univariate(tree, setNames(dat$lnBMR, tree$tip.label), pred = pred)

missing <- which(is.na(cache$pred[,3]))

getPreValues <- function(cache){
  V <- vcvPhylo(cache$phy, anc.nodes=FALSE)
  X <- cache$pred[,3]
  unknown <- is.na(X)
  known <- !unknown
  Vkk <- V[known, known]
  Vuu <- V[unknown, unknown]
  Vku <- V[known, unknown]
  Vuk <- V[unknown, known]
  iVkk <- solve(Vkk)
  sigmabar <- as.matrix(forceSymmetric(Vuu - Vuk%*%iVkk%*%Vku))
  cholSigmabar <- chol(sigmabar)
  mubarmat <- Vuk%*%iVkk
  return(list(V=V, X=X, unknown=unknown, known=known, Vkk=Vkk, Vuu=Vuu, Vku=Vku, Vuk=Vuk, iVkk=iVkk, sigmabar=sigmabar, mubarmat=mubarmat, cholSigmabar=cholSigmabar))
}
pv <- getPreValues(cache)

cMVNorm <- function(cache, pars, prevalues=pv, known=FALSE){
  X <- prevalues$X
  known <- prevalues$known
  unknown <- prevalues$unknown
  mu <- rep(pars$pred.root, cache$n)
  muk <- mu[known]
  muu <- mu[unknown]
  mubar <- t(muu + prevalues$mubarmat%*%(X[known]-muk))
  #sigmabar <- pars$pred.sig2*prevalues$sigmabar
  myChol <-sqrt(pars$pred.sig2)*prevalues$cholSigmabar
  res <- dmvn(pars$missing.pred, mu=mubar, sigma = myChol, log=TRUE, isChol=TRUE)
  return(res)
}

## Proposal function to simulate conditional draws from a multivariate normal distribution
.imputePredBM <- function(cache, pars, d, move,ct=NULL, prevalues=pv){
  #(tree, dat, sig2, plot=TRUE, ...){
  X <- prevalues$X
  Vuk <- pars$pred.sig2*prevalues$Vuk
  iVkk <- (1/pars$pred.sig2)*prevalues$iVkk
  Vku <- pars$pred.sig2*prevalues$Vku
  Vuu <- pars$pred.sig2*prevalues$Vuu
  known <- prevalues$known
  unknown <- prevalues$unknown
  mu <- rep(pars$pred.root, cache$n)
  muk <- mu[known]
  muu <- mu[unknown]
  mubar <- t(muu + Vuk%*%iVkk%*%(X[known]-muk))
  sigmabar <- Vuu - Vuk%*%iVkk%*%Vku
  res <- MASS::mvrnorm(1, mubar, sigmabar)
  pars.new <- pars
  pars.new$missing.pred <- res
  hr=Inf
  type="impute"
  return(list(pars=pars.new, hr=hr, decision = type))
}

## We're going to define a custom likelihood function; this is nearly identical to bayou.lik (bayou's standard
## OU function), except we use pars$beta1 to calculate the residuals after accounting for lnMass.
custom.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  pred[is.na(pred[,3]),3] <- pars$missing.pred
  tipreg <- rev(bayou:::.pars2map(pars,cache)$theta)
  ntipreg <- names(tipreg)
  dups <- !duplicated(ntipreg) & ntipreg %in% (1:nrow(cache$edge))[cache$externalEdge]
  tipreg <- tipreg[dups]
  ntipreg <- ntipreg[dups]
  o <- order(cache$edge[as.numeric(ntipreg), 2])
  betaID <- tipreg[o]
  #betaID <- sapply(tipreglist[cache$externalEdge][o], function(x) x[length(x)])
  #beta <- cbind(sapply(c("beta1", "beta2"), function(x) pars[[x]][betaID]), pars$beta3)
  X = X - pars$beta1[betaID]*pred[,1] + pars$beta2[betaID]*pred[,2]+pars$beta3*pred[,3]
  cache$dat <- X
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
  llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN)
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}


plotSimmap(pars2simmap(pars, cache$phy)$tree, colors=pars2simmap(pars, cache$phy)$col, ftype="off")

startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(31, 0.7, 0.1), beta2=rnorm(31,0,0.01), beta3=rnorm(1,0,0.1), k=30, ntheta=31, theta=rnorm(31, 0, 1), sb=c(120, 390, 389, 395, 676, 675, 841, 1696, 1695, 1701, 1703, 1702, 1704, 1708, 119, 578, 262, 840, 1258, 1392, 1331, 342, 196, 1254, 1500, 380, 113, 788, 575, 378), loc=rep(0,30), t2=2:31)

## Get optimized starting values and trait evolutionary models for genome size. 
require(phylolm)
tdgs <- filter(td, !is.na(lnGenSize))
fits <- lapply(c("BM", "OUrandomRoot", "OUfixedRoot", "EB"), function(x) phylolm(lnGenSize~1, data=tdgs$dat, phy=tdgs$phy, model=x))
aics <- sapply(fits, function(x) x$aic)
bestfit <- fits[[which(aics == min(aics))]]
#phenogram(tdgs$phy, setNames(tdgs$dat$lnGenSize, tdgs$phy$tip.label), fsize=0.5)

startpar$pred.root <- unname(bestfit$coeff)
startpar$pred.sig2 <- unname(bestfit$sigma2)
#phenogram(td$phy, setNames(tmppred3[[1]], td$phy$tip.label), ftype="off", colors = makeTransparent("black", 0))
startpar <- .imputePredBM(cache, startpar, d=1, NULL, ct=NULL, prevalues=pv)$pars
tmppred3 <- as.vector(td$dat[,3])
tmppred3[is.na(tmppred3), ] <- startpar$missing.pred
#phenogram(td$phy, setNames(tmppred3[[1]], td$phy$tip.label), ftype="off", colors = makeTransparent("black", 5), add=TRUE)

## This is a function to monitor the output of bayou for our custom model
BetaBMR.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1","rbeta2","rbeta3", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$beta2[2], pars$beta3, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.BetaBMR <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", beta1=".vectorMultiplier", beta2=".vectorSlidingWindow",beta3=".vectorSlidingWindow", 
                                   k=".splitmergebd", theta=".adjustTheta", slide=".slide",  pred.sig2=".multiplierProposal", 
                                   pred.root=".slidingWindowProposal", missing.pred=".imputePredBM"),
                      control.weights = list("alpha"=4,"sig2"=2,"beta1"=4, "beta2"=4, "beta3"=4, "theta"=4,"slide"=2,"k"=10, "pred.root"=2, "pred.sig2"=2, "missing.pred"=3),
                      D = list(alpha=1, sig2= 0.75, beta1=0.75, beta2=0.1, beta3=0.2, k = c(1, 1, 1), theta=2, slide=1, pred.sig2=1, pred.root=2, missing.pred=1),
                      parorder = c("alpha", "sig2", "beta1", "beta2", "beta3","pred.sig2", "pred.root", "missing.pred", "k", "ntheta", "theta"),
                      rjpars = c("beta1","beta2", "theta"),
                      shiftpars = c("sb", "loc", "t2"),
                      monitor.fn = BetaBMR.monitor,
                      lik.fn = custom.lik)

## Now we define the prior:
prior <- make.prior(tree, plot.prior = FALSE, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",
                                     dbeta1="dnorm",dbeta2="dnorm",dbeta3="dnorm", dsb="dsb", dk="cdpois", dtheta="dnorm",
                                     dpred.sig2="dhalfcauchy", dpred.root="dnorm"), 
                          param=list(dalpha=list(scale=1), dsig2=list(scale=1), 
                               dbeta1=list(mean=0.7, sd=0.1), dbeta2=list(mean=0, sd=0.01),dbeta3=list(mean=0, sd=0.5), dk=list(lambda=50, kmax=500), dsb=list(bmax=1,prob=1), 
                               dtheta=list(mean=0, sd=2.5), dpred.sig2=list(scale=1), dpred.root=list(mean=1, sd=1)), model="ffancova")

prior(startpar)

mymcmc <- bayou.makeMCMC(tree, lnBMR, pred=pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)

## Now run the chain for 50,000 generations. This will write to a set of files in your temporary directory.

#sink(paste(mymcmc$dir,mymcmc$outname,".log",sep=""), append=TRUE)

mymcmc$run(args[[2]])

#sink(NULL)

chain <- mymcmc$load()
chain <- set.burnin(chain, 0.3)
out <- summary(chain)
print(out)

saveRDS(chain, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".chain.rds",sep=""))
saveRDS(mymcmc, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_",args[[3]],".mcmc.rds",sep=""))
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

