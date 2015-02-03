## Load in command args:
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
print(args)
#args <- list("")

## Load tree and data into R
require(ape)
require(aRbor)
require(devtools)
#load_all("~/repos/bayou/bayou_1.0")
#install_github("uyedaj/bayou", ref="dev")
require(bayou)

td <- readRDS(paste("../output/data/", args[[1]], ".rds", sep=""))
tree <- td$phy
dat <- td$dat
tree <- multi2di(tree)
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps
td <- make.treedata(tree, dat)
td <- reorder(td, "postorder")
td <- mutate(td, lnMass = log(mean.mass), lnBMR = log(q10smr))
td <- filter(td, !is.na(lnMass), !is.na(lnBMR))
tree <- td$phy
dat <- td$dat

#par(mfrow=c(1,2))
#plot(tree, show.tip.label=FALSE)
#plot(dat$lnMass, dat$lnBMR)
lnBMR <- setNames(dat[['lnBMR']], tree$tip.label)
lnMass <- setNames(dat[['lnMass']], tree$tip.label)
pred <- cbind(setNames(dat[['lnMass']], tree$tip.label),setNames(dat[,'lnMass']^2, tree$tip.label))

cache <- bayou:::.prepare.ou.univariate(tree, setNames(dat$lnBMR, tree$tip.label), pred = pred)



## We're going to define a custom likelihood function; this is nearly identical to bayou.lik (bayou's standard
## OU function), except we use pars$beta1 to calculate the residuals after accounting for lnMass.
custom.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  ##THiS IS INEFFICIENT AND HORRIBLE, REPLACE
  nodedesc <- cache$edge[,2]
  tipedge <- which(cache$edge[,2] <= cache$ntips)
  tipnode <- as.character(nodedesc[tipedge])
  tipreg <- bayou:::.pars2map(pars,cache)$theta
  betaID <- sapply(1:length(tipnode), function(x) tail(tipreg[tipnode[x]],1))
  betaID <- betaID[order(as.numeric(names(betaID)))]
  X = X - sapply(1:cache$ntips, function(x) sum(pars$beta1[betaID[x]]*cache$pred[x,1]+pars$beta2[betaID[x]]*cache$pred[x,2]))
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
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

## This is a function to monitor the output of bayou for our custom model
BetaBMR.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1","rbeta2", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$beta2[2], pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.BetaBMR <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", beta1=".vectorMultiplier", beta2=".vectorMultiplier", k=".splitmergebd", "theta"=".adjustTheta", slide=".slide"),
                      control.weights = list("alpha"=4,"sig2"=2,"beta1"=4, "beta2"=4, "theta"=4,"slide"=2,"k"=10),
                      D = list(alpha=1, sig2= 1, beta1=0.05, beta2=0.05, k = c(1, 1, 1), theta=1, slide=1),
                      parorder = c("alpha", "sig2", "beta1", "beta2", "k", "ntheta", "theta"),
                      rjpars = c("beta1","beta2", "theta"),
                      shiftpars = c("sb", "loc", "t2"),
                      monitor.fn = BetaBMR.monitor,
                      lik.fn = custom.lik)

## Now we define the prior:
prior <- make.prior(tree, plot.prior = FALSE, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",
                                     dbeta1="dnorm",dbeta2="dnorm", dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                          param=list(dalpha=list(scale=1), dsig2=list(scale=1), 
                               dbeta1=list(mean=0.7, sd=0.1), dbeta2=list(mean=0, sd=0.01), dk=list(lambda=50, kmax=500), dsb=list(bmax=1,prob=1), 
                               dtheta=list(mean=0, sd=2.5)), model="ffancova")


startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(31, 0.7, 0.1), beta2=rnorm(31,0,0.01), k=30, ntheta=31, theta=rnorm(31, 0, 1), sb=c(120, 390, 389, 395, 676, 675, 841, 1696, 1695, 1701, 1703, 1702, 1704, 1708, 119, 578, 262, 840, 1258, 1392, 1331, 342, 196, 1254, 1500, 380, 113, 788, 575, 378), loc=rep(0,30), t2=2:31)
mymcmc <- bayou.makeMCMC(tree, lnBMR, pred=pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=100, samp = 50)

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

