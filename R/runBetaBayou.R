## Load in command args:
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
print(args)

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
cache <- bayou:::.prepare.ou.univariate(tree, setNames(dat$lnBMR, tree$tip.label))

lnBMR <- setNames(dat[,'lnBMR'], tree$tip.label)
lnMass <- setNames(dat[,'lnMass'], tree$tip.label)
pred <- cbind(setNames(dat[,'lnMass'], tree$tip.label))


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
  X = X - sapply(1:cache$ntips, function(x) sum(pars$beta1[betaID[x]]*cache$pred[x]))
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
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.BetaBMR <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", beta1=".vectorMultiplier", k=".splitmergebd", "theta"=".adjustTheta", slide=".slide"),
                      control.weights = list("alpha"=4,"sig2"=2,"beta1"=4, "theta"=4,"slide"=2,"k"=10),
                      D = list(alpha=2, sig2= 1, beta1=0.05, k = c(1, 1), theta=1, slide=1),
                      parorder = c("alpha", "sig2", "beta1", "k", "ntheta", "theta"),
                      rjpars = c("beta1", "theta"),
                      shiftpars = c("sb", "loc", "t2"),
                      monitor.fn = BetaBMR.monitor,
                      lik.fn = custom.lik)

## Now we define the prior:
prior <- make.prior(tree, plot.prior = FALSE, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",
                                     dbeta1="dnorm", dsb="dsb", dk="cdpois", dtheta="dunif"), 
                          param=list(dalpha=list(scale=1), dsig2=list(scale=1), 
                               dbeta1=list(mean=0.7, sd=0.1), dk=list(lambda=15, kmax=200), dsb=list(bmax=1,prob=1), 
                               dtheta=list(min=-10, max=100)), model="ffancova")


startpar <- list(alpha=0.1, sig2=3, beta1=c(0.66, 0.75), k=1, ntheta=2, theta=c(0,0), sb=sample(1:nrow(dat), 1), loc=0, t2=2)
mymcmc <- bayou.makeMCMC(tree, lnBMR, pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=args[[1]], plot.freq=NULL, ticker.freq=10000, samp = 50)

## Now run the chain for 50,000 generations. This will write to a set of files in your temporary directory.

#sink(paste(mymcmc$dir,mymcmc$outname,".log",sep=""), append=TRUE)
mymcmc$run(args[[2]])
#sink(NULL)

chain <- mymcmc$load()
chain <- set.burnin(chain, 0.3)
out <- summary(chain)

saveRDS(chain, file=paste("../output/runs/",args[[1]],"/",args[[1]],"chain.rds",sep=""))
saveRDS(mymcmc, file=paste("../output/runs/",args[[1]],"/",args[[1]],"mcmc.rds",sep=""))
#plotSimmap.mcmc(tree, chain, burnin=0.3)
