#devtools::install_github("uyedaj/bayou", ref="dev")
require(bayou)
require(devtools)
require(aRbor)
require(Matrix)
source("./betaBayouFunctions.R")
args <- list("tetrapods_ei", "100000", "_RR000_ME2", "sumpars_u6", "RR000")
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
source("./tetrapods_ei_Models.R")

prior <- priors[paste("prior", args[[5]], sep=".")][[1]]
model <- models[paste("model", args[[5]], sep=".")][[1]]
startpar <- startpars[paste("model", args[[5]], sep=".")][[1]]
#startpar$theta <- rep(tmp$coef[1], length(startpar$theta))
#startpar$beta_lnMass <- rep(tmp$coef[2], length(startpar$beta_lnMass))
#startpar$beta_TempK <- rep(tmp$coef[4])
#startpar$beta_lnMass2 <- tmp$coef[3]
#startpar$alpha <- 0.06
#prior(startpar)
#model$lik.fn(startpar, cache, cache$dat)$loglik
SEs <- c(0.05, 0.1, 0.3)

mymcmc1 <- bayou.makeMCMC(cache$phy, cache$dat, pred=cache$pred, SE=SEs[1], model=model, prior=prior, startpar=startpar, new.dir=paste("../output/runs/n",args[[1]],sep=""), outname=paste(args[[1]],args[[3]], "ME1",sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)
mymcmc1$run(args[[2]])

mymcmc2 <- bayou.makeMCMC(cache$phy, cache$dat, pred=cache$pred, SE=SEs[2], model=model, prior=prior, startpar=startpar, new.dir=paste("../output/runs/n",args[[1]],sep=""), outname=paste(args[[1]],args[[3]], "ME2",sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)
mymcmc2$run(args[[2]])

mymcmc3 <- bayou.makeMCMC(cache$phy, cache$dat, pred=cache$pred, SE=SEs[3], model=model, prior=prior, startpar=startpar, new.dir=paste("../output/runs/n",args[[1]],sep=""), outname=paste(args[[1]],args[[3]], "ME3",sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)
mymcmc3$run(args[[2]])

pdf("../output/mcmcwSE.pdf", height=12, width=8)
require(coda)
chain1 <- mymcmc1$load()
chain1 <- set.burnin(chain1, 0.3)
chain2 <- mymcmc2$load()
chain2 <- set.burnin(chain2, 0.3)
chain3 <- mymcmc3$load()
chain3 <- set.burnin(chain3, 0.3)
plot(chain1)
plot(chain2)
plot(chain3)
dev.off()

chain <- mymcmc$load()
chain <- set.burnin(chain, 0.4)
saveRDS(chain1, file=paste("../output/runs/n",args[[1]], "/", args[[1]], args[[3]],"_chain_ME1.rds",sep=""))
saveRDS(mymcmc1, file=paste("../output/runs/n",args[[1]], "/", args[[1]], args[[3]],"_mcmc_ME1.rds",sep=""))

saveRDS(chain2, file=paste("../output/runs/n",args[[1]], "/", args[[1]], args[[3]],"_chain_ME2.rds",sep=""))
saveRDS(mymcmc2, file=paste("../output/runs/n",args[[1]], "/", args[[1]], args[[3]],"_mcmc_ME2.rds",sep=""))

saveRDS(chain3, file=paste("../output/runs/n",args[[1]], "/", args[[1]], args[[3]],"_chain_ME3.rds",sep=""))
saveRDS(mymcmc3, file=paste("../output/runs/n",args[[1]], "/", args[[1]], args[[3]],"_mcmc_ME3.rds",sep=""))

require(foreach)
require(doParallel)
registerDoParallel(cores=5)
Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)
ss <- mymcmc$steppingstone(args[[2]], chain, Bk, burnin=0.3, plot=TRUE)
plot(ss)
print(ss$lnr)
saveRDS(ss, file=paste("../output/runs/n",args[[1]],"/",args[[1]],"_",args[[3]],".ss.rds",sep=""))

#plot(chain)
#sumstats <- summary(chain)


#cutoff <- 0.1
#sumpars <- list(sb = which(sumstats$branch.posteriors$pp > cutoff))
#sumpars$k <- length(sumpars$sb)
#sumpars$ntheta <- length(sumpars$sb)+1
#sumpars$loc <- rep(0, sumpars$k)
#sumpars$t2 <- 2:sumpars$ntheta
#sumpars$betaC <- sumstats$statistics['betaC', 1]
#sumpars$betaH <- median(chain$betaH[200:500])
#sumpars$betaT <- median(chain$betaT[200:500])
#sumpars$theta <- sumstats$statistics['root.theta',1]
#sumpars$beta1 <- sumstats$statistics['root.beta1',1]
#tr <- pars2simmap(sumpars, tree)

#plotSimmap(tr$tree, colors=tr$col, fsize=0.25)
#summarizeDerivedState <- function(branch, chain){
#  if(branch==0){
#    Th <- sapply(chain$theta, function(x) x[1])
#    B1 <- sapply(chain$beta1, function(x) x[1])
#    B2 <- sapply(chain$beta2, function(x) x[1])
#    #B3 <- unlist(chain$beta3)
#  } else {
#    SB <- unlist(chain$sb)
#    gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
#    ind <- which(SB==branch)
#    gen <- gen[ind]
#    T2 <- unlist(chain$t2)[ind]
#    B1 <- sapply(1:length(T2), function(x) chain$beta1[[gen[x]]][T2[x]])
#    Th <- sapply(1:length(T2), function(x) chain$theta[[gen[x]]][T2[x]])
#    B2 <- unlist(chain$beta2)[gen]
#    #B3 <- unlist(chain$beta3)[gen]
#  }
 # medians = list(theta=median(Th), beta1=median(B1), beta2=median(B2))
#  densities = list(theta=density(Th), beta1=density(B1), beta2=density(B2))
#  return(list(medians=medians, densities=densities))
#}


#tipregs <- bayou:::.tipregime(sumpars, tree)
#descendents <- lapply(1:(length(sb)+1), function(x) names(tipregs[tipregs==x])) 

#par(mfrow=c(1,2))
#plot(cache$pred$log.no, cache$dat, pch=21, bg=cache$pred$growthForm, col=cache$pred$growthForm, main="Regressions", ylab="log Size")
#curve(sumpars$theta+sumpars$betaC+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col=1)
#curve(sumpars$theta+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col=3)
#curve(sumpars$theta+sumpars$betaH+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col=2)
#curve(sumpars$theta+sumpars$betaT+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col="blue")


#plot(density(sapply(chain$theta[200:500], function(x) x[1])), col="green", lwd=3, xlim=c(-6, 0), main="Densities", xlab="Intercept")
#lines(density(chain$betaH[200:500]+sapply(chain$theta[200:500], function(x) x[1])), col="red", lwd=3)
#lines(density(chain$betaC[200:500]+sapply(chain$theta[200:500], function(x) x[1])), col="black", lwd=3)
#lines(density(chain$betaT[200:500]+sapply(chain$theta[200:500], function(x) x[1])), col="blue", lwd=3)

