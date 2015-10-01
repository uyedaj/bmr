require(ape)
require(geiger)
require(phytools)
pars <- list(beta1=0.75, beta1within=0.82, beta0=-1.2, sig2=1, alpha=0.5, sig2within=0.4, E=0.65)
n <- 1000
k <- 8.617332478E-5

0.62/(k*(38+273.15))

samples <- rpois(n, 5)+1

tree <- sim.bdtree(1, d=0, stop="taxa", n=n)
lnMass <- fastBM(tree, 4, 1)
lnMassind <- lapply(1:n, function(x) rnorm(samples[x], lnMass[x], 1))
tempind <- lapply(1:n, function(x) runif(samples[x], 10,40)+274.15)
MRind <- lapply(1:n, function(x) rnorm(samples[x], pars$beta0+pars$beta1within*lnMassind[[x]]-pars$E*(1/(k*tempind[[x]])), sd=sqrt(pars$sig2within)))

resid <- fastBM(tree, 0, pars$sig2)
resid+pars$beta0+pars$beta1


plot(Temp1, E*(Temp1-Temp2)/(k*Temp2*Temp1) + log(Temp2))
#plot(Temp1, E*(Temp2-Temp1)/(k*T1) + log(Temp1))
lines(0:40, E*(0:40+273.15-Temp2)/(k*Temp2*(0:40+273.15)) + log(Temp2)-5, col="red", lwd=2)


i<-27
lnMassind[[i]]
plot(lnMassind[[i]], MRind[[i]],cex=(tempind[[x]]-290)/10)

T1 <- 320
T0 <- 279
T2 <- 293.15

(pars$E/k)*(1/T1- 1/T2)
(pars$E/k)*(1/T0- 1/T2)


MTstand <- function(lnBMRs, lnMasses, TempKs, E=0.65, T2=293.15){
  stds <- (lnBMRs + E/k*(1/TempKs-1/T2))/lnMasses*mean(lnMasses)
  return(list(x=mean(stds), raw=stds, original=lnBMRs))
}

MTstand2 <- function(lnBMRs, lnMasses, TempKs, E=0.65, T2=293.15){
  stds <- log((exp(lnBMRs)*exp(E/k*(1/TempKs-1/T2)))/exp(lnMasses)*exp(mean(lnMasses)))
  return(list(x=mean(stds), raw=stds, original=lnBMRs))
}

MTstand3 <- function(lnBMRs, lnMasses, TempKs, E=0.65, T2=293.15){
  stds <-lnBMRs + E/k*(1/TempKs-1/T2)+pars$beta1within*(mean(lnMasses)- lnMasses)
  return(list(x=mean(stds), raw=stds, original=lnBMRs))
}


MTstand(MRind[[1]], lnMassind[[1]], tempind[[1]])

pars$beta1within <-0.75
tBs <- 2
tM <- 4
N <- 100000
lnMasses <- c(rnorm(N/2, tM, sd=1), rnorm(N/2, tM+2, sd=1))
temps <- runif(N, 293.15+1, 293.15+10)
vals <- rnorm(N, tBs, sd=0.2)
etBs <- mean(vals)
eTvals <- vals - pars$E/k*(1/temps - 1/293.15)
#par(mfrow=c(2,2))
#plot(rep(293.15, N), vals, col=1:N, ylim=c(0, 3))
#plot(temps, eTvals, col=1:N, ylim=c(0, 3))
eMTvals <- eTvals + pars$beta1within*(lnMasses - mean(lnMasses))
#plot(lnMasses, eMTvals, col=1:N)

MTstand(eMTvals, lnMasses, temps, E=pars$E)$x
MTstand2(eMTvals, lnMasses, temps, E=pars$E)$x
MTstand3(eMTvals, lnMasses, temps, E=pars$E)$x
etBs

var(MTstand3(eMTvals, lnMasses, temps, E=pars$E)$raw)
var(vals)
(var(eMTvals) + pars$beta1within^2*var(mean(lnMasses)-lnMasses)+pars$E/k*var(1/temps))

## DO THE FULL MULTIVARIATE NORMAL SHIITEEEEE
## First let's simulate data under the assumed process:
require(bayou)
require(dplyr)
k <- 8.617332478E-5; T2 <- 293.15
alpha <- log(2)/0.12
sig2 <- 0.1
ntaxa <- 100
nshifts <- 5
theta <- rnorm(nshifts+1, -1, 1)
betaR <- runif(nshifts+1, 0.6, 0.8)
eR <- runif(nshifts+1, 0.58, 0.72)
nwithin <- 1+rnbinom(ntaxa, 1, prob=0.1)
hist(nwithin)
betaij <- runif(ntaxa, 0.2, 1)
tree <- sim.bdtree(1.0, d=0, stop="taxa", n=ntaxa, seed=1) %>% reorder("postorder")
tree$edge.length <- tree$edge.length/max(tree$edge.length)
plot(tree)
if(!(exists("pars"))){
pars <- identifyBranches(tree, nshifts)
}
pars$alpha <- alpha
pars$sig2 <- sig2
pars$t2 <- 2:(nshifts+1)
pars$theta <- theta
pars$betaR <- betaR
pars$k <- nshifts
pars$betaij <- betaij
pars$ntheta <- nshifts+1
pars$eR <- eR
sig2w <- 0.17
Mj <- fastBM(tree, sig2=3)+5
phenogram(tree, Mj)
Mij <- lapply(1:ntaxa, function(x) log(rnorm(nwithin[x], mean=exp(Mj[x]), sd=0.25*exp(Mj[x]))))
eMj <- sapply(Mij, mean)
plot(eMj, Mj)
dat0 <- dataSim(pars, model="OU", tree)
dat1 <-  dat0$dat + pars$betaR[bayou:::.tipregime(pars, tree)]*eMj

edat0 <- lapply(1:ntaxa, function(x) rep(dat0$E.th[x], nwithin[x]))
edat1 <- lapply(1:ntaxa, function(x) rep(dat1[x], nwithin[x]))
Tij <- lapply(1:ntaxa, function(x) runif(nwithin[x], 15, 25)+273.15)
dat <- lapply(1:ntaxa, function(x) sapply(1:nwithin[x], function(y)
    edat1[[x]][y] - pars$betaij[x]*(eMj[x]-Mij[[x]][y]) - eR[bayou:::.tipregime(pars, tree)][x]/k*(1/(Tij[[x]][y])- 1/293.15)
  ))
plot(Tij[[1]], dat[[1]])
plot(Mij[[1]], dat[[1]])

errw <- lapply(1:ntaxa, function(x) rnorm(nwithin[x], 0, sqrt(sig2w)))
Bij <- lapply(1:ntaxa, function(x) dat[[x]] + errw[[x]]) 

plot(eMj, dat1, pch=21, bg=bayou:::.tipregime(pars,tree))
lapply(1:pars$ntheta, function(x) abline(pars$theta[x], pars$betaR[x], lty=2))
tipregime <- bayou:::.tipregime(pars, tree)

plot(eMj, sapply(Bij, mean), type="n", xlim=c(0,max(unlist(Mij))+1), ylim=c(min(unlist(Bij))-1, max(unlist(Bij))+1))
eBij <- sapply(Bij, mean)
lapply(1:pars$ntheta, function(x) abline(pars$theta[x], pars$betaR[x], lty=2))
lapply(1:ntaxa, function(x) lapply(1:nwithin[x], function(y) lines(c(eMj[x], Mij[[x]][y]), c(eBij[x], Bij[[x]][y]), pch=".", col=tipregime[x])))
points(eMj, eBij,  pch=21, bg=tipregime)

## Data has now been simulated. Now let's check to see if we can calculate the correct likelihood.
##  The expectation of each individual data point is then.
## Let's figure out how to make the Regime matrix first...

expand <- function(x) {
  unlist(lapply(1:ntaxa, function(y) rep(x[y], nwithin[y])), F,F)
}
tipregimeij <- expand(tipregime)
species <- unlist(lapply(1:ntaxa, function(x) rep(x, nwithin[x])))
R <- matrix(0, nrow=length(species), ncol=pars$ntheta)
vtipregime <- sapply(1:length(species), function(x) x + (tipregimeij[x]-1)*length(tipregimeij))
R[vtipregime] <- 1

EBj <- expand(parmap.W(tree, pars)%*%pars$theta) +  R %*% pars$betaR * expand(eMj) 
EBjm <- EBj -  expand(pars$betaij)* (expand(eMj) - unlist(Mij))
EBjmt <- EBjm - R%*%pars$eR/k *(1/unlist(Tij) - 1/T2)
par(mfrow=c(2,2))
var(unlist(Bij) - EBj)
var(unlist(Bij) - EBjm)
var(unlist(Bij) - EBjmt)
plot(density((unlist(Bij) - EBj)))
lines(density(unlist(Bij)- EBjm))


ExpFn <- function(pars, eMj, Mij, Tij) {
  expand(parmap.W(tree, pars)%*%pars$theta) +  R %*% pars$betaR * expand(eMj) -  expand(pars$betaij)* (expand(eMj) - unlist(Mij)) -  R%*%pars$eR/k *(1/unlist(Tij) - 1/T2)
}


resid <- (unlist(Bij) - ExpFn(pars, eMj, Mij, Tij))
etree <- tree





