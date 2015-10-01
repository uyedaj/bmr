args <- list("tetrapods", "10000", "r1")
setwd("~/repos/bmr/R")

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
#tree <- multi2di(tree)
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps
td <- make.treedata(tree, dat)
td <- reorder(td, "postorder")
#td <- mutate(td, lnMass = log(mean.mass), lnBMR = log(q10smr))
td <- filter(td, !is.na(lnMass), !is.na(lnBMR))
tree <- td$phy
dat <- td$dat
lnBMR <- setNames(dat[,'lnBMR'], tree$tip.label)
lnMass <- setNames(dat[,'lnMass'], tree$tip.label)
pred <- cbind(setNames(dat[,'lnMass'], tree$tip.label))

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
  X = X - sapply(1:cache$ntips, function(x) sum(pars$beta1[betaID[x]]*cache$pred[x,1]))
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
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.BetaBMR <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", beta1=".vectorMultiplier",  k=".splitmergebd", "theta"=".adjustTheta", slide=".slide"),
                      control.weights = list("alpha"=4,"sig2"=2,"beta1"=4, "theta"=4,"slide"=2,"k"=10),
                      D = list(alpha=1, sig2= 1, beta1=0.05, k = c(1, 1), theta=1, slide=1),
                      parorder = c("alpha", "sig2", "beta1",  "k", "ntheta", "theta"),
                      rjpars = c("beta1","theta"),
                      shiftpars = c("sb", "loc", "t2"),
                      monitor.fn = BetaBMR.monitor,
                      lik.fn = custom.lik)

## Now we define the prior:
prior <- make.prior(tree, plot.prior = FALSE, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",
                                                         dbeta1="dnorm", dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                    param=list(dalpha=list(scale=1), dsig2=list(scale=1), 
                               dbeta1=list(mean=0.7, sd=0.1),  dk=list(lambda=50, kmax=500), dsb=list(bmax=1,prob=1), 
                               dtheta=list(mean=0, sd=2.5)), model="ffancova")


startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(31, 0.7, 0.1), k=30, ntheta=31, theta=rnorm(31, 0, 1), sb=c(120, 390, 389, 395, 676, 675, 841, 1696, 1695, 1701, 1703, 1702, 1704, 1708, 119, 578, 262, 840, 1258, 1392, 1331, 342, 196, 1254, 1500, 380, 113, 788, 575, 378), loc=rep(0,30), t2=2:31)


require(ggplot2)
plotPosteriors <- function(tree, dat, pred, chain, burnin=0.3, cutoff=0.1, thin=1){
  getBranchPosterior <- function(sb, chain, pars=c("beta1", "theta"), burnin=0.3, thin=100){
    postburn <- round(seq(burnin*length(chain$gen), length(chain$gen), length.out=length(chain$gen)/thin),0)
    pbchain <- lapply(chain, function(x) x[postburn])
    if(sb==0){
      out <- lapply(pars, function(x) sapply(pbchain[[x]], function(y) y[1]))
    } else {
      sbs <- unlist(pbchain$sb)
      out <- lapply(pars, function(x) unlist(lapply(pbchain[[x]], function(y) y[-1]),F,F)[sbs==sb])
    }
    return(out)
  }
  LP <- Lposterior(chain, tree, burnin=burnin)
  focalSB <- c(0, which(LP[,1] > cutoff))
  PostSB <- lapply(focalSB, function (y) getBranchPosterior(y, chain, burnin=burnin, thin=thin))
  postBetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[1]]))
  postThetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[2]]))
  postBetas <- data.frame(do.call(rbind, postBetas))
  postBetas$var <- rep("Beta", nrow(postBetas))
  postThetas <- data.frame(do.call(rbind, postThetas))
  postThetas$var <- rep("Theta", nrow(postThetas))
  allpar <- rbind(postBetas, postThetas)
  allpar$x <- factor(allpar$x)
  focalpars <- list(k=length(focalSB[-1]), 
                    ntheta=length(focalSB), 
                    sb=focalSB[-1], 
                    loc=rep(0, length(focalSB[-1])), 
                    t2=2:length(focalSB))
  tr <- pars2simmap(focalpars, tree)
  postburn <- round(seq(burnin*length(chain$gen), length(chain$gen), length.out=length(chain$gen)/thin),0)
  pbchain <- lapply(chain, function(x) x[postburn])
  p <- ggplot(allpar, aes(V2, fill=x))
  p <- p + geom_density(alpha=0.5)
  p <- p + facet_wrap(~ var, scales="free")
  p <- p + scale_fill_manual(values=tr$col)
  #p <- p + ylim(c(0,5))
  #p <- p + xlim(c(0.6, 0.8))
  print(p + geom_vline(x=c(3/4, 2/3), lty=2))
  medianBetas <- tapply(postBetas[,2], postBetas[,1], median)
  medianThetas <- tapply(postThetas[,2], postThetas[,1], median)
  est.pars <- list(beta1=medianBetas, 
                   thetas=medianThetas, 
                   k=length(focalSB)-1, 
                   ntheta=length(focalSB), 
                   sb=focalSB[-1], 
                   loc=rep(0, length(focalSB)-1), 
                   t2=2:length(focalSB))
  tipreg <- bayou:::.tipregime(est.pars, tree)
  plotSimmap.mcmc(tr$tree, pbchain, colors=tr$col, burnin=burnin, fsize=0.1)
  plot(pred[,1], dat, type="n", main=i, 
       xlab="Mass", ylab="BMR", 
       xlim=c(min(pred), max(pred))*c(0.9,1.1), 
       ylim=c(min(dat), max(dat))*c(0.9,1.1))
  for(i in 1:est.pars$ntheta){
    points(pred[tipreg==i,],dat[tipreg==i], pch=21, 
           bg=bayou:::makeTransparent(tr$col[i],50), col=bayou:::makeTransparent(tr$col[i],50))
    abline(a=est.pars$theta[i], b=est.pars$beta1[i], lwd=2, col=tr$col[i],lty=2)
  }
}

#taxa <- list.files("../output/runs")

#for(i in taxa){
#  chain <- readRDS(paste("../output/runs/",i, "/",i,"chain.rds", sep=""))
#  mymcmc <- readRDS(paste("../output/runs/",i, "/",i,"mcmc.rds", sep=""))
#  tree <- mymcmc$tree
#  dat <- mymcmc$dat
#  pred <- mymcmc$pred
#  dir <- paste("../output/figures/", sep="")
#  #dir.create(dir)
#  pdf(paste(dir ,i,"posteriors.pdf", sep=""))
#  try(plotPosteriors(tree, dat, pred, chain, burnin=0.3, cutoff=0.1))
#  dev.off()
#}
truncate <- function(chain){
  mn <- min(sapply(chain, length))-10
  chain <- lapply(chain, function(x) x[1:mn])
  chain
}
args[[3]] <-"r1"
mymcmc1 <- bayou.makeMCMC(tree, lnBMR, pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=10000, samp = 50)
chain1 <- mymcmc1$load()
#chain1 <- truncate(chain1)
args[[3]] <- "r2"
mymcmc2 <- bayou.makeMCMC(tree, lnBMR, pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=10000, samp = 50)
chain2 <- mymcmc2$load()
#chain2 <- truncate(chain2)
args[[3]] <- "r3"
mymcmc3 <- bayou.makeMCMC(tree, lnBMR, pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=10000, samp = 50)
chain3 <- mymcmc3$load()
chain3 <- truncate(chain3)
args[[3]] <- "r4"
mymcmc4 <- bayou.makeMCMC(tree, lnBMR, pred, SE=0, model=model.BetaBMR, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=10000, samp = 50)
chain4 <- mymcmc4$load()
chain4 <- truncate(chain4)
chain1 <- set.burnin(chain1, 0.3)
chain2 <- set.burnin(chain2, 0.3)
chain3 <- set.burnin(chain3, 0.6)
chain4 <- set.burnin(chain4, 0.4)
chaina <- combine.chains(chain1, chain2, burnin.prop=0.3)
chainb <- combine.chains(chain3, chain4, burnin.prop=0.5)
chain <- combine.chains(chaina, chainb, burnin.prop=0)
#chain <- readRDS(paste("../output/runs/",i, "/",i,"chain.rds", sep=""))
#mymcmc <- readRDS(paste("../output/runs/",i, "/",i,"mcmc.rds", sep=""))
#tree <- mymcmc$tree
#dat <- mymcmc$dat
#pred <- mymcmc$pred
dir <- paste("../output/figures/", sep="")
#dir.create(dir)
co <- 0.3
pdf(paste(dir ,i,"posteriors.pdf", sep=""), height=10, width=10)
#try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chain1, burnin=0.3, cutoff=co))
#try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chain2, burnin=0.3, cutoff=co))
#try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chain3, burnin=0.3, cutoff=co))
#try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chain4, burnin=0.3, cutoff=co))
#try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chaina, burnin=0.3, cutoff=co))
#try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chainb, burnin=0.3, cutoff=co))
try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chain, burnin=0.3, cutoff=co))
#try(plotPosteriors(tree, setNames(dat$lnBMR, rownames(dat)), pred, chain2, burnin=0.3, cutoff=0.3))
dev.off()

L1 <- Lposterior(chain1,tree, burnin=0.3)
L2 <- Lposterior(chain2, tree, burnin=0.3)
L3 <- Lposterior(chain3, tree, burnin=0.3)
L4 <- Lposterior(chain4, tree, burnin=0.3)
pp <- cbind(L1$pp, L2$pp, L3$pp, L4$pp)
pairs((pp))
pairs(log(pp))

require(corrplot)
corrplot(cor(pp), method = "ellipse")

plot(density(chain1$lnL+chain1$prior), xlim=c(-800, -400))
lines(density(chain2$lnL+chain2$prior), col="blue")
lines(density(chain3$lnL+chain3$prior), col='red')
lines(density(chain4$lnL+chain4$prior), col="green")
allsb <- c(unlist(chain1$sb), unlist(chain2$sb), unlist(chain3$sb), unlist(chain4$sb))
plot(c(0, length(allsb)+2e5), c(0, nrow(tree$edge)+1), type="n")
abline(h=clades, col="yellow", lwd=2)
points(allsb, pch=".")
text(rep(length(allsb), length(clades)), clades, labels=names(clades), pos=4)
abline(v=cumsum(c(length(unlist(chain1$sb)), length(unlist(chain2$sb)),length(unlist(chain3$sb)), length(unlist(chain4$sb)))))
#identifyBranches(tree, 10)

clades <- c("Plethodontidae"=344, "Serpentes"=578, "Boidae"=501, "Microtus/Lemmus"=1222, "Chiroptera"=1618, "Passeriformes"=822, "Phrynosomatidae"=526)
combChains <- function(chains, burnins){
  applyburnins <- lapply(1:length(chains), function(x) lapply(chains[[x]], function(y) y[(round(length(y)*burnins[x], 0):length(y))]))
  nc <- names(chains[[1]])
  chain <- lapply(1:length(nc), function(x) do.call(c, lapply(applyburnins, function(y) y[[nc[x]]])))
  names(chain) <- nc
  chain
}

chain <- combChains(list(chain1, chain2, chain3, chain4), burnins=c(0.3, 0.3, 0.5, 0.5))
saveRDS(chain, file="../output/runs/tetrapods/tetrapodsCombinedChains.rds")
plot(chain$lnL)

tree <- reorder(tree, "postorder")
majorclades <- identifyBranches(tree, 5)
majorcladetips <- lapply(majorclades$sb, function(x) {node <- tree$edge[x, 2]; cache$desc$tips[[node]]})

.ancestorBranches <- function(branch, cache){
  ancbranches <- which(sapply(cache$bdesc, function(x) branch %in% x))
  sort(ancbranches, decreasing=FALSE)
}
.branchRegime <- function(branch, abranches, chain, seqx, summary=FALSE){
  ancs <- c(branch, abranches[[branch]])
  ancshifts <- lapply(1:length(seqx), function(x) chain$t2[[seqx[x]]][which(chain$sb[[seqx[x]]] == ancs[min(which(ancs %in% chain$sb[[seqx[x]]]))])])
  ancshifts <- sapply(ancshifts, function(x) ifelse(length(x)==0, 1, x))
  thetas <- sapply(1:length(ancshifts), function(x) chain$theta[[seqx[x]]][ancshifts[x]])
  beta1s <- sapply(1:length(ancshifts), function(x) chain$beta1[[seqx[x]]][ancshifts[x]])
  res <- cbind(thetas, beta1s)
  if(summary){
    return(apply(res, 2, median))
  } else {
    return(res)
  }
}


abranches <- lapply(1:nrow(tree$edge), .ancestorBranches, cache=cache)
i=5
cl1 <- "white"
cl2 <- "green"
branch = majorclades$sb[i]
nn <- 10000
seq1 <- floor(seq(1,length(chain[[1]]), length.out=nn))

res <- lapply(majorclades$sb, function(x) .branchRegime(x, abranches, chain, seq1))
lwd=2

pdf("../output/figures/Slopes&IntMajorClades.pdf")
for(i in 1:length(majorclades$sb)){
  par(bg="black")
  plot(lnMass, lnBMR, pch=21, bg=makeTransparent(cl1, 50), col=makeTransparent(cl1,50), xlab="ln Mass", ylab="ln BMR")
  axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue")
  axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue")
  points(lnMass[majorcladetips[[i]]], lnBMR[majorcladetips[[i]]], pch=21, bg=makeTransparent(cl2,150), col=makeTransparent(cl2,150))
  lapply(1:nrow(res[[i]]), function(x) abline(res[[i]][x,1], res[[i]][x, 2], lty=1, col=makeTransparent(cl2, 1)))
}
dev.off()

nn <- 10000
seq1 <- floor(seq(1,length(chain[[1]]), length.out=nn))




res <- lapply(majorclades$sb, function(x) .branchRegime(x, abranches, chain, seq1))
res2 <- lapply(clades, function(x) .branchRegime(x, abranches, chain, seq1))
lwd=2

par(mfrow=c(1,2), bg="black")
majorcladecolors <- c("lightblue", "green", "orange", "purple", "red")
plot(density(res[[i]][,1]), type="n", xlim=c(-1, 2), ylim=c(0,7), main="Intercept", xlab="", ylab="")
lapply(1:length(majorclades$sb), function(i) polygon(density(res[[i]][,1]), col = makeTransparent(majorcladecolors[i], 100), border =makeTransparent(majorcladecolors[i], 250) ))
mtext("Intercept", side=3, cex = 2, col="skyblue")
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue")
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue")
plot(density(res[[i]][,1]), type="n", xlim=c(0.6, 0.9), ylim=c(0,38), main="Intercept", xlab="", ylab="")
abline(v=c(2/3, 0.75), lty=2, lwd=lwd, col="skyblue")
lapply(1:length(majorclades$sb), function(i) polygon(density(res[[i]][,2]), col = makeTransparent(majorcladecolors[i], 100), border =makeTransparent(majorcladecolors[i], 250) ))
mtext("Slope", side=3, cex = 2, col="skyblue")
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue")
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue")

pdf("../output/figures/densityplotsMinorClades.pdf", height=6, width=10)
for(j in 1:length(clades)){
  par(mfrow=c(1,2), mar=c(5.1, 2.1, 4.1, 1), bg="black", ask=FALSE)
  majorcladecolors <- c("lightblue", "green", "orange", "purple", "red")
  plot(density(res[[i]][,1]), type="n", xlim=c(-2, 2), ylim=c(0,7), main="Intercept", xlab="", ylab="")
  lapply(1:length(majorclades$sb), function(i) polygon(density(res[[i]][,1]), col = makeTransparent("white", 100), border =makeTransparent("white", 250) ))
  polygon(density(res2[[j]][,1]), col = makeTransparent("green", 100), border =makeTransparent("green", 250) )
  mtext("Intercept", side=3, cex = 2, col="skyblue")
  axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue")
  axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue")
  par(mar=c(5.1, 1, 4.1, 2.1))
  plot(density(res[[i]][,1]), type="n", xlim=c(0.6, 0.9), ylim=c(0,38), main="Intercept", xlab="", ylab="")
  abline(v=c(2/3, 0.75), lty=2, lwd=lwd, col="skyblue")
  lapply(1:length(majorclades$sb), function(i) polygon(density(res[[i]][,2]), col = makeTransparent("white", 100), border =makeTransparent("white", 250) ))
  polygon(density(res2[[j]][,2]), col = makeTransparent("green", 100), border =makeTransparent("green", 250) )
  mtext("Slope", side=3, cex = 2, col="skyblue")
  axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue")
  axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue")
  print(names(clades)[j])
}
dev.off()  


polygon(density(res[[i]][,1]), col = makeTransparent(cl2, 100), border =makeTransparent(cl2, 250) )
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue")
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue")

  
allbranches <- sapply(1:nrow(tree$edge), function(x) .branchRegime(x, abranches, chain, seq1, summary=TRUE))
allbranches <- t(allbranches)
colorRamp <- function(trait, pal,nn){
  strait <- (trait-min(trait))/max(trait-min(trait))
  itrait <- round(strait*nn, 0)+1
  return(pal(nn+1)[itrait])
}

pdf("../output/figures/BMRslopeIntecept.pdf")
par(bg="black")
plot(tree, edge.color=colorRamp(allbranches[,1], rainbow, 100), cex=0.3,edge.width=1.5)
plot(tree, edge.color=colorRamp(allbranches[,2], rainbow, 100), cex=0.3,edge.width=1.5)
dev.off()

nn <- 10000
burnin <- 0.3


taxa <- c("fish","amphibians", "squamates", "birds",  "newmammals")
Allbranches <- list()

for(i in taxa){
  chain <- readRDS(paste("../output/runs/",i, "/",i,"chain.rds", sep=""))
  mymcmc <- readRDS(paste("../output/runs/",i, "/",i,"mcmc.rds", sep=""))
  seq1 <- floor(seq(burnin*length(chain[[1]]),length(chain[[1]]), length.out=nn))
  tree <- mymcmc$tree
  dat <- mymcmc$dat
  pred <- mymcmc$pred
  cache <- bayou:::.prepare.ou.univariate(tree, dat, SE=0, pred=pred)
  abranches <- lapply(1:nrow(cache$edge), .ancestorBranches, cache=cache)
  Allbranches[[i]] <- sapply(1:nrow(tree$edge), function(x) .branchRegime(x, abranches, chain, seq1, summary=TRUE))
  Allbranches[[i]] <- t(Allbranches[[i]])
  dir <- paste("../output/figures/", sep="")
  #dir.create(dir)
  pdf(paste(dir ,i,"Int&Slope.pdf", sep=""))
  par(bg="black")
  plot(tree, edge.color=colorRamp(Allbranches[[i]][,1], cm.colors, 100), cex=0.3,edge.width=1.5)
  plot(tree, edge.color=colorRamp(Allbranches[[i]][,2], cm.colors, 100), cex=0.3,edge.width=1.5)
  dev.off()
}

### Plethodon, Snakes and Boids
tmp <- identifyBranches(tree, 3)
minorclades <- c("Plethodontidae"=344, "Serpentes"=578, "Colubridae"=456, "Boidae"=501, Passeriformes=822, "Chiroptera"=1618, Phrynosomatidae=528, Dendrobatidae=138)

plotPosteriors2 <- function(tree, dat, pred, chain, burnin=0.3, c() thin=1){
  getBranchPosterior <- function(sb, chain, pars=c("beta1", "theta"), burnin=0.3, thin=100){
    postburn <- round(seq(burnin*length(chain$gen), length(chain$gen), length.out=length(chain$gen)/thin),0)
    pbchain <- lapply(chain, function(x) x[postburn])
    if(sb==0){
      out <- lapply(pars, function(x) sapply(pbchain[[x]], function(y) y[1]))
    } else {
      sbs <- unlist(pbchain$sb)
      out <- lapply(pars, function(x) unlist(lapply(pbchain[[x]], function(y) y[-1]),F,F)[sbs==sb])
    }
    return(out)
  }
  LP <- Lposterior(chain, tree, burnin=burnin)
  focalSB <- c(minorclades)
  PostSB <- lapply(focalSB, function (y) getBranchPosterior(y, chain, burnin=burnin, thin=thin))
  postBetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[1]]))
  postThetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[2]]))
  postBetas <- data.frame(do.call(rbind, postBetas))
  postBetas$var <- rep("Beta", nrow(postBetas))
  postThetas <- data.frame(do.call(rbind, postThetas))
  postThetas$var <- rep("Theta", nrow(postThetas))
  allpar <- rbind(postBetas, postThetas)
  allpar$x <- factor(allpar$x)
  focalpars <- list(k=length(focalSB[-1]), 
                    ntheta=length(focalSB), 
                    sb=focalSB[-1], 
                    loc=rep(0, length(focalSB[-1])), 
                    t2=2:length(focalSB))
  tr <- pars2simmap(focalpars, tree)
  postburn <- round(seq(burnin*length(chain$gen), length(chain$gen), length.out=length(chain$gen)/thin),0)
  pbchain <- lapply(chain, function(x) x[postburn])
  p <- ggplot(allpar, aes(V2, fill=x))
  p <- p + geom_density(alpha=0.5)
  p <- p + facet_wrap(~ var, scales="free")
  p <- p + scale_fill_manual(values=tr$col)
  #p <- p + ylim(c(0,5))
  #p <- p + xlim(c(0.6, 0.8))
  print(p + geom_vline(x=c(3/4, 2/3), lty=2))
  medianBetas <- tapply(postBetas[,2], postBetas[,1], median)
  medianThetas <- tapply(postThetas[,2], postThetas[,1], median)
  est.pars <- list(beta1=medianBetas, 
                   thetas=medianThetas, 
                   k=length(focalSB)-1, 
                   ntheta=length(focalSB), 
                   sb=focalSB[-1], 
                   loc=rep(0, length(focalSB)-1), 
                   t2=2:length(focalSB))
  tipreg <- bayou:::.tipregime(est.pars, tree)
  plotSimmap.mcmc(tr$tree, pbchain, colors=tr$col, burnin=burnin, fsize=0.1)
  plot(pred[,1], dat, type="n", main=i, 
       xlab="Mass", ylab="BMR", 
       xlim=c(min(pred), max(pred))*c(0.9,1.1), 
       ylim=c(min(dat), max(dat))*c(0.9,1.1))
  for(i in 1:est.pars$ntheta){
    points(pred[tipreg==i,],dat[tipreg==i], pch=21, 
           bg=bayou:::makeTransparent(tr$col[i],50), col=bayou:::makeTransparent(tr$col[i],50))
    abline(a=est.pars$theta[i], b=est.pars$beta1[i], lwd=2, col=tr$col[i],lty=2)
  }
}


allpar

majorcladeDensities <- lapply(res, function(x) apply(x, 2, density))
minorcladeDensities <- lapply(1:length(minorclades), function(x) list(thetas=density(allpar[allpar[,1]==x & allpar[,3]=="Beta", 2]), beta1s=density(allpar[allpar[,1]==x & allpar[,3]=="Theta", 2])))
minorcladecolors <- c("skyblue", "yellow", "white")

par(bg="black", mfrow=c(1,2))
plot(c(0.65, 0.83), c(0, 40), type="n")
lapply(1:length(majorclades$sb), function(i) polygon(majorcladeDensities[[i]][[2]], col = makeTransparent(majorcladecolors[i], 100), border =makeTransparent(majorcladecolors[i], 250) ))
lapply(1:length(minorcladeDensities),function(i) polygon(minorcladeDensities[[i]][[1]], col = makeTransparent(minorcladecolors[i], 100), border=makeTransparent(minorcladecolors[i], 250) ))

plot(c(-2, 2), c(0, 6), type="n")
lapply(1:length(majorclades$sb), function(i) polygon(majorcladeDensities[[i]][[1]], col = makeTransparent(majorcladecolors[i], 100), border =makeTransparent(majorcladecolors[i], 250) ))
lapply(1:length(minorcladeDensities),function(i) polygon(minorcladeDensities[[i]][[2]], col = makeTransparent(minorcladecolors[i], 100), border=makeTransparent(minorcladecolors[i], 250) ))

require(MASS)
pdf("../output/figures/Plethodontidae.pdf")
i=1
par(bg="black")
plot(do.call(rbind, res), type="n", ylim=c(0.6, 0.85), xlim=c(-2, 2))
minorclades
z <- kde2d(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], n=60)
contour(z, col="skyblue", add=TRUE,nlevels=6, lwd=4, zlim = c(2,max(z$z)), drawlabels = FALSE, method = "simple")
#points(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], pch=21, col=makeTransparent("skyblue", 5) ,bg=makeTransparent("skyblue", 5))

mtext(side=3, text=paste(names(minorclades)[i], ":", round(LP[minorclades[i],1],2)), col="skyblue", cex=2)
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue", cex=2)
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue", cex=2)
for(j in 1:5){
  z <- kde2d(res[[j]][,1], res[[j]][,2])
  contour(z, col=majorcladecolors[j], add=TRUE, lwd=4, nlevels=4, levels=c(10, 20, 30, 50, 80,120, 200), drawlabels=FALSE)
}

dev.off()


plotPosteriorRegressions <- function(node, root=FALSE, nn=1000, chain, alpha=5,alphapts=200, add=FALSE, col="green",linecol=NULL){
  if(!add){
    plot(lnMass, lnBMR, pch=21, col=makeTransparent("white", 25), bg=makeTransparent("white",25))
    axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue", cex=2)
    axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue", cex=2)
  }
  id <- tree$tip.label[cache$desc$tips[[tree$edge[node,2]]]]
  if(root){
    thetas <- sapply(chain$theta, function(x) x[1])
    beta1s <- sapply(chain$beta1, function(x) x[1])
  } else {
  tmpid <- sapply(chain$sb, function(x) which(node==x))
  t2s <- sapply(1:length(tmpid), function(x) ifelse(length(tmpid[[x]]) > 0, chain$t2[[x]][tmpid[[x]]], NA))
  thetas <- sapply(1:length(t2s), function(x) chain$theta[[x]][t2s[x]])
  beta1s <- sapply(1:length(t2s), function(x) chain$beta1[[x]][t2s[x]])
  }
  res <- cbind(thetas, beta1s)[!is.na(thetas),]
  sq <- floor(seq(1, nrow(res), length.out=nn))
  if(length(sq)> nrow(res)){sq <- 1:nrow(res)}
  if(is.null(linecol)){
    linecol=col
  }
  lapply(sq, function(x) abline(res[x,1], res[x,2], lwd=1, col=makeTransparent(linecol, alpha)))
  points(lnMass[id], lnBMR[id], pch=21, col="white", bg=makeTransparent(col,alphapts))
  
}

pdf("../output/figures/MinorCladeRegressions.pdf")
par(bg="black")
plotPosteriorRegressions(1706, TRUE, 1000, chain, 2,alphapts=200, add=FALSE, col="green")
plotPosteriorRegressions(minorclades[1],FALSE, 2000, chain, 5,alphapts=250, add=TRUE, col="skyblue")

plotPosteriorRegressions(1706, TRUE, 1000, chain, 2,alphapts=200, add=FALSE, col="green")
plotPosteriorRegressions(minorclades[8],FALSE, 2000, chain, 5,alphapts=250, add=TRUE, col="skyblue")

plotPosteriorRegressions(842, FALSE, 5000, chain, 2,alphapts=200, add=FALSE, col="green")
plotPosteriorRegressions(minorclades[3],FALSE, 2000, chain, 5,alphapts=250, add=TRUE, col="skyblue")

plotPosteriorRegressions(842, FALSE, 5000, chain, 2,alphapts=200, add=FALSE, col="green")
plotPosteriorRegressions(minorclades[3],FALSE, 1000, chain, 5,alphapts=250, add=TRUE, col="skyblue")
plotPosteriorRegressions(462,FALSE, 1000, chain, 5,alphapts=250, add=TRUE, col="skyblue")

plotPosteriorRegressions(842, FALSE, 5000, chain, 2,alphapts=200, add=FALSE, col="green")
plotPosteriorRegressions(minorclades[2],FALSE, 1, chain, 1,alphapts=250, add=TRUE, col="green", linecol="skyblue")
plotPosteriorRegressions(minorclades[4],FALSE, 1000, chain, 5,alphapts=250, add=TRUE, col="skyblue")
dev.off()



pdf("../output/figures/Dendrobatidae.pdf")
i=8


par(bg="black")
plot(do.call(rbind, res), type="n", ylim=c(0.6, 0.85), xlim=c(-2, 2))
minorclades
z <- kde2d(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], n=60)
contour(z, col="skyblue", add=TRUE,nlevels=6, lwd=4, zlim = c(2,max(z$z)), drawlabels = FALSE, method = "simple")
#points(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], pch=21, col=makeTransparent("skyblue", 5) ,bg=makeTransparent("skyblue", 5))

mtext(side=3, text=paste(names(minorclades)[i], ":", round(LP[minorclades[i],1],2)), col="skyblue", cex=2)
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue", cex=2)
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue", cex=2)
for(j in 1:5){
  z <- kde2d(res[[j]][,1], res[[j]][,2])
  contour(z, col=majorcladecolors[j], add=TRUE, lwd=4, nlevels=4, levels=c(10, 20, 30, 50, 80,120, 200), drawlabels=FALSE)
}
dev.off()




supported <- which(LP[,1] > 0.2)
desctips <- lapply(supported, function(x) tree$tip.label[cache$desc$tips[[tree$edge[x, 2]]]])
focalSB <- supported
PostSB <- lapply(focalSB, function (y) getBranchPosterior(y, chain, burnin=burnin, thin=thin))
postBetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[1]]))
postThetas <- lapply(1:length(PostSB), function(x) cbind(x, PostSB[[x]][[2]]))
postBetas <- data.frame(do.call(rbind, postBetas))
postBetas$var <- rep("Beta", nrow(postBetas))
postThetas <- data.frame(do.call(rbind, postThetas))
postThetas$var <- rep("Theta", nrow(postThetas))
allpar <- rbind(postBetas, postThetas)
allpar$x <- factor(allpar$x)


for(i in 1:length(supported)){
par(bg="black", ask=TRUE)
plot(do.call(rbind, res), type="n", ylim=c(0.6, 0.85), xlim=c(-2, 2))
minorclades
z <- kde2d(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], n=60)
contour(z, col="skyblue", add=TRUE,nlevels=6, lwd=4, zlim = c(0,max(z$z)), drawlabels = FALSE, method = "simple")
print(desctips[[i]])
#points(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], pch=21, col=makeTransparent("skyblue", 5) ,bg=makeTransparent("skyblue", 5))

mtext(side=3, text=paste(supported[i], ":", round(LP[supported[i],1],2)), col="skyblue", cex=2)
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue", cex=2)
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue", cex=2)
for(j in 1:5){
  z <- kde2d(res[[j]][,1], res[[j]][,2])
  contour(z, col=majorcladecolors[j], add=TRUE, lwd=4, nlevels=4, levels=c(10, 20, 30, 50, 80,120, 200), drawlabels=FALSE)
}

}

Colubridae=456
ColubridaeElapidaeViperidae=462
Anguidae = 572




pdf("../output/figures/Colubridae+.pdf")
par(bg="black")
plot(do.call(rbind, res), type="n", ylim=c(0.6, 0.85), xlim=c(-2, 2))
z <- kde2d(allpar[allpar[,1] %in% c(10,11) & allpar[,3]=="Theta", 2], allpar[allpar[,1] %in% c(10,11) & allpar[,3]=="Beta", 2], n=60)
contour(z, col="skyblue", add=TRUE,nlevels=6, lwd=4, zlim = c(2,max(z$z)), drawlabels = FALSE, method = "simple")
#points(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], pch=21, col=makeTransparent("skyblue", 5) ,bg=makeTransparent("skyblue", 5))

mtext(side=3, text=paste("Colubridae+" ,":", round(sum(cumserpentes)/length(cumserpentes),2)), col="skyblue", cex=2)
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue", cex=2)
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue", cex=2)
for(j in 1:5){
  z <- kde2d(res[[j]][,1], res[[j]][,2])
  contour(z, col=majorcladecolors[j], add=TRUE, lwd=4, nlevels=4, levels=c(10, 20, 30, 50, 80,120, 200), drawlabels=FALSE)
}
dev.off()

pdf("../output/figures/BigFive.pdf")
for(k in 1:5){
par(bg="black")
plot(do.call(rbind, res), type="n", ylim=c(0.6, 0.85), xlim=c(-2, 2))
#z <- kde2d(allpar[allpar[,1] %in% c(10,11) & allpar[,3]=="Theta", 2], allpar[allpar[,1] %in% c(10,11) & allpar[,3]=="Beta", 2], n=60)
#contour(z, col="skyblue", add=TRUE,nlevels=6, lwd=4, zlim = c(2,max(z$z)), drawlabels = FALSE, method = "simple")
#points(allpar[allpar[,1]==i & allpar[,3]=="Theta", 2], allpar[allpar[,1]==i & allpar[,3]=="Beta", 2], pch=21, col=makeTransparent("skyblue", 5) ,bg=makeTransparent("skyblue", 5))

#mtext(side=3, text=paste("Colubridae+" ,":", round(sum(cumserpentes)/length(cumserpentes),2)), col="skyblue", cex=2)
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue", cex=2)
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue", cex=2)
for(j in 1:k){
  z <- kde2d(res[[j]][,1], res[[j]][,2])
  contour(z, col=majorcladecolors[j], add=TRUE, lwd=4, nlevels=4, levels=c(10, 20, 30, 50, 80,120, 200), drawlabels=FALSE)
}
}
dev.off()


pdf("../output/figures/allBMR.pdf")
par(bg="black")
plot(lnMass, lnBMR, pch=21, bg = makeTransparent("green",150), col= makeTransparent("green",150))
axis(side = 1, lwd = lwd, col="skyblue", col.axis="skyblue", cex=2)
axis(side = 2, lwd = lwd,col="skyblue", col.axis="skyblue", cex=2)

plot(tree, edge.color="green", show.tip.label=FALSE, lwd=2)
dev.off()




