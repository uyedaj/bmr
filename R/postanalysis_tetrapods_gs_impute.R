## Load tree and data into R
require(ape) 
require(aRbor)
require(devtools)
#load_all("~/repos/bayou/bayou_1.0")
#install_github("uyedaj/bayou", ref="dev")
require(bayou)
setwd("~/repos/bmr/R")

rr <- sapply(1:8, function(x) paste("impute_r", x, sep=""))

chain1 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[1], ".chain.rds", sep=""))
chain2 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[2], ".chain.rds", sep=""))
chain3 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[3], ".chain.rds", sep=""))
chain4 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[4], ".chain.rds", sep=""))
chain5 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[5], ".chain.rds", sep=""))
chain6 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[6], ".chain.rds", sep=""))
chain7 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[7], ".chain.rds", sep=""))
chain8 <- readRDS(paste("../output/runs/tetrapods_gs/tetrapods_gs_", rr[8], ".chain.rds", sep=""))
mymcmc <- readRDS("../output/runs/tetrapods_gs/tetrapods_gs_impute_r1.mcmc.rds")
tree <- mymcmc$tree
dat <- mymcmc$dat
pred <- mymcmc$pred
chains <- list(chain1, chain2, chain3, chain4, chain5, chain6, chain7, chain8)
Lpps <- lapply(chains, function(x) Lposterior(x, tree, burnin=0.3))
pairs(sapply(Lpps, function(x) x$pp))

chainA <- combine.chains(chain1, chain2, burnin.prop=0.3)
chainB <- combine.chains(chain3, chain4, burnin.prop=0.3)
chainC <- combine.chains(chain5, chain6, burnin.prop=0.3)
chainD <- combine.chains(chain7, chain8, burnin.prop=0.3)
chains1 <- combine.chains(chainA, chainB, burnin.prop=0)
chains2 <- combine.chains(chainC, chainD, burnin.prop=0)
chains <- combine.chains(chains1, chains2)
rm("chainsA", "chainsB", "chainsC", "chainsD", "chains1","chains2")
chainsmp <- chains
chains <- chains[-which(names(chains)=="missing.pred")]
names(chains) <- names(chain1)
chains$gen <- 1:length(chains$gen)
attributes(chains)$model <- attributes(chain1)$model; attributes(chains)$model.pars <- attributes(chain1)$model.pars; attributes(chains)$tree <- attributes(chain1)$tree; attributes(chains)$dat <- attributes(chain1)$dat; attributes(chains)$class <- attributes(chain1)$class; attributes(chains)$burnin <-0
sumstats <- summary(chains)
plot(chains)

plot(unlist(chains$beta2))

plotSimmap.mcmc(tree, chains, burnin=0)

#plot(unlist(chain_gs$beta3))
#plot(density(unlist(chain_gs$beta3)), main="Coefficient on genome size")

cutoff <- 0.3
sumpars <- list(sb = which(sumstats$branch.posteriors$pp > cutoff))
sumpars$k <- length(sumpars$sb)
sumpars$ntheta <- length(sumpars$sb)+1
sumpars$loc <- rep(0, sumpars$k)
sumpars$t2 <- 2:sumpars$ntheta
tr <- pars2simmap(sumpars, tree)
plotSimmap(tr$tree, colors=tr$col, fsize=0.25)

summarizeDerivedState <- function(branch, chain){
  if(branch==0){
    Th <- sapply(chain$theta, function(x) x[1])
    B1 <- sapply(chain$beta1, function(x) x[1])
    B2 <- sapply(chain$beta2, function(x) x[1])
    #B3 <- unlist(chain$beta3)
  } else {
    SB <- unlist(chain$sb)
    gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
    ind <- which(SB==branch)
    gen <- gen[ind]
    T2 <- unlist(chain$t2)[ind]
    B1 <- sapply(1:length(T2), function(x) chain$beta1[[gen[x]]][T2[x]])
    Th <- sapply(1:length(T2), function(x) chain$theta[[gen[x]]][T2[x]])
    B2 <- unlist(chain$beta2)[gen]
    #B3 <- unlist(chain$beta3)[gen]
  }
  medians = list(theta=median(Th), beta1=median(B1), beta2=median(B2))
  densities = list(theta=density(Th), beta1=density(B1), beta2=density(B2))
  return(list(medians=medians, densities=densities))
}

sb <- sumpars$sb
cladesummaries <- lapply(c(0, sb), function(x) summarizeDerivedState(x, chains))
regressions <- t(sapply(cladesummaries, function(x) unlist(x$medians)))
rownames(regressions) = c("root", sumpars$sb[sb])
cache <- bayou:::.prepare.ou.univariate(tree, dat, SE=0, pred)
tipregs <- bayou:::.tipregime(sumpars, tree)
descendents <- lapply(1:(length(sb)+1), function(x) names(tipregs[tipregs==x])) 
#descendents <- lapply(sb, function(x) na.exclude(cache$tip.label[cache$edge[c(cache$bdesc[[x]], x), 2]] ))
nodesc <- which(sapply(descendents, length)==0)


.ancestorBranches <- function(branch, cache){
  ancbranches <- which(sapply(cache$bdesc, function(x) branch %in% x))
  sort(ancbranches, decreasing=FALSE)
}
.branchRegime <- function(branch, abranches, chain, parameter, seqx, summary=FALSE){
  ancs <- c(branch, abranches[[branch]])
  ancshifts <- lapply(1:length(seqx), function(x) chain$t2[[seqx[x]]][which(chain$sb[[seqx[x]]] == ancs[min(which(ancs %in% chain$sb[[seqx[x]]]))])])
  ancshifts <- sapply(ancshifts, function(x) ifelse(length(x)==0, 1, x))
  ests <- sapply(1:length(ancshifts), function(x) chain[[parameter]][[seqx[x]]][ancshifts[x]])
  res <- cbind(ests)
  if(summary){
    return(apply(res, 2, median))
  } else {
    return(res)
  }
}


colorRamp <- function(trait, pal, nn){
  strait <- (trait-min(trait))/max(trait-min(trait))
  itrait <- round(strait*nn, 0)+1
  return(pal(nn+1)[itrait])
}
addColorBar <- function(x, y, height, width, pal, trait, ticks, adjx=0, n=100,cex.lab=1,pos=2, text.col="black"){
  legend_image <- as.raster(matrix(rev(pal(n)),ncol = 1))
  #text(x = 1.5, y = round(seq(range(ave.Div)[1], range(ave.Div)[2], l = 5), 2), labels = seq(range(ave.Div)[1], range(ave.Div)[2], l = 5))
  seqtrait <- seq(min(trait), max(trait), length.out=nrow(legend_image))
  mincut <- n-which(abs(seqtrait - min(ticks))==min(abs(seqtrait-min(ticks))))
  maxcut <- n-which(abs(seqtrait - max(ticks))==min(abs(seqtrait-max(ticks))))
  legend_cut <- legend_image[maxcut:mincut,]
  legend_cut <- rbind(matrix(rep(legend_image[1,1],round(0.05*n,0)),ncol=1), legend_cut)
  rasterImage(legend_cut, x, y, x+width, y+height)
  ticklab <- format(ticks, digits=2, trim=TRUE)
  ticklab[length(ticklab)] <- paste(">", ticklab[length(ticklab)], sep="")
  text(x+adjx, y=seq(y, y+height, length.out=length(ticks)), labels=ticklab, pos=pos,cex=cex.lab, col=text.col)
}

plotBranchHeatMap <- function(tree, chain, variable, burnin=0, nn=NULL, pal, legend_ticks, ...){
  seq1 <- floor(max(seq(burnin*length(chain$gen),1), length(chain$gen), 1))
  if(is.null(nn)) nn <- length(seq1) else { seq1 <- floor(seq(max(burnin*length(chain$gen),1), length(chain$gen), length.out=nn))}
  if(length(nn) > length(chain$gen)) stop("Number of samples greater than chain length, lower nn")
  abranches <- lapply(1:nrow(tree$edge), .ancestorBranches, cache=cache)
  allbranches <- sapply(1:nrow(tree$edge), function(x) .branchRegime(x, abranches, chain, variable, seq1, summary=TRUE))
  plot(tree, edge.color=colorRamp(allbranches, pal, 100), ...)
  addColorBar(x=470, y=100, height=150, width=10, pal=pal, n=100, trait=allbranches, ticks=legend_ticks,adjx=25, cex.lab=.5, text.col="white")
}
require(wesanderson)
palB <- wes_palette("FantasticFox")
pal <- colorRampPalette(rev(c(palB[5],palB[4],"white", palB[2], "yellow2")))

pdf("../output/figures/tetrapods_gs_impute_ch1-8_heatposteriors.pdf", width=15, height=8)
par(mfrow=c(1,3), mar=c(0,0,0,0), bg="black")
#lapply(chains, function(x) plotBranchHeatMap(tree, x, "beta1", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(0.5, 1.3, 0.2)))
plotBranchHeatMap(tree, chains, "beta1", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(0.65, 0.85, 0.05), edge.width=2)
text(50, 100, labels = "Beta1", cex=2, col = "white")
plotBranchHeatMap(tree, chains, "beta2", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(-0.005, 0.005, 0.002), edge.width=2)
text(50, 100, labels = "Beta2", cex=2, col = "white")
plotBranchHeatMap(tree, chains, "theta", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(-3, 2.5, 0.5), edge.width=2)
text(50, 100, labels = "Theta", cex=2, col = "white")
dev.off()


pdf("../output/figures/tetrapods_gs.pdf", height=12, width=8)
plot(chain_gs)

par(mfrow=c(2,2), mar=c(3, 3, 1, 0))
plot(c(-4, 4), c(0, 4), type="n", main="Theta", ylab="density", xlab="Theta")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$theta, col=x))
plot(c(0.4, 1), c(0, 15), type="n", main="Beta1", ylab="density", xlab="Beta1")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$beta1, col=x))
plot(c(-0.05, 0.05), c(0, 110), type="n", main="Beta2", ylab="density", xlab="Beta2")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$beta2, col=x))
plot(c(-0.1, 0.5), c(0, 10), type="n", main="Beta3", ylab="density", xlab="Beta3")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$beta3, col=x))

pdf("../output/figures/regressions_tetrapods_gs_impute_ch1-8.pdf", height=8, width=12)
par(mfrow=c(2,3), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
c1 <- "pink"
palx <- rainbow
beta3 <- median(chains$beta3)
missing.pred <- do.call(rbind, chainsmp$missing.pred)
mps <- apply(missing.pred, 2, median)
impPred <- pred
impPred[is.na(pred[,3]),3] <- mps
for(i in (1:nrow(regressions))){
  plotBayoupars(sumpars, tree, colors=setNames(c(palx(nrow(regressions))[i], rep("gray80", nrow(regressions)-1)), c(i, (1:nrow(regressions))[-i])), fsize=0.2)
  if(i==1){
    plot(density(chains$beta3), lwd=3, main="Beta3", col="gray80", xlim=c(-0.25, 0.25))
    abline(v=0, lty=2, lwd=2, col="gray80")
  } else {
    barplot(setNames(sapply(Lpps, function(x) x[sumpars$sb[i-1],1]), c("chain1", "chain2", "chain3", "chain4")), ylim=c(0, 1), col="gray80", main="Posterior Probabilities")
  }
  plot(pred[,1], dat, xlab="lnMass", ylab="lnBMR", pch=21, bg=makeTransparent("gray20", 100), col =makeTransparent("gray20", 100) )
  include <- which(names(dat) %in% descendents[[i]])
  text(pred[include, 1], dat[include], labels=names(dat[include]), col="white", cex=0.4, pos = 2)
  points(pred[include,1], dat[include], pch=21, bg=palx(nrow(regressions))[i])
  print(descendents[[i]])
  expected <- regressions[i,1]+regressions[i,2]*pred[include,1]+regressions[i,3]*pred[include,2]+beta3*impPred[include,3]
  o <- order(pred[include,1])
  lines(pred[include,1][o], expected[o], col=palx(nrow(regressions))[i], lwd=2)
  plot(cladesummaries[[i]]$densities$beta2, col=palx(nrow(regressions))[i], xlim=c(-0.05, 0.05), lwd=3, main="Beta2")
  abline(v=0, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(cladesummaries[[i]]$densities$beta1, col=palx(nrow(regressions))[i], xlim=c(0.5, 1.25), lwd=3, main="Beta1")
  abline(v=0.75, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(cladesummaries[[i]]$densities$theta, col=palx(nrow(regressions))[i], xlim=c(-2.5,3), lwd=3, main="Theta")
  
}
dev.off()

#legend(-2, 10, legend=c("Root", "Plethodontids", "Birds", "Mammals"), col=1:4, lwd=2)

plotSimmap.mcmc(tree, chain_gs, burnin=0, fsize=0.5)
dev.off()

tmp <- lm(dat~ pred[,1]+pred[,2])
plot(log(pred[,3]), tmp$resid)
curve(-0.0247*x, add=TRUE)


pdf("tetrapods_gs_pruned_PosteriorShifts.pdf", height=12, width=16)
par(mfrow=c(2,3))
lapply(chains, function(x) plotSimmap.mcmc(tree, x, burnin=0.3, ftype="off"))
lapply(chains, function(x) plot(density(x$beta3), xlim=c(-0.1, 0.5)))

dev.off()
