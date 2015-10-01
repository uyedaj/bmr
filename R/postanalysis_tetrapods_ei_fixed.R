## Load tree and data into R
require(ape) 
require(aRbor)
require(devtools)
#load_all("~/repos/bayou/bayou_1.0")
#install_github("uyedaj/bayou", ref="dev")
require(bayou)
setwd("~/repos/bmr/R")

rr <- c("_fixed_r001", "_fixed_r002", "_fixed_r003", "_fixed_r004")#, "_u005", "_u006")

chain1 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[1], ".chain.rds", sep=""))
chain2 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[2], ".chain.rds", sep=""))
chain3 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[3], ".chain.rds", sep=""))
chain4 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[4], ".chain.rds", sep=""))
#chain5 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[5], ".chain.rds", sep=""))
#chain6 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[6], ".chain.rds", sep=""))
mymcmc_gs <- readRDS("../output/runs/tetrapods_ei/tetrapods_ei__t001.mcmc.rds")
tree <- mymcmc_gs$tree
tree <- reorder(tree, "postorder")
dat <- mymcmc_gs$dat
pred <- mymcmc_gs$pred
chains <- list(chain1, chain2, chain3, chain4)#, chain5, chain6)
Lpps <- lapply(chains, function(x) Lposterior(x, tree, burnin=0.3))
pairs(sapply(Lpps, function(x) x$pp))

chain_gs1 <- combine.chains(chain1, chain2, burnin.prop=0.3)
chain_gs2 <- combine.chains(chain3, chain4, burnin.prop=0.3)
#chain_gs3 <- combine.chains(chain5, chain6, burnin.prop=0.3)
chain <- combine.chains(chain_gs1, chain_gs2, burnin.prop=0)
#chain <- lapply(1:length(chain), function(x) c(chain[[x]], chain_gs3[[x]]))
names(chain) <- names(chain_gs1)
#chain_gs <- combine.chains(chain_gs, chain_gs3, burnin.prop=0)
thinseq <- floor(seq(1, length(chain$gen), length.out=10000))
chain$gen <- 1:length(chain$gen)
attributes(chain)$model <- attributes(chain1)$model; attributes(chain)$model.pars <- attributes(chain1)$model.pars; attributes(chain)$tree <- attributes(chain1)$tree; attributes(chain)$dat <- attributes(chain1)$dat; attributes(chain)$class <- attributes(chain1)$class; attributes(chain)$burnin <-0
plot(chain)
sumstats <- summary(chain)

#pdf("../output/figures/tetrapods_gs_simmapmcmc.pdf", height=15, width=10)
#plotSimmap.mcmc(tree, chain, burnin=0.3, fsize=0.5)
#dev.off()

cutoff <- 0.2
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
    #B2 <- sapply(1:length(T2), function(x) chain$beta2[[gen[x]]][T2[x]])
    #B3 <- unlist(chain$beta3)[gen]
  }
  medians = list(theta=median(Th), beta1=median(B1))#, beta2=median(B2), beta3=median(B3))
  densities = list(theta=density(Th), beta1=density(B1))#, beta2=density(B2), beta3=density(B3))
  return(list(medians=medians, densities=densities))
}
sb <- sumpars$sb
cladesummaries <- lapply(c(0, sb), function(x) summarizeDerivedState(x, chain))
regressions <- t(sapply(cladesummaries, function(x) unlist(x$medians)))
rownames(regressions) = c("root", sumpars$sb[sb])
cache <- bayou:::.prepare.ou.univariate(tree, dat, SE=0, pred)
tipregs <- bayou:::.tipregime(sumpars, tree)
descendents <- lapply(1:(length(sb)+1), function(x) names(tipregs[tipregs==x])) 
#descendents <- lapply(sb, function(x) na.exclude(cache$tip.label[cache$edge[c(cache$bdesc[[x]], x), 2]] ))

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

pdf("../output/figures/tetrapods_gs_pruned_ch1-6_heatposteriors.pdf", width=15, height=8)
par(mfrow=c(1,3), mar=c(0,0,0,0), bg="black")
#lapply(chains, function(x) plotBranchHeatMap(tree, x, "beta1", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(0.5, 1.3, 0.2)))
plotBranchHeatMap(tree, chain, "beta1", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(0.65, 0.85, 0.05), edge.width=2)
text(50, 20, labels = "Beta1", cex=2, col = "white")
#plotBranchHeatMap(tree, chain, "beta2", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(-0.005, 0.005, 0.002), edge.width=2)
#text(50, 20, labels = "Beta2", cex=2, col = "white")
plotBranchHeatMap(tree, chain, "theta", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(-2.5, 2.5, 0.5), edge.width=2)
text(50, 20, labels = "Theta", cex=2, col = "white")
dev.off()

pdf("../output/figures/tetrapods_gs.pdf", height=12, width=8)
plot(chain)

par(mfrow=c(2,2), mar=c(3, 3, 1, 0), col.axis="white", col.lab="white")
plot(c(-4, 4), c(0, 4), type="n", main="Theta", ylab="density", xlab="Theta")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$theta, col=x))
plot(c(0.4, 1), c(0, 15), type="n", main="Beta1", ylab="density", xlab="Beta1")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$beta1, col=x))
#plot(c(-0.05, 0.05), c(0, 110), type="n", main="Beta2", ylab="density", xlab="Beta2")
#lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$beta2, col=x))
#plot(c(-0.1, 0.5), c(0, 10), type="n", main="Beta3", ylab="density", xlab="Beta3")
#lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$beta3, col=x))

#par(mfrow=c(1,1))
#c1 <- "pink"
#plot(pred[,1], dat, xlab="lnMass", ylab="lnBMR", pch=21, bg=makeTransparent(c1, 100), col =makeTransparent(c1, 100) )
#for(i in 1:nrow(regressions)){
#  expected <- regressions[i,1]+regressions[i,2]*pred[,1]+regressions[i,3]*pred[,2]+regressions[i,4]*pred[,3]
#  o <- order(pred[,1])
#  lines(pred[o,1], expected[o], col=i, lwd=2)
#}
#legend(-2, 10, legend=c("Root", "Plethodontids", "Birds", "Mammals"), col=1:4, lwd=2)



pdf("../output/figures/regressions_tetrapods_ei_ch1_6.pdf", height=8, width=12)
par(mfrow=c(2,3), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
c1 <- "pink"
palx <- rainbow
chains <- chain
endo <- median(chain$endo)
#beta3 <- median(chains$beta3)
#missing.pred <- do.call(rbind, chainsmp$missing.pred)
#mps <- apply(missing.pred, 2, median)
#impPred <- pred
#impPred[is.na(pred[,3]),3] <- mps
for(i in (1:nrow(regressions))){
  plotBayoupars(sumpars, tree, colors=setNames(c(palx(nrow(regressions))[i], rep("gray80", nrow(regressions)-1)), c(i, (1:nrow(regressions))[-i])), fsize=0.2)
  if(i==1){
    barplot(setNames(rep(1,6), c("chain1", "chain2", "chain3", "chain4")), ylim=c(0, 1), col="gray80", main="Posterior Probabilities")
  } else {
    barplot(setNames(sapply(Lpps, function(x) x[sumpars$sb[i-1],1]), c("chain1", "chain2", "chain3", "chain4")), ylim=c(0, 1), col="gray80", main="Posterior Probabilities")
  }
  plot(pred[,1], dat, xlab="lnMass", ylab="lnBMR", pch=21, bg=makeTransparent("gray20", 100), col =makeTransparent("gray20", 100) )
  include <- which(names(dat) %in% descendents[[i]])
  text(pred[include, 1], dat[include], labels=names(dat[include]), col="white", cex=0.4, pos = 2)
  points(pred[include,1], dat[include], pch=21, bg=palx(nrow(regressions))[i])
  print(descendents[[i]])
  expected <- regressions[i,1]+regressions[i,2]*pred[include,1]+endo*pred[include,3]
  o <- order(pred[include,1])
  lines(pred[include,1][o], expected[o], col=palx(nrow(regressions))[i], lwd=2)
  #plot(cladesummaries[[i]]$densities$beta2, col=palx(nrow(regressions))[i], xlim=c(-0.05, 0.05), lwd=3, main="Beta2")
  #abline(v=0, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(density(chain$endo), col=palx(nrow(regressions))[i], xlim=c(3.5, 5), lwd=3, main="Endothermy")
  plot(cladesummaries[[i]]$densities$beta1, col=palx(nrow(regressions))[i], xlim=c(0.5, 1.25), lwd=3, main="Beta1")
  abline(v=0.75, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(cladesummaries[[i]]$densities$theta, col=palx(nrow(regressions))[i], xlim=c(-7,-1), lwd=3, main="Theta")
  
}
dev.off()

saveRDS(sumpars, file="../output/runs/tetrapods_ei/sumpars_u6.rds")

par(mfrow=c(2,2), mar=c(3,3,5,1), bg="black", ask=TRUE, col.axis="white", col.lab="white", col.main="white")
c1 <- "pink"
palx <- rainbow
for(i in 1:nrow(regressions)){
  plot(pred[,1], dat, type="n", xlab="lnMass", ylab="lnBMR", pch=21, bg=palx(nrow(regressions))[1], col =palx(nrow(regressions))[1] )
  include <- which(names(dat) %in% descendents[[i]])
  points(pred[include,1], dat[include], pch=21, bg=palx(nrow(regressions))[i])
  print(descendents[[i]])
  expected <- regressions[i,1]+regressions[i,2]*pred[,1]+regressions[i,3]*pred[,2]+regressions[i,4]*pred[,3]
  o <- order(pred[,1])
  lines(pred[o,1], expected[o], col=palx(nrow(regressions))[i], lwd=2)
  
  plot(cladesummaries[[i]]$densities$beta2, col=palx(nrow(regressions))[i], xlim=c(-0.05, 0.05), lwd=3, main="Beta2")
  plot(cladesummaries[[i]]$densities$beta1, col=palx(nrow(regressions))[i], xlim=c(0.5, 1.25), lwd=3, main="Beta1")
  plot(cladesummaries[[i]]$densities$theta, col=palx(nrow(regressions))[i], xlim=c(-2.5,3), lwd=3, main="Theta")
  
}






plotSimmap.mcmc(tree, chain, burnin=0, fsize=0.5)
dev.off()

tmp <- lm(dat~ pred[,1]+pred[,2])
plot(log(pred[,3]), tmp$resid)
curve(-0.0247*x, add=TRUE)


pdf("tetrapods_gs_pruned_PosteriorShifts.pdf", height=12, width=16)
par(mfrow=c(2,3))
lapply(chains, function(x) plotSimmap.mcmc(tree, x, burnin=0.3, ftype="off"))
lapply(chains, function(x) plot(density(x$beta3), xlim=c(-0.1, 0.5)))

dev.off()

ldesc <- sapply(descendents, length)

thetas <- do.call(rbind, chain$theta)
rws <- sample(1:nrow(thetas), 3000, replace=FALSE)
theta.medians <- apply(thetas, 2, median)
o <- order(theta.medians)
boxplot(thetas[rws,o], pch=".")

beta1s <- do.call(rbind, chain$beta1)
rws <- sample(1:nrow(beta1s), 5000, replace=FALSE)
beta1.medians <- apply(beta1s, 2, median)
boxplot(beta1s[rws,ldesc==1])
#o <- order(beta1.medians)
#boxplot(beta1s[rws,o], pch=".")

i=27
plot(beta1s[,27])
plot(thetas[,27])
range(thetas)
range(beta1s)
for(i in 1:ncol(thetas)){
  plot(1,1, type="n", xlim=c(-8, 5), ylim=c(0.3, 1.1), main=ldesc[i]) 
  rws <- sample(1:nrow(beta1s), 5000, replace=FALSE)
  points(thetas[rws,i], beta1s[rws,i], pch=".", col=tr$col[i])
}
apply(coda::mcmc(thetas), 2, effectiveSize)
apply(coda::mcmc(beta1s), 2, effectiveSize)
plot(beta1s[,29])
iqrbeta1 <- apply(beta1s, 2, function(x) diff(quantile(x, c(0.25,0.75))))
plot(log(ldesc), iqrbeta1)

boxplot(beta1s[rws,], width=log(ldesc)+5, pch=".")

require(rotl)
require(vioplot)
require(rphylopic)
thetas <- do.call(rbind, chain$theta)
beta1s <- do.call(rbind, chain$beta1)
licas <- sapply(1:length(descendents), function(x){
  tmp <- tnrs_match_names(gsub("_"," ", descendents[[x]]), do_approximate_matching = FALSE)$ott_id;
  rotl::taxonomy_lica(ott_ids = as.numeric(tmp)[!is.na(tmp)])$lica$'ot:ottTaxonName'
  }
)
licas[licas=="Bifurcata"] <- "Squamata"
licas[licas=="Euteleostomi"] <- "Root"
licas[licas=="Sauria"] <- "Aves"
names(sb) <- licas[2:length(licas)]
colnames(thetas) <- colnames(beta1s) <- names(descendents) <- licas

bigclades <- which(ldesc > 3)

thetas[,c('Aves', 'Mammalia', 'Chiroptera')] <- thetas[,c('Aves', 'Mammalia', 'Chiroptera')]+chain$endo


pdf('../output/figures/boxplots_fixed_ei.pdf', width=8, height=10)
par(bg="black")
pal <- colorRampPalette(c("#ecf0f1", "#2ecc71", "#3498db", "#9b59b6","#e74c3c",  "#e67e22", "#f1c40f"))
set.seed(5)
cols <- setNames(pal(length(sb)+1), c(1, sample(2:(length(sb)+1), length(sb), replace=FALSE)))
cols <- cols[order(as.numeric(names(cols)))]
plotSimmap(tr$tree, colors=cols, ftype="off", mar=c(0.1,0.1,0.1,6))

par(mfrow=c(2,1), bg="black", col.axis="white", lwd=2,col.lab="white", col.main="white", col.sub="white")
par(mar=c(1,4,9,1))
o <-c(1,order(apply(beta1s[,bigclades[-1]], 2, median))+1)
#boxplot(thetas[,bigclades[o]], las=2, names = rep("", 10), pch=".", col=tr$col[bigclades[o]])
plot(0,0, type="n", xlab="", ylab="Intercept", xlim=c(0,11), ylim=c(-5,4.5), xaxt="n", cex.lab=1.25, yaxt="n")
axis(2, at=seq(-6,6,2), col="white", col.ticks="white", cex=1.25)
abline(h=seq(-6,5,1), col="gray50")
#axis(1, at=1:10, licas[bigclades[o]], las=2)
lapply(1:length(bigclades), function(x) vioplot(thetas[,bigclades[o]][,x], at=x, add=TRUE, col=cols[bigclades[o]][x], border="white"))
box(col="white")
par(mar=c(9,4,1,1))

#boxplot(beta1s[,bigclades[o]], las=2, pch=".",col=tr$col[bigclades[o]])

plot(0,0, type="n", xlab="", ylab="Slope", xlim=c(0,11), ylim=c(0.3,1.05), xaxt="n", cex.lab=1.25, yaxt="n")
abline(h=seq(0,2,0.1), col="gray50")
axis(1, at=1:10, paste(licas[bigclades[o]], " (", ldesc[bigclades[o]],")", sep=""), las=2, col="white", col.ticks = "white", cex=1.5)
axis(2, at=seq(0,1.1,0.3), col="white", col.ticks="white", cex=1.25)
abline(h=c(2/3, 0.75), lty=2, col="white", lwd=3)
lapply(1:length(bigclades), function(x) vioplot(beta1s[,bigclades[o]][,x], at=x, add=TRUE, col=cols[bigclades[o]][x], border="white", cex=2))
box(col="white")

dev.off()


pdf("../output/figures/singletons.pdf")
par(mfrow=c(1,1), mar=c(3,3,5,1),type="n", bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
plot(pred[,1], dat, xlab="lnMass", ylab="lnBMR",xaxt="n", yaxt="n", pch=21, bg="gray20", col ="gray20" )
axis(1, at=seq(-2,14,2), col="white", col.ticks="white", cex.axis=1.5)
axis(2, at=seq(-12,12,4), col="white", col.ticks="white", cex.axis=1.5)
box(lwd=1.25, col="white")
points(pred[,1], dat, pch=21, cex=1.5, bg=makeTransparent("gray90", 200), col=makeTransparent("gray50", 255))

par(mfrow=c(1,1), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
cols <- setNames(rep(makeTransparent("gray30",150), length(ldesc)), 1:30)
cols[ldesc<=3] <- "white"
plotSimmap(tr$tree, colors=cols, ftype="off",lwd=2)
singletons <- unlist(descendents[which(ldesc<=3)])
singletons <- unlist(singletons)
par(mfrow=c(1,1), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
plot(pred[,1], dat, xlab="lnMass", ylab="lnBMR",xaxt="n", yaxt="n", pch=21, bg="gray20", col ="gray20" )
axis(1, at=seq(-2,14,2), col="white", col.ticks="white", cex.axis=1.5)
axis(2, at=seq(-12,12,4), col="white", col.ticks="white", cex.axis=1.5)
box(lwd=1.25, col="white")
include <- which(names(dat) %in% singletons)
text(pred[include, 1], jitter(dat[include], 10), labels=gsub("_", " ", names(dat[include])), col="white", cex=0.3, pos = 4)
points(pred[include,1], dat[include], pch=21, cex=1.25, bg="white")
dev.off()

