pdf("../output/figures/fixedchains27.pdf")

fixed.pars <- list(k=length(sb), ntheta=length(sb)+1, sb=unname(sb), t2=2:(length(sb)+1), loc=rep(0, length(sb)))
tr <- pars2simmap(fixed.pars, cache$phy)
cols <- setNames(rainbow(length(sb))[sample(1:(length(sb)),length(sb), replace=FALSE)], 1:length(sb))
plotSimmap(tr$tree, colors=cols, fsize=0.5, ftype="off")


par(mar=c(15, 3, 1,1))
beta1 <- do.call(rbind, chains$beta1)
colnames(beta1) <- c("root", names(sb))
o <- order(apply(beta1,2, mean))
beta1 <- beta1[,o]
boxplot(beta1, las=2, pch=".", col = cols[o])
abline(h=seq(0, 1.5, 0.1), lty=2, col="gray80")
boxplot(beta1, las=2, pch=".", col = cols[o], add=TRUE)


par(mar=c(15, 3, 1,1))
theta <- do.call(rbind, chains$theta)
colnames(theta) <- c("root", names(sb))
o <- order(apply(theta,2, mean))
theta <- theta[,o]
boxplot(theta, las=2, pch=".", col=cols[o])
abline(h=seq(-20,20, 1), lty=2, col="gray80")
boxplot(theta, las=2, pch=".", col=cols[o], add=TRUE)

#plot(c(0, length(chains$gen)), c(-5, 5), type="n")
#lapply(1:ncol(theta), function(x) lines(theta[,x], col=x))
##summary(coda::mcmc(theta))
#summary(coda::mcmc(beta1))
#apply(coda::mcmc(theta), 2, effectiveSize)
#apply(coda::mcmc(beta1), 2, effectiveSize)
dev.off()



pdf("../output/figures/fixedchains27_regressions.pdf")
par(mfrow=c(2,2), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
c1 <- "pink"
palx <- rainbow
#beta3 <- median(chains$beta3)
#missing.pred <- do.call(rbind, chainsmp$missing.pred)
#mps <- apply(missing.pred, 2, median)
#impPred <- pred
#impPred[is.na(pred[,3]),3] <- mps
for(i in (1:nrow(regressions))){
  plotBayoupars(sumpars, tree, colors=setNames(c(palx(nrow(regressions))[i], rep("gray80", nrow(regressions)-1)), c(i, (1:nrow(regressions))[-i])), fsize=0.2)
  plot(pred[,1], dat, xlab="lnMass", ylab="lnBMR", pch=21, bg=makeTransparent("gray20", 100), col =makeTransparent("gray20", 100) )
  include <- which(names(dat) %in% descendents[[i]])
  text(pred[include, 1], dat[include], labels=names(dat[include]), col="white", cex=0.4, pos = 2)
  points(pred[include,1], dat[include], pch=21, bg=palx(nrow(regressions))[i])
  print(descendents[[i]])
  expected <- regressions[i,1]+regressions[i,2]*pred[include,1]#+regressions[i,3]*pred[include,2]+beta3*impPred[include,3]
  o <- order(pred[include,1])
  lines(pred[include,1][o], expected[o], col=palx(nrow(regressions))[i], lwd=2)
  #plot(cladesummaries[[i]]$densities$beta2, col=palx(nrow(regressions))[i], xlim=c(-0.05, 0.05), lwd=3, main="Beta2")
  #abline(v=0, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(cladesummaries[[i]]$densities$beta1, col=palx(nrow(regressions))[i], xlim=c(0.5, 1.25), lwd=3, main="Beta1")
  abline(v=0.75, col=palx(nrow(regressions))[i], lty=2, lwd=2)
  plot(cladesummaries[[i]]$densities$theta, col=palx(nrow(regressions))[i], xlim=c(-2.5,3), lwd=3, main="Theta")
  
}
dev.off()

require(wesanderson)
palB <- wes_palette("FantasticFox")
pal <- colorRampPalette(rev(c(palB[5],palB[4],"white", palB[2], "yellow2")))


pdf("../output/figures/fixedchains27_heatposteriors.pdf", width=12, height=8)
par(mfrow=c(1,2), mar=c(0,0,0,0), bg="black")
#lapply(chains, function(x) plotBranchHeatMap(tree, x, "beta1", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(0.5, 1.3, 0.2)))
plotBranchHeatMap(tree, chains, "beta1", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(0.65, 0.85, 0.05), edge.width=2)
text(50, 100, labels = "Beta1", cex=2, col = "white")
#plotBranchHeatMap(tree, chains, "beta2", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(-0.005, 0.005, 0.002), edge.width=2)
#text(50, 100, labels = "Beta2", cex=2, col = "white")
plotBranchHeatMap(tree, chains, "theta", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(-3, 2.5, 0.5), edge.width=2)
text(50, 100, labels = "Theta", cex=2, col = "white")
dev.off()