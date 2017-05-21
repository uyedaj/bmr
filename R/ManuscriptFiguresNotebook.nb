---
title: "Final Figures Notebook"
output: html_notebook
---

Creating final figures for the manuscript.

```{r}
## Load in fixed chains
require(ape) 
require(aRbor)
require(devtools)
#load_all("~/repos/bayou/bayou_1.0")
#install_github("uyedaj/bayou", ref="dev")
require(bayou)
setwd("~/repos/bmr/R")

rr <- c("_fixed_r001", "_fixed_r002", "_fixed_r003", "_fixed_r004")

chains <- lapply(rr, function(x) readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", x, ".chain.rds", sep="")))

chain <- combine.chains(chains, burnin.prop=0.3)

mymcmc_gs <- readRDS("../output/runs/tetrapods_ei/tetrapods_ei__t001.mcmc.rds")

tree <- mymcmc_gs$tree
tree <- reorder(tree, "postorder")
dat <- mymcmc_gs$dat
pred <- mymcmc_gs$pred

cutoff <- 0.2
sumpars <- list(sb = which(sumstats$branch.posteriors$pp > cutoff))
sumpars$k <- length(sumpars$sb)
sumpars$ntheta <- length(sumpars$sb)+1
sumpars$loc <- rep(0, sumpars$k)
sumpars$t2 <- 2:sumpars$ntheta
tr <- pars2simmap(sumpars, tree)
descendents <- lapply(1:(length(sumpars$sb)+1), function(x) names(tipregs[tipregs==x])) 

```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.
```{r}
## Figure 2 boxplots
betas <- do.call(rbind, chain$beta1)
thetas <- do.call(rbind, chain$theta)
endo <- chain$endo
ids <- which(sapply(descendents, length)>3)
cladenames <- c("Root", "Pleuronectinae", "Plethodontidae", "Salamandroidea", "Caudata", "Serpentes", "Squamata", "Chiroptera", "Mammalia", "Aves")
o <- c(1, 4, 5, 9, 10, 2, 3, 8, 6, 7)
ids <- ids[o]
cladenames <- cladenames[o]
betas <- betas[,ids]
thetas <- thetas[,ids]
thetas[,c(4,5,8)] <- thetas[,c(4,5,8)] + cbind(endo, endo, endo)
ldesc <- sapply(ids, function(x) length(descendents[[x]]))


pal <- colorRampPalette(c("#ecf0f1", "#2ecc71", "#3498db", "#9b59b6","#e74c3c",  "#e67e22", "#f1c40f"))
set.seed(5)
cols <- setNames(pal(length(sb)+1), c(1, sample(2:(length(sb)+1), length(sb), replace=FALSE)))
cols <- cols[order(as.numeric(names(cols)))]
#plotSimmap(tr$tree, colors=cols, ftype="off", mar=c(0.1,0.1,0.1,6))

require(vioplot)
pdf("../output/figures/Fig2_boxplots.pdf", height=10, width=8, useDingbats = FALSE)
par(mfrow=c(2,1), lwd=2)
par(mar=c(1, 5, 8, 1))
plot(0,0, type="n", xlab="", ylab="Intercept", xlim=c(0,11), ylim=c(-5,4.5), xaxt="n", cex.lab=1.25, yaxt="n")
axis(2, at=seq(-6,6,2), cex=1.25)
abline(h=seq(-6,5,1), col="gray50")
lapply(1:length(ids), function(x) vioplot(thetas[,x], at=x, add=TRUE, col=cols[ids[x]], border="black"))
box()

par(mar=c(9, 5, 0, 1))
plot(0,0, type="n", xlab="", ylab="Slope", xlim=c(0,11), ylim=c(0.3,1.05), xaxt="n", cex.lab=1.25, yaxt="n")
abline(h=seq(0,2,0.1), col="gray50")
axis(1, at=1:10, paste(cladenames, " (", ldesc,")", sep=""), las=2, cex=1.5)
axis(2, at=seq(0,1.1,0.3),  cex=1.25)
abline(h=c(2/3, 0.75), lty=2, lwd=3)
lapply(1:length(ids), function(x) vioplot(betas[,x], at=x, add=TRUE, col=cols[ids[x]], cex=2))
box()
dev.off()



```

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).
