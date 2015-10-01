## # Figures for Evolution 2015 in Brazil
## Load tree and data into R
require(ape) 
require(aRbor)
require(devtools)
#load_all("~/repos/bayou/bayou_1.0")
#install_github("uyedaj/bayou", ref="dev")
require(bayou)
setwd("~/repos/bmr/R")

td <- readRDS("../output/data/tetrapods_ei.rds")
lnMass <- td$dat$lnMass
lnBMR <- td$dat$lnBMR
rr <- c("_t001", "_t002", "_t003", "_t004", "_t005", "_t006")

chain1 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[1], ".chain.rds", sep=""))
chain2 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[2], ".chain.rds", sep=""))
chain3 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[3], ".chain.rds", sep=""))
chain4 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[4], ".chain.rds", sep=""))
chain5 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[5], ".chain.rds", sep=""))
chain6 <- readRDS(paste("../output/runs/tetrapods_ei/tetrapods_ei_", rr[6], ".chain.rds", sep=""))
mymcmc_gs <- readRDS("../output/runs/tetrapods_ei/tetrapods_ei__t001.mcmc.rds")
tree <- mymcmc_gs$tree
tree <- reorder(tree, "postorder")
dat <- mymcmc_gs$dat
pred <- mymcmc_gs$pred
chains <- list(chain1, chain2, chain3, chain4, chain5, chain6)
Lpps <- lapply(chains, function(x) Lposterior(x, tree, burnin=0.3))
pairs(sapply(Lpps, function(x) x$pp))

chain_gs1 <- combine.chains(chain1, chain2, burnin.prop=0.3)
chain_gs2 <- combine.chains(chain3, chain4, burnin.prop=0.3)
chain_gs3 <- combine.chains(chain5, chain6, burnin.prop=0.3)
chain <- combine.chains(chain_gs1, chain_gs2, burnin.prop=0)
chain <- lapply(1:length(chain), function(x) c(chain[[x]], chain_gs3[[x]]))
names(chain) <- names(chain_gs1)
#chain_gs <- combine.chains(chain_gs, chain_gs3, burnin.prop=0)
thinseq <- floor(seq(1, length(chain$gen), length.out=10000))
chain$gen <- 1:length(chain$gen)
attributes(chain)$model <- attributes(chain1)$model; attributes(chain)$model.pars <- attributes(chain1)$model.pars; attributes(chain)$tree <- attributes(chain1)$tree; attributes(chain)$dat <- attributes(chain1)$dat; attributes(chain)$class <- attributes(chain1)$class; attributes(chain)$burnin <-0
plot(chain)
sumstats <- summary(chain)


pt.col <- "#2ecc71" #"#3498db"#"#1abc9c"#"#27ae60""#ecf0f1
axis.col <- "#ecf0f1"#"#2980b9"##"#3498db"
lwd=5
phy.lwd <- 1.5


pdf("../output/figures/allBMR.pdf")
par(bg="black")
plot(lnMass, lnBMR, pch=21, cex=1, bg = makeTransparent(pt.col,150), col= makeTransparent(pt.col,250))
axis(side = 1, lwd = lwd, col=axis.col, col.axis=axis.col, cex=2, cex.axis=1.5)
axis(side = 2, lwd = lwd,col=axis.col, col.axis=axis.col, cex=2, cex.axis=1.5)
plot(tree, edge.color=pt.col, show.tip.label=FALSE, edge.width=phy.lwd, direction="leftwards")
dev.off()

cache <- bayou:::.prepare.ou.univariate(tree, dat, SE=0, pred)
pdf("../output/figures/rjallometry.pdf")
tipregs <- bayou:::.tipregime(list(k=2, ntheta=3, sb=c(845,1707), loc=c(0,0), t2=2:3), tree)
endotherms <- unlist(lapply(2:3, function(x) names(tipregs[tipregs==x])))
for(i in seq(20000,40000, 200)){
  par(mfrow=c(1,2), bg="black",mar=c( 5.1, 4.1, 4.1, 2.1))
  pars <- pull.pars(i, chain, model=mymcmc_gs$model.pars)
  tr <- pars2simmap(pars, tree)
  tr$col[1] <- pt.col
  sb <- pars$sb
  tipregs <- bayou:::.tipregime(pars, tree)
  descendents <- lapply(1:(length(sb)+1), function(x) names(tipregs[tipregs==x])) 
  
  regsize <- sapply(descendents, length) 
  thermy <- sapply(1:length(descendents), function(x) sum(descendents[[x]] %in% endotherms)/regsize[x])
  #regsize <- c(length(tree$tip.label)-sum(regsize), regsize)
  ptcols <- tr$col[tipregs]
  plot(lnMass, lnBMR, pch=21, cex=0.3, bg = makeTransparent(ptcols,150), col= makeTransparent(ptcols,250), ylim=c(-6,12), xlim=c(-2, 12))
  axis(side = 1, lwd = lwd, col=axis.col, col.axis=axis.col, cex=2, cex.axis=1.5)
  axis(side = 2, lwd = lwd,col=axis.col, col.axis=axis.col, cex=2, cex.axis=1.5)
  lapply(1:length(tr$col), function(x) abline(pars$theta[x]+pars$endo*round(thermy[x],0), pars$beta1[x], lwd=min((regsize[x]+10)/100*3, 3), col=makeTransparent(tr$col[x],min((regsize[x]+10)/100*255, 255))))
  plotSimmap(tr$tree, col=tr$col, ftype="off", direction="leftwards")  
  print(i)
}
dev.off()

