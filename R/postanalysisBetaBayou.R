require(bayou)
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
  plotSimmap.mcmc(tr$tree, pbchain, burnin=burnin, fsize=0.2)
  plot(pred[,1], dat, type="n", main=i, 
       xlab="Mass", ylab="BMR", 
       xlim=c(min(pred), max(pred))*c(0.9,1.1), 
       ylim=c(min(dat), max(dat))*c(0.9,1.1))
  for(i in 1:est.pars$ntheta){
    points(pred[tipreg==i,],dat[tipreg==i], pch=21, 
           bg=bayou:::makeTransparent(i,50), col=bayou:::makeTransparent(i,50))
    abline(a=est.pars$theta[i], b=est.pars$beta1[i], lwd=2, col=i,lty=2)
  }
}

taxa <- list.files("../output/runs")

for(i in taxa){
  chain <- readRDS(paste("../output/runs/",i, "/",i,"chain.rds", sep=""))
  mymcmc <- readRDS(paste("../output/runs/",i, "/",i,"mcmc.rds", sep=""))
  tree <- mymcmc$tree
  dat <- mymcmc$dat
  pred <- mymcmc$pred
  dir <- paste("../output/figures/", sep="")
  #dir.create(dir)
  pdf(paste(dir ,i,"posteriors.pdf", sep=""))
  try(plotPosteriors(tree, dat, pred, chain, burnin=0.3, cutoff=0.1))
  dev.off()
}


