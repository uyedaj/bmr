getPreValues <- function(cache){
  V <- vcvPhylo(cache$phy, anc.nodes=FALSE)
  X <- cache$pred[,3]
  unknown <- is.na(X)
  known <- !unknown
  Vkk <- V[known, known]
  Vuu <- V[unknown, unknown]
  Vku <- V[known, unknown]
  Vuk <- V[unknown, known]
  iVkk <- solve(Vkk)
  sigmabar <- as.matrix(forceSymmetric(Vuu - Vuk%*%iVkk%*%Vku))
  cholSigmabar <- chol(sigmabar)
  mubarmat <- Vuk%*%iVkk
  return(list(V=V, X=X, unknown=unknown, known=known, Vkk=Vkk, Vuu=Vuu, Vku=Vku, Vuk=Vuk, iVkk=iVkk, sigmabar=sigmabar, mubarmat=mubarmat, cholSigmabar=cholSigmabar))
}

cMVNorm <- function(cache, pars, prevalues=pv, known=FALSE){
  X <- prevalues$X
  known <- prevalues$known
  unknown <- prevalues$unknown
  mu <- rep(pars$pred.root, cache$n)
  muk <- mu[known]
  muu <- mu[unknown]
  mubar <- t(muu + prevalues$mubarmat%*%(X[known]-muk))
  #sigmabar <- pars$pred.sig2*prevalues$sigmabar
  myChol <-sqrt(pars$pred.sig2)*prevalues$cholSigmabar
  res <- dmvn(pars$missing.pred, mu=mubar, sigma = myChol, log=TRUE, isChol=TRUE)
  return(res)
}

## Proposal function to simulate conditional draws from a multivariate normal distribution
.imputePredBM <- function(cache, pars, d, move,ct=NULL, prevalues=pv){
  #(tree, dat, sig2, plot=TRUE, ...){
  X <- prevalues$X
  Vuk <- pars$pred.sig2*prevalues$Vuk
  iVkk <- (1/pars$pred.sig2)*prevalues$iVkk
  Vku <- pars$pred.sig2*prevalues$Vku
  Vuu <- pars$pred.sig2*prevalues$Vuu
  known <- prevalues$known
  unknown <- prevalues$unknown
  mu <- rep(pars$pred.root, cache$n)
  muk <- mu[known]
  muu <- mu[unknown]
  mubar <- t(muu + Vuk%*%iVkk%*%(X[known]-muk))
  sigmabar <- Vuu - Vuk%*%iVkk%*%Vku
  res <- MASS::mvrnorm(1, mubar, sigmabar)
  pars.new <- pars
  pars.new$missing.pred <- res
  hr=Inf
  type="impute"
  return(list(pars=pars.new, hr=hr, decision = type))
}

make.monitorFn <- function(model, noMonitor=c("missing.pred", "ntheta"), integers=c("gen","k")){
  parorder <- model$parorder
  rjpars <- model$rjpars
  exclude <- which(parorder %in% noMonitor)
  if(length(exclude) > 0){
    pars2monitor <- parorder[-exclude]
  } else {pars2monitor <- parorder}
  if(length(rjpars) > 0){
    rjp <- which(pars2monitor %in% rjpars)
    pars2monitor[rjp] <- paste("r", pars2monitor[rjp], sep="")
  }
  pars2monitor <- c("gen", "lnL", "prior", pars2monitor)
  type <- rep(".2f", length(pars2monitor))
  type[which(pars2monitor %in% integers)] <- "i"
  string <- paste(paste("%-8", type, sep=""), collapse="")
  monitor.fn = function(i, lik, pr, pars, accept, accept.type, j){
    names <- pars2monitor
    #names <- c("gen", "lnL", "prior", "alpha" , "sig2", "rbeta1", "endo", "k")
    #string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
    acceptratios <- tapply(accept, accept.type, mean)
    names <- c(names, names(acceptratios))
    if(j==0){
      cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
    }
    cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$endo, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
  }
}

getTipMap <- function(pars, cache){
  map <- bayou:::.pars2map(pars,cache)
  tipreg <- rev(map$theta)
  ntipreg <- rev(map$branch)
  #ntipreg <- names(map$theta)
  dups <- !duplicated(ntipreg) & ntipreg %in% (1:nrow(cache$edge))[cache$externalEdge]
  tipreg <- tipreg[which(dups)]
  ntipreg <- ntipreg[which(dups)]
  o <- order(cache$edge[as.numeric(ntipreg), 2])
  betaID <- tipreg[o]
}

liks <- list(
  "fixed.shift_lnMass_lnMass2.fixed_endo_lnGS.impute_lnGS"=
    function(pars, cache, X, model="Custom"){
      n <- cache$n
      X <- cache$dat
      pred <- cache$pred
      pred[is.na(pred[,3]),3] <- pars$missing.pred #$impute
      betaID <- getTipMap(pars, cache)
      ## Specify the model here
      X = X - pars$beta1[betaID]*pred[,1] - pars$beta2[betaID]*pred[,2] - pars$beta3*pred[,3]
      cache$dat <- X
      ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
      dpars <- pars
      dpars$theta[dpars$t2[which(dpars$sb %in% c(841, 1703))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(841, 1703))]]+dpars$endo
      ### The part below mostly does not change
      X.c <- bayou:::C_weightmatrix(cache, dpars)$resid
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
      llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
      return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
    },
  "rj.shift_lnMass.fixed_endo.missing_drop"=
    function(pars, cache, X, model="Custom"){
      n <- cache$n
      X <- cache$dat
      pred <- cache$pred
      #pred[is.na(pred[,3]),3] <- pars$missing.pred #$impute
      betaID <- getTipMap(pars, cache)
      ## Specify the model here
      X = X - pars$beta1[betaID]*pred[,1]# - pars$beta2[betaID]*pred[,2] - pars$beta3*pred[,3]
      cache$dat <- X
      ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
      dpars <- pars
      dpars$theta[dpars$t2[which(dpars$sb %in% c(841, 1703))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(841, 1703))]]+dpars$endo
      ### The part below mostly does not change
      X.c <- bayou:::C_weightmatrix(cache, dpars)$resid
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
      #llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
      return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
    },
  "rj.shift_lnMass.fixed_endo.impute_TempQ10"
  )

startpars <- list(
  "fixed.shift_lnMass_lnMass2.fixed_endo_lnGS.impute_lnGS"=
    function(sb, k) {
      k <- length(sb)
      startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(k+1, 0.7, 0.1), beta2=rnorm(k+1, 0, 0.01), beta3=rnorm(1, 0, 0.05), endo=2, k=k, ntheta=k+1, theta=rnorm(k+1, 0, 1), sb=sb, loc=rep(0, k), t2=2:(k+1))
      return(startpar)
    },
  "rj.shift_lnMass.fixed_endo.missing_drop"=
    function(sb, k) {
      sb <- sample(1:length(cache$bdesc), k, replace=FALSE, prob = sapply(cache$bdesc, length))
      startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(k+1, 0.7, 0.1), endo=2, k=k, ntheta=k+1, theta=rnorm(k+1, 0, 1), sb=sb, loc=rep(0, k), t2=2:(k+1))
      return(startpar)
    }
)
monitors <- list(
  "fixed.shift_lnMass_lnMass2.fixed_endo_lnGS.impute_lnGS"=
    function(i, lik, pr, pars, accept, accept.type, j){
      names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1", "endo", "k")
      string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
      acceptratios <- tapply(accept, accept.type, mean)
      names <- c(names, names(acceptratios))
      if(j==0){
        cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
      }
      cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$endo, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
    },
  "rj.shift_lnMass.fixed_endo.missing_drop"=
    function(i, lik, pr, pars, accept, accept.type, j){
      names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1", "endo", "k")
      string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
      acceptratios <- tapply(accept, accept.type, mean)
      names <- c(names, names(acceptratios))
      if(j==0){
        cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
      }
      cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$endo, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
    }
  )

models <- list(
  "fixed.shift_lnMass_lnMass2.fixed_endo_lnGS.impute_lnGS"=
    function(startpar, monitor, lik){
      list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                        beta1=".vectorMultiplier", beta2=".vectorSlidingWindow",
                        beta3=".slidingWindowProposal", endo=".slidingWindowProposal", 
                        theta=".adjustTheta", slide=".slide", 
                        pred.sig2=".multiplierProposal" , pred.root=".slidingWindowProposal", #$impute 
                        missing.pred=".imputePredBM" #$impute
      ),
      control.weights = list(alpha=4, sig2=2, beta1=10, 
                             beta2=8, beta3=4, endo=3, 
                             theta=10, slide=2, k=0,
                             pred.sig2=1, pred.root=1, missing.pred=3 #$impute
      ),
      D = list(alpha=1, sig2= 0.75, beta1=0.75, beta2=0.05, beta3=0.05, endo=0.25, theta=2, slide=1, 
               pred.sig2=1, pred.root=1, missing.pred=1 #$impute
      ),
      parorder = names(startpar),
      rjpars = c("theta"),
      shiftpars = c("sb", "loc", "t2"),
      monitor.fn = monitor,
      lik.fn = lik)
    },
  "rj.shift_lnMass.fixed_endo.missing_drop"=
    function(startpar, monitor, lik){
      list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                        beta1=".vectorMultiplier", endo=".slidingWindowProposal", 
                        k=".splitmergebd", theta=".adjustTheta", slide=".slide"                       
      ),
      control.weights = list(alpha=4, sig2=2, beta1=5, 
                             endo=3,  k=10, theta=5, slide=2
      ),
      D = list(alpha=1, sig2= 0.75, beta1=0.75, k=c(1,1), endo=0.25, theta=2, slide=1
      ),
      parorder = names(startpar),
      rjpars = c("theta", "beta1"),
      shiftpars = c("sb", "loc", "t2"),
      monitor.fn = monitor,
      lik.fn = lik)
    }
  )


#priors <- list(
#  "fixed.shift_lnMass_lnMass2.fixed_endo_lnGS.impute_lnGS"=
#    make.prior(tree, plot.prior = FALSE, 
#               dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta1="dnorm",
#                          dbeta2="dnorm", dbeta3="dnorm", dendo="dnorm", 
#                          dsb="fixed", dk="fixed", dtheta="dnorm",
#                          dpred.sig2="dhalfcauchy", dpred.root="dnorm" #$impute
#               ), 
#               param=list(dalpha=list(scale=1), dsig2=list(scale=1), dbeta1=list(mean=0.7, sd=0.1), 
#                          dbeta2=list(mean=0, sd=0.05),dbeta3=list(mean=0, sd=0.05), 
#                          dendo=list(mean=0, sd=4),  dk="fixed", dsb="fixed", 
#                          dtheta=list(mean=0, sd=4),
#                          dpred.sig2=list(scale=1), dpred.root=list(mean=1, sd=1) #$impute
#               ), 
#               fixed=list(sb=startpar$sb, k=startpar$k)
#    ),
#  "rj.shift_lnMass.fixed_endo.missing_drop"=
#    make.prior(tree, plot.prior = FALSE, 
#               dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta1="dnorm",
#                          dendo="dnorm", dsb="dsb", dk="cdpois", dtheta="dnorm"
#                ), 
#               param=list(dalpha=list(scale=1), dsig2=list(scale=1), 
#                          dbeta1=list(mean=0.7, sd=0.1), dendo=list(mean=0, sd=4),
#                          dk=list(lambda=50, kmax=500), dsb=list(bmax=1,prob=1), 
#                          dtheta=list(mean=0, sd=2.5)
#               ), 
#               model="ffancova"
#    )
#  )


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



