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

