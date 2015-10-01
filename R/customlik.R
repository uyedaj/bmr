custom.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  #pred[is.na(pred[,3]),3] <- pars$missing.pred
  tipreg <- rev(bayou:::.pars2map(pars,cache)$theta)
  ntipreg <- names(tipreg)
  dups <- !duplicated(ntipreg) & ntipreg %in% (1:nrow(cache$edge))[cache$externalEdge]
  tipreg <- tipreg[dups]
  ntipreg <- ntipreg[dups]
  o <- order(cache$edge[as.numeric(ntipreg), 2])
  betaID <- tipreg[o]
  #betaID <- sapply(tipreglist[cache$externalEdge][o], function(x) x[length(x)])
  #beta <- cbind(sapply(c("beta1", "beta2"), function(x) pars[[x]][betaID]), pars$beta3)
  X = X - pars$beta1[betaID]*pred[,1]
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
  #llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN)
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}