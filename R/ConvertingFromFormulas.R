f <- formula(lnBMR ~ lnMass + I((lnMass)^2) + lnGS)
terms.formula(f)
model.matrix(f, dat=pred)

model <- model.BetaBMR

impute <- function(x){
  x[is.na(x)] <- pars[[paste("missing", lazy(x)$expr, sep=".")]]
  return(x)
}

shiftmap <- function(x){
  pars[[paste("beta",lazy(x)$expr, sep=".")]][betaID] * x
}

lapply(rhs, parse)
microbenchmark(model.frame(formula, data=pred))
head(cbind(pred$lnMass*pars$beta.lnMass[betaID], pred$lnMass^2*pars[['beta.I(lnMass^2)']][betaID]))

rhs <- attr(tf, "term.labels")
fla <- lapply(rhs, substitute)

make.bayouModel <- function(tree, dat, pred, formula, parameterization="OU"){
  parnames <- switch(parameterization, "OU"=c("alpha", "sig2", "k", "ntheta"), "OUrepar"=c("halflife", "Vy", "k", "ntheta"), "QG"=c("h2", "P", "Ne", "omega2", "k", "ntheta"))
  cache <- bayou:::.prepare.ou.univariate(tree, dat, SE, pred)
  tf <- terms(formula)
  rhs <- attr(tf, "term.labels")

  lik.fn <- function(pars, cache, X, model="Custom"){
    n <- cache$n
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
  }

}

model.frame.bayou <- function(formula, ...){
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", 
                 "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) 
      env <- parent.frame()
    eval(fcall, env)
  }
  else formula$model
  
}

f <- formula(lnBMR ~ shiftmap(1) + shiftmap(lnMass) + I(lnMass^2) + impute(lnGS))


theta + beta1*lnMass + beta2*I(lnMass^2) + beta3*impute(lnGS)
lnBMR ~ theta + lnMass + I(lnMass^2) + impute(lnGS)
