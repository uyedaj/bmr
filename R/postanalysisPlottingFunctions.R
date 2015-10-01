f <- dat ~ beta1 * lnMass + endo * endo + theta*1
summaryPlots <- function(f, chain, pred, x.var="lnMass", cutoff=0.2){
  tree <- attributes(chain)$tree
  model.pars <- attributes(chain)$model.pars
  dat <- attributes(chain)$dat
  burnin <- attributes(chain)$burnin
  fc <- as.character(f)
  rhs <- fc[3]
  rhs <- gsub(" ", "", rhs)
  terms <- strsplit(rhs, "+", fixed=TRUE)[[1]]
  terms <- lapply(terms, function(x) strsplit(x, "*", fixed=TRUE)[[1]])
  coefs <- sapply(terms, function(x) x[1])
  vars <- sapply(terms, function(x) x[2])
  pp <- Lposterior(chain, reorder(tree, "postorder"), burnin)
  sumpars <- list(sb = which(pp$pp > cutoff))
  sumpars$k <- length(sumpars$sb)
  sumpars$ntheta <- length(sumpars$sb)+1
  sumpars$loc <- rep(0, sumpars$k)
  sumpars$t2 <- 2:sumpars$ntheta
  tr <- pars2simmap(sumpars, tree)
  plotSimmap(tr$tree, colors=tr$col, fsize=0.25)
  rj <- coefs %in% model.pars$rjpars
  summarizeDerivedState <- function(branch, chain, coefs, rj){
  if(sum(!rj)>0){
     globalsums <- lapply(coefs[!rj], function(x) unlist(chain[[x]]))
     names(globalsums) <- coefs[!rj]
  } else { globalsums <- list() }
    if(branch==0){
      if(sum(rj)>0 & length(coefs)>0){
        rjsums <- lapply(coefs[rj], function(x) sapply(chain[[x]], function(x) x[1]))
        names(rjsums) <- coefs[rj]
      } else {rjsums <- list()}
    } else {
      SB <- unlist(chain$sb)
      gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
      ind <- which(SB==branch)
      gen <- gen[ind]
      T2 <- unlist(chain$t2)[ind]
      rjsums <- lapply(coefs[rj], function(cf) sapply(1:length(T2), function(x) chain[[cf]][[gen[x]]][T2[x]]))
      names(rjsums) <- coefs[rj]
    }
  
    sums <- c(globalsums, rjsums)
    medians = lapply(sums, median)
    densities <- lapply(sums, density)  
  return(list(sums, medians=medians, densities=densities))
  }
  summaries <- lapply(c(0, sumpars$sb), function(x) summarizeDerivedState(x, chain, coefs, rj))
}
