## Preliminary data preparation
## Match and prepare the datasets
td <- readRDS(paste("../output/data/", args[[1]], ".rds", sep=""))
tree <- td$phy
dat <- td$dat
rownames(dat) <- attributes(td)$tip.label
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches
td <- make.treedata(tree, dat) # Rematch data and tree
td <- reorder(td, "postorder") # Reorder tree
td <- filter(td, !is.na(lnBMR.k), !is.na(lnMass.k))
#td_gs <- filter(td, !is.na(TempK)) # Produce a 2nd dataset of complete temperature data
tree <- td$phy
dat <- td$dat
## BM likelihood function for genome size
#gs.lik <- bm.lik(td_gs$phy, setNames(td_gs$dat[[3]], td_gs$phy$tip.label), SE=0, model="BM")
lnBMR <- setNames(dat[['lnBMR.k']], tree$tip.label)
lnMass <- setNames(dat[['lnMass.k']], tree$tip.label)
.pred <- cbind(setNames(dat[['lnMass.k']], tree$tip.label), setNames(dat[['lnMass.k']]^2, tree$tip.label), setNames(1/(dat[['Temperature.C']]+273.15), tree$tip.label))#,setNames(dat[['lnGenSize']], tree$tip.label))
colnames(.pred) <- c("lnMass", "lnMass2", "TempK")#, "lnGS")

## Create a bayou cache object
cache <- bayou:::.prepare.ou.univariate(tree, setNames(dat$lnBMR.k, tree$tip.label), pred = .pred)
pred <- cache$pred

## Get optimized starting values and trait evolutionary models for genome size. 
require(phylolm)#$impute 
td <- mutate(td, TempK = 1/(Temperature.C+273.15))
tdgs <- filter(td, !is.na(TempK))#$impute 
fits <- lapply(c("BM"), function(x) phylolm(TempK~1, data=tdgs$dat, phy=tdgs$phy, model=x))#$impute 
aics <- sapply(fits, function(x) x$aic)#$impute 
bestfit <- fits[[which(aics == min(aics))]]#$impute 


missing <- which(is.na(cache$pred[,3])) #$impute
pv <- getPreValues(cache) #$impute

lnBMR <- setNames(dat[['lnBMR.k']], tree$tip.label)
lnMass <- setNames(dat[['lnMass.k']], tree$tip.label)
.pred <- cbind(setNames(dat[['lnMass.k']], tree$tip.label), setNames(dat[['lnMass.k']]^2, tree$tip.label), setNames(1/(dat[['Temperature.C']]+273.15), tree$tip.label))#,setNames(dat[['lnGenSize']], tree$tip.label))
colnames(.pred) <- c("lnMass", "lnMass2", "TempK")#, "lnGS")

#pred <- select(pred, lnBMR, lnMass, lnBMRR, lnLs, growthC, humidity, elevation,  tempK, tcoldq.me, precip,pdryq.me, alt.me)

cache <- bayou:::.prepare.ou.univariate(tree, lnBMR, pred = .pred)
tmp <- lm(lnBMR ~ lnMass + lnMass2 + TempK, data=data.frame(.pred))
summary(tmp)
sumpars <- readRDS(paste("../output/data/", args[[4]], ".rds", sep=""))

## Metabolic rate models
## Priors for different parameters, so that these are easily changed for all models
param.alpha <- list(scale=1)
param.sig2 <- list(scale=1)
param.beta_lnMass <- list(mean=0.7, sd=0.3)
param.beta_lnMass2 <- list(mean=0, sd=0.01)
param.beta_TempK <- list(mean=0, sd=10000)
param.k <- list(lambda=26.15, kmax=53)
param.sb=list(bmax=1, prob=1)
fixed.k <- sumpars$k
fixed.sb <- sumpars$sb
fixed.loc <- rep(0, sumpars$k)
param.theta <- list(mean=0, sd=15)
param.pred.root <- unname(bestfit$coeff)#$impute 
param.pred.sig2 <- unname(bestfit$sigma2)#$impute 
bmax <- prob <- rep(1, nrow(cache$edge))
#bmax[which(cache$edge.length == .Machine$double.eps)] <- prob[which(cache$edge.length == .Machine$double.eps)] <- 0


## Proposal widths
D.XXX <- function(nrj) list(alpha=0.5, sig2= 0.5, beta_lnMass=0.1, beta_lnMass2=0.004, beta_TempK=75, k=rep(1,nrj), theta=2, slide=1, missing.pred=1)
D.1XX <- function(nrj) list(alpha=0.5, sig2= 0.5, beta_lnMass=0.1, beta_lnMass2=0.004, beta_TempK=75, k=rep(1,nrj), theta=0.15, slide=1, missing.pred=1)
D.XNX <- function(nrj) list(alpha=0.5, sig2= 0.5, beta_lnMass=0.2, beta_lnMass2=0.004, beta_TempK=50, k=rep(1,nrj), theta=0.7, slide=1, missing.pred=1)


## Simplest regression; 1101
prior.1101 <- make.prior(tree, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                               dbeta_TempK="dnorm",
                               dsb="fixed", dk="fixed", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                    param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                               dbeta_TempK=param.beta_TempK, dk="fixed", dsb="fixed", 
                               dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                    fixed=list(k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0), pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)

model.1101 <- makeBayouModel(lnBMR ~ lnMass + TempK, rjpars = numeric(0), cache, prior.1101, impute = "TempK", D=D.1XX(1))
prior.1101(model.1101$startpar)
model.1101$model$lik.fn(model.1101$startpar, cache, cache$dat)$loglik

## Simplest regression; N101
prior.N101 <- make.prior(tree, plot.prior = FALSE, 
                         dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                    dbeta_TempK="dnorm",
                                    dsb="fixed", dk="fixed", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                         param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                    dbeta_TempK=param.beta_TempK, dk="fixed", dsb="fixed", 
                                    dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                         fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)

model.N101 <- makeBayouModel(lnBMR ~ lnMass + TempK, rjpars = c("theta"), cache, prior.N101, impute = "TempK", D=D.XXX(1))
prior.N101(model.N101$startpar)
model.N101$model$lik.fn(model.N101$startpar, cache, cache$dat)$loglik

## 1111
prior.1111 <- make.prior(tree, plot.prior = FALSE, 
                        dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                   dbeta_lnMass2="dnorm", dbeta_TempK="dnorm",
                                   dsb="fixed", dk="fixed", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                        param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                   dbeta_lnMass2=param.beta_lnMass2,
                                   dbeta_TempK=param.beta_TempK, dk="fixed", dsb="fixed", 
                                   dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                        fixed=list(k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0), pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)


model.1111 <- makeBayouModel(lnBMR ~ lnMass + lnMass2 + TempK, rjpars = numeric(0), cache, prior.1111, impute = "TempK", D=D.1XX(1))
prior.1111(model.1111$startpar)
model.1111$model$lik.fn(model.1111$startpar, cache, cache$dat)$loglik

## NN01
prior.NN01 <- make.prior(tree, plot.prior = FALSE, 
                        dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                   #dbeta_lnMass2="dnorm", 
                                   dbeta_TempK="dnorm",
                                   dsb="fixed", dk="fixed", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                        param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                   #dbeta_lnMass2=param.beta_lnMass2,
                                   dbeta_TempK=param.beta_TempK, dk="fixed", dsb="fixed", 
                                   dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                        fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)

model.NN01 <- makeBayouModel(lnBMR ~ lnMass + TempK, rjpars = c("theta", "lnMass"), cache, prior.NN01, impute = "TempK", D=D.XNX(2))
prior.NN01(model.NN01$startpar)
model.NN01$model$lik.fn(model.NN01$startpar, cache, cache$dat)$loglik

## N111
prior.N111 <- make.prior(tree, plot.prior = FALSE, 
                         dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                    dbeta_lnMass2="dnorm", 
                                    dbeta_TempK="dnorm",
                                    dsb="fixed", dk="fixed", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                         param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                    dbeta_lnMass2=param.beta_lnMass2,
                                    dbeta_TempK=param.beta_TempK, dk="fixed", dsb="fixed", 
                                    dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                         fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)

model.N111 <- makeBayouModel(lnBMR ~ lnMass + lnMass2 + TempK, rjpars = c("theta"), cache, prior.N111, impute = "TempK", D=D.XXX(1))
prior.N111(model.N111$startpar)
model.N111$model$lik.fn(model.N111$startpar, cache, cache$dat)$loglik


## NN11
prior.NN11 <- make.prior(tree, plot.prior = FALSE, 
                        dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                   dbeta_lnMass2="dnorm", 
                                   dbeta_TempK="dnorm",
                                   dsb="fixed", dk="fixed", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                        param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                   dbeta_lnMass2=param.beta_lnMass2,
                                   dbeta_TempK=param.beta_TempK, dk="fixed", dsb="fixed", 
                                   dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                        fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)

model.NN11 <- makeBayouModel(lnBMR ~ lnMass + lnMass2 + TempK, rjpars = c("theta", "lnMass"), cache, prior.NN11, impute = "TempK", D=D.XNX(2))
prior.NN11(model.NN11$startpar)
model.NN11$model$lik.fn(model.NN11$startpar, cache, cache$dat)$loglik

priors <- lapply(objects()[grep("prior.", objects(), fixed=TRUE)], function(x) get(x)); names(priors) <- objects()[grep("prior.", objects(), fixed=TRUE)]
models <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$model); names(models) <- objects()[grep("model.", objects(), fixed=TRUE)]
startpars <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$startpar); names(startpars) <- objects()[grep("model.", objects(), fixed=TRUE)]

## RR01
prior.RR01 <- make.prior(tree, plot.prior = FALSE, 
                         dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                    #dbeta_lnMass2="dnorm", 
                                    dbeta_TempK="dnorm",
                                    dsb="dsb", dk="cdpois", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                         param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                    #dbeta_lnMass2=param.beta_lnMass2,
                                    dbeta_TempK=param.beta_TempK, dk=param.k, dsb=list(bmax=bmax, prob=prob), 
                                    dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                         fixed=list(pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)

model.RR01 <- makeBayouModel(lnBMR ~ lnMass + TempK, rjpars = c("theta", "lnMass"), cache, prior.RR01, impute = "TempK", D=D.XNX(2))
prior.RR01(model.RR01$startpar)
model.RR01$model$lik.fn(model.RR01$startpar, cache, cache$dat)$loglik

## RR11
prior.RR11 <- make.prior(tree, plot.prior = FALSE, 
                         dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                    dbeta_lnMass2="dnorm", 
                                    dbeta_TempK="dnorm",
                                    dsb="dsb", dk="cdpois", dtheta="dnorm", dpred.sig2="fixed", dpred.root="fixed"), 
                         param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                    dbeta_lnMass2=param.beta_lnMass2,
                                    dbeta_TempK=param.beta_TempK, dk=param.k, dsb=list(bmax=bmax, prob=prob), 
                                    dtheta=param.theta, dpred.sig2="fixed", dpred.root="fixed"),
                         fixed=list(pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)

model.RR11 <- makeBayouModel(lnBMR ~ lnMass + lnMass2 + TempK, rjpars = c("theta", "lnMass"), cache, prior.RR11, impute = "TempK", D=D.XNX(2))
prior.RR11(model.RR11$startpar)
model.RR11$model$lik.fn(model.RR11$startpar, cache, cache$dat)$loglik


priors <- lapply(objects()[grep("prior.", objects(), fixed=TRUE)], function(x) get(x)); names(priors) <- objects()[grep("prior.", objects(), fixed=TRUE)]
models <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$model); names(models) <- objects()[grep("model.", objects(), fixed=TRUE)]
startpars <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$startpar); names(startpars) <- objects()[grep("model.", objects(), fixed=TRUE)]

