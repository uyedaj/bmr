## Preliminary data preparation
## Match and prepare the datasets
#td <- readRDS(paste("../output/data/", args[[1]], ".rds", sep=""))
tdgs <-  readRDS(paste("../output/data/", "tetrapods_gs", ".rds", sep=""))
td$dat$lnGenSize <- tdgs$dat$lnGenSize[match(td$phy$tip.label, tdgs$phy$tip.label)]
tree <- td$phy
dat <- td$dat
rownames(dat) <- attributes(td)$tip.label
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches
td <- make.treedata(tree, dat) # Rematch data and tree
td <- reorder(td, "postorder") # Reorder tree
td <- filter(td, !is.na(lnBMR), !is.na(lnMass))
#td_gs <- filter(td, !is.na(TempK)) # Produce a 2nd dataset of complete temperature data
tree <- td$phy
dat <- td$dat
## BM likelihood function for genome size
#gs.lik <- bm.lik(td_gs$phy, setNames(td_gs$dat[[3]], td_gs$phy$tip.label), SE=0, model="BM")
lnBMR <- setNames(dat[['lnBMR']], tree$tip.label)
lnMass <- setNames(dat[['lnMass']], tree$tip.label)
.pred <- cbind(setNames(dat[['lnMass']], tree$tip.label), setNames(dat[['lnMass']]^2, tree$tip.label),setNames(dat[['lnGenSize']], tree$tip.label))
colnames(.pred) <- c("lnMass", "lnMass2", "lnGS")

## Create a bayou cache object
cache <- bayou:::.prepare.ou.univariate(tree, setNames(dat$lnBMR, tree$tip.label), pred = .pred)
pred <- cache$pred

## Get optimized starting values and trait evolutionary models for genome size. 
require(phylolm)#$impute 
#td <- mutate(td, TempK = 1/(Temperature.C+273.15))
tdgs <- filter(td, !is.na(lnGenSize))#$impute 
fits <- lapply(c("BM"), function(x) phylolm(lnGenSize~1, data=tdgs$dat, phy=tdgs$phy, model=x))#$impute 
aics <- sapply(fits, function(x) x$aic)#$impute 
bestfit <- fits[[which(aics == min(aics))]]#$impute 

missing <- which(is.na(cache$pred[,3])) #$impute
cache$pred <- as.data.frame(cache$pred)
pv <- getPreValues(cache, 3) #$impute

lnBMR <- setNames(dat[['lnBMR']], tree$tip.label)
lnMass <- setNames(dat[['lnMass']], tree$tip.label)
.pred <- cbind(setNames(dat[['lnMass']], tree$tip.label), setNames(dat[['lnMass']]^2, tree$tip.label), setNames(dat[['lnGenSize']], tree$tip.label))
colnames(.pred) <- c("lnMass", "lnMass2", "lnGS")
endotherms <- lapply(cache$desc$tips[cache$phy$edge[c(845, 1707), 2]], function(x) cache$phy$tip.label[x])
pred <- data.frame(pred, endo=as.numeric(cache$phy$tip.label %in% unlist(endotherms)))
.pred <- pred

#pred <- select(pred, lnBMR, lnMass, lnBMRR, lnLs, growthC, humidity, elevation,  tempK, tcoldq.me, precip,pdryq.me, alt.me)

cache <- bayou:::.prepare.ou.univariate(tree, lnBMR, pred = .pred)
tmp <- lm(lnBMR ~ lnMass + lnMass2 + lnGS+endo, data=data.frame(.pred))
summary(tmp)
sumpars <- readRDS(paste("../output/data/", "sumpars_u6", ".rds", sep=""))

## Metabolic rate models
## Priors for different parameters, so that these are easily changed for all models
param.alpha <- list(scale=1)
param.sig2 <- list(scale=1)
param.beta_lnMass <- list(mean=0.7, sd=0.1)
param.beta_lnMass2 <- list(mean=0, sd=0.01)
param.beta_lnGS <- list(mean=0, sd=0.5)
param.beta_lnGSxlnMass <- list(mean=0, sd=0.25)
param.beta_endo <- list(mean=4.5, sd=0.5)
#param.beta_TempK <- list(mean=0, sd=10000)
param.k <- list(lambda=42.85, kmax=86)
param.sb=list(bmax=1, prob=1)
fixed.k <- sumpars$k
fixed.sb <- sumpars$sb
fixed.loc <- rep(0, sumpars$k)
fixed.k.nopleth <- fixed.k - 1
fixed.sb.nopleth <- fixed.sb[!(fixed.sb %in% c(350))]
fixed.loc.nopleth <- rep(0, fixed.k.nopleth)
fixed.k.nosals <- fixed.k - 3
fixed.sb.nosals <- fixed.sb[!(fixed.sb %in% c(350,394,397))]
fixed.loc.nosals <- rep(0, fixed.k.nosals)
param.theta <- list(mean=-2.5, sd=1.75)
param.pred.root <- unname(bestfit$coeff)#$impute 
param.pred.sig2 <- unname(bestfit$sigma2)#$impute 
bmax <- prob <- rep(1, nrow(cache$edge))
#bmax[which(cache$edge.length == .Machine$double.eps)] <- prob[which(cache$edge.length == .Machine$double.eps)] <- 0


## Proposal widths
D.XXX <- function(nrj) list(alpha=0.7, sig2= 0.5, beta_lnMass=0.125, beta_lnMass2=0.004, beta_lnGS=0.2, beta_lnMassxlnGS=0.05, beta_endo=0.5, k=rep(1,nrj), theta=3, slide=1, missing.pred=1)
D.1XX <- function(nrj) list(alpha=0.7, sig2= 0.5, beta_lnMass=0.1, beta_lnMass2=0.004, beta_lnGS=0.2, beta_lnMassxlnGS=0.05, beta_endo=0.5, k=rep(1,nrj), theta=0.15, slide=1, missing.pred=1)
D.XNX <- function(nrj) list(alpha=0.7, sig2= 0.5, beta_lnMass=0.2, beta_lnMass2=0.004, beta_lnGS=0.2, beta_lnMassxlnGS=0.05, beta_endo=0.5, k=rep(1,nrj), theta=1, slide=1, missing.pred=1)


## Code explanation. 5 numbers given either R - RJ; N - Fixed multiple; 1 - Fixed global; 0 - Absent
## 1: B0 (i.e. theta) 2: B1 (i.e. slope lnMass) 3: beta_lnMass2 4: beta_lnGS 5: beta_lnGSxlnMass

## Simplest regression; 11000
prior.11000 <- make.prior(tree, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                               #dbeta_TempK="dnorm",
                               dbeta_endo="dnorm",
                               dsb="fixed", dk="fixed", dtheta="dnorm"), 
                    param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                               dbeta_endo=list(mean=4.5, sd=0.5),
                               dk="fixed", dsb="fixed", 
                               dtheta=param.theta),
                    fixed=list(k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0))
)

model.11000 <- makeBayouModel(lnBMR ~ lnMass + endo, tree=cache$phy, dat=cache$dat, pred=cache$pred, rjpars = numeric(0), prior=prior.11000, impute = NULL, D=D.1XX(1))
prior.11000(model.11000$startpar)
model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik

## Quadratic regression; 11100
prior.11100 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     dbeta_lnMass2="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm"), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta),
                          fixed=list(k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0))
                          )

model.11100 <- makeBayouModel(lnBMR ~ lnMass + lnMass2 + endo, rjpars = numeric(0), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.11100, impute = NULL, D=D.1XX(1))
prior.11100(model.11100$startpar)
model.11100$model$lik.fn(model.11100$startpar, cache, cache$dat)$loglik

## Quadratic regression; mamm11100
tmp <-  getTipMap(list(k=1, ntheta=2, sb=1707, loc=0, t2=2), cache)
mammtree <- drop.tip(cache$phy, cache$phy$tip.label[tmp==1])
mammdat <- dat[tmp==2,]
mammpred <- .pred[tmp==2,]
mammcache <- bayou:::.prepare.ou.univariate(mammtree, setNames(mammdat$lnBMR, mammtree$tip.label), pred=mammpred[mammtree$tip.label,])
 
prior.mamm11100 <- make.prior(mammtree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     dbeta_lnMass2="dnorm",
                                     #dbeta_TempK="dnorm",
                                     #dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm"), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                     dbeta_lnMass=param.beta_lnMass,
                                     dbeta_lnMass2=param.beta_lnMass2,
                                     dtheta=list(mean=2.0, sd=1.75),
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta),
                          fixed=list(k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0))
)

model.mamm11100 <- makeBayouModel(lnBMR ~ lnMass + lnMass2, rjpars = numeric(0), tree=mammcache$phy, dat=mammcache$dat, pred=mammcache$pred, prior=prior.mamm11100, impute = NULL, D=D.1XX(1))
prior.mamm11100(model.mamm11100$startpar)
model.mamm11100$model$lik.fn(model.mamm11100$startpar, mammcache, mammcache$dat)$loglik

## Simple regression, mammals only; mamm11000
prior.mamm11000 <- make.prior(mammtree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                         #dbeta_TempK="dnorm",
                                         #dbeta_endo="fixed",
                                         dsb="fixed", dk="fixed", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                         #dbeta_endo="fixed",
                                         dk="fixed", dsb="fixed", 
                                         dtheta=list(mean=2.0, sd=1.75)),
                              fixed=list(k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0))
)

model.mamm11000 <- makeBayouModel(lnBMR ~ lnMass, rjpars = numeric(0), tree=mammcache$phy, dat=mammcache$dat, pred=mammcache$pred, prior=prior.mamm11000, impute = NULL, D=D.1XX(1))
prior.mamm11000(model.mamm11000$startpar)
model.mamm11000$model$lik.fn(model.mamm11000$startpar, mammcache, mammcache$dat)$loglik

## Separate intercepts + quadratic 
prior.N1100 <- make.prior(tree, plot.prior = FALSE, 
                        dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                              dbeta_lnMass2="dnorm",
                                              #dbeta_TempK="dnorm",
                                              dbeta_endo="dnorm",
                                              dsb="fixed", dk="fixed", dtheta="dnorm"), 
                                   param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                              dbeta_lnMass2=param.beta_lnMass2,
                                              dbeta_endo=param.beta_endo,
                                              dk="fixed", dsb="fixed", 
                                              dtheta=param.theta),
                                   fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc)
)


model.N1100 <- makeBayouModel(lnBMR ~ lnMass + lnMass2 + endo, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.N1100, impute = NULL, D=D.XXX(1))
prior.N1100(model.N1100$startpar)
model.N1100$model$lik.fn(model.N1100$startpar, cache, cache$dat)$loglik

## NN100
prior.NN100 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     dbeta_lnMass2="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm"), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta),
                          fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc)
)
model.NN100 <- makeBayouModel(lnBMR ~ lnMass + lnMass2 + endo, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.NN100, impute = NULL, D=D.XXX(1))
prior.NN100(model.NN100$startpar)
model.NN100$model$lik.fn(model.NN100$startpar, cache, cache$dat)$loglik


#### Genome Size Models
## Simple regression, genome size
prior.11010 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=0, sb=numeric(0), loc=numeric(0), pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.11010 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo, rjpars = numeric(0), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.11010, impute = c("lnGS"), D=D.1XX(1))
prior.11010(model.11010$startpar)
model.11010$model$lik.fn(model.11010$startpar, cache, cache$dat)$loglik

## Simple regression, genome size + interaction
prior.11011 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     dbeta_lnMassxlnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dbeta_lnMassxlnGS=param.beta_lnGSxlnMass,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=0, sb=numeric(0), loc=numeric(0), pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.11011 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo + lnGS*lnMass, rjpars = numeric(0), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.11011, impute = c("lnGS"), D=D.1XX(1))
prior.11011(model.11011$startpar)
model.11011$model$lik.fn(model.11011$startpar, cache, cache$dat)$loglik

## Separate intercepts, genome size
prior.N1010 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.N1010 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.N1010, impute = c("lnGS"), D=D.XXX(1))
prior.N1010(model.N1010$startpar)
model.N1010$model$lik.fn(model.N1010$startpar, cache, cache$dat)$loglik

## Separate intercepts, separate slopes, genome size
prior.NN010 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.NN010 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.NN010, impute = c("lnGS"), D=D.XXX(2))
prior.NN010(model.NN010$startpar)
model.NN010$model$lik.fn(model.NN010$startpar, cache, cache$dat)$loglik

## Separate intercepts, genome size + interaction
prior.N1011 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     dbeta_lnMassxlnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dbeta_lnMassxlnGS=param.beta_lnGSxlnMass,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.N1011 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo + lnGS*lnMass, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.N1011, impute = c("lnGS"), D=D.XXX(1))
prior.N1011(model.N1011$startpar)
model.N1011$model$lik.fn(model.N1011$startpar, cache, cache$dat)$loglik

## Separate intercepts, separate slopes, genome size + interaction
prior.NN011 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     dbeta_lnMassxlnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dbeta_lnMassxlnGS=param.beta_lnGSxlnMass,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.NN011 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo + lnGS*lnMass, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.NN011, impute = c("lnGS"), D=D.XXX(2))
prior.NN011(model.NN011$startpar)
model.NN011$model$lik.fn(model.NN011$startpar, cache, cache$dat)$loglik

### No Salamanders
## Separate intercepts, genome size
prior.S1010 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nosals, sb=fixed.sb.nosals, loc=fixed.loc.nosals, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.S1010 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.S1010, impute = c("lnGS"), D=D.XXX(1))
prior.S1010(model.S1010$startpar)
model.S1010$model$lik.fn(model.S1010$startpar, cache, cache$dat)$loglik

## Separate intercepts, separate slopes, genome size
prior.SS010 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nosals, sb=fixed.sb.nosals, loc=fixed.loc.nosals, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.SS010 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.SS010, impute = c("lnGS"), D=D.XXX(2))
prior.SS010(model.SS010$startpar)
model.SS010$model$lik.fn(model.SS010$startpar, cache, cache$dat)$loglik

## Separate intercepts, genome size + interaction
prior.S1011 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     dbeta_lnMassxlnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dbeta_lnMassxlnGS=param.beta_lnGSxlnMass,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nosals, sb=fixed.sb.nosals, loc=fixed.loc.nosals, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.S1011 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo + lnGS*lnMass, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.S1011, impute = c("lnGS"), D=D.XXX(1))
prior.S1011(model.S1011$startpar)
model.S1011$model$lik.fn(model.S1011$startpar, cache, cache$dat)$loglik

## Separate intercepts, separate slopes, genome size + interaction
prior.SS011 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     dbeta_lnMassxlnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dbeta_lnMassxlnGS=param.beta_lnGSxlnMass,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nosals, sb=fixed.sb.nosals, loc=fixed.loc.nosals, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.SS011 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo + lnGS*lnMass, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.SS011, impute = c("lnGS"), D=D.XXX(2))
prior.SS011(model.SS011$startpar)
model.SS011$model$lik.fn(model.SS011$startpar, cache, cache$dat)$loglik

### No Plethodon
## Separate intercepts, genome size
prior.P1010 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nopleth, sb=fixed.sb.nopleth, loc=fixed.loc.nopleth, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.P1010 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.P1010, impute = c("lnGS"), D=D.XXX(1))
prior.P1010(model.P1010$startpar)
model.P1010$model$lik.fn(model.P1010$startpar, cache, cache$dat)$loglik

## Separate intercepts, separate slopes, genome size
prior.PP010 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nopleth, sb=fixed.sb.nopleth, loc=fixed.loc.nopleth, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.PP010 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.PP010, impute = c("lnGS"), D=D.XXX(2))
prior.PP010(model.PP010$startpar)
model.PP010$model$lik.fn(model.PP010$startpar, cache, cache$dat)$loglik

## Separate intercepts, genome size + interaction
prior.P1011 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     dbeta_lnMassxlnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dbeta_lnMassxlnGS=param.beta_lnGSxlnMass,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nopleth, sb=fixed.sb.nopleth, loc=fixed.loc.nopleth, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.P1011 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo + lnGS*lnMass, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.P1011, impute = c("lnGS"), D=D.XXX(1))
prior.P1011(model.P1011$startpar)
model.P1011$model$lik.fn(model.P1011$startpar, cache, cache$dat)$loglik

## Separate intercepts, separate slopes, genome size + interaction
prior.PP011 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     dbeta_lnGS="dnorm",
                                     dbeta_lnMassxlnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm",
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dbeta_lnMassxlnGS=param.beta_lnGSxlnMass,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta,
                                     dpred.sig2="fixed", dpred.root="fixed"
                          ),
                          fixed=list(k=fixed.k.nopleth, sb=fixed.sb.nopleth, loc=fixed.loc.nopleth, pred.sig2=param.pred.sig2, pred.root=param.pred.root)
)
model.PP011 <- makeBayouModel(lnBMR ~ lnMass + lnGS + endo + lnGS*lnMass, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.PP011, impute = c("lnGS"), D=D.XXX(2))
prior.PP011(model.PP011$startpar)
model.PP011$model$lik.fn(model.PP011$startpar, cache, cache$dat)$loglik

## NN000
prior.NN000 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     #dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     #dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta
                          ),
                          fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc)
)
model.NN000 <- makeBayouModel(lnBMR ~ lnMass + endo, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.NN000, D=D.XXX(2))
prior.NN000(model.NN000$startpar)
model.NN000$model$lik.fn(model.NN000$startpar, cache, cache$dat)$loglik

## N1000
prior.N1000 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     #dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     #dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk="fixed", dsb="fixed", 
                                     dtheta=param.theta
                          ),
                          fixed=list(k=fixed.k, sb=fixed.sb, loc=fixed.loc)
)
model.N1000 <- makeBayouModel(lnBMR ~ lnMass + endo, rjpars = c("theta"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.N1000, D=D.XXX(1))
prior.N1000(model.N1000$startpar)
model.N1000$model$lik.fn(model.N1000$startpar, cache, cache$dat)$loglik

## RR000
prior.RR000 <- make.prior(tree, plot.prior = FALSE, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                     #dbeta_lnMass2="dnorm",
                                     #dbeta_lnGS="dnorm",
                                     #dbeta_TempK="dnorm",
                                     dbeta_endo="dnorm",
                                     dsb="dsb", dk="cdpois", dtheta="dnorm"
                          ), 
                          param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnMass=param.beta_lnMass,
                                     #dbeta_lnMass2=param.beta_lnMass2,
                                     #dbeta_lnGS=param.beta_lnGS,
                                     dbeta_endo=param.beta_endo,
                                     dk=param.k, dsb=param.sb, 
                                     dtheta=param.theta
                          )
)
model.RR000 <- makeBayouModel(lnBMR ~ lnMass + endo, rjpars = c("theta", "lnMass"), tree=cache$phy, dat=cache$dat, pred=cache$pred, prior=prior.RR000, D=D.XXX(2))
prior.RR000(model.RR000$startpar)
model.RR000$model$lik.fn(model.RR000$startpar, cache, cache$dat)$loglik


## N1010
## NN010
## 11011
## N1011
## NN011


priors <- lapply(objects()[grep("prior.", objects(), fixed=TRUE)], function(x) get(x)); names(priors) <- objects()[grep("prior.", objects(), fixed=TRUE)]
models <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$model); names(models) <- objects()[grep("model.", objects(), fixed=TRUE)]
startpars <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$startpar); names(startpars) <- objects()[grep("model.", objects(), fixed=TRUE)]

