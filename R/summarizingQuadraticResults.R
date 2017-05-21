## tetrapods_ei model summaries
require(coda)
rj_modelnames <- c("ntetrapods_ei/kolokotronesNN01_chain_r1-6", "ntetrapods_ei/tetrapods_eiRR00_chain_r1-6")

quad_modelnames <- c("ntetrapods_ei/tetrapods_ei_fixed_chain_r1-4.rds", "ntetrapods_ei/tetrapods_ei_11100_r001", "ntetrapods_ei/tetrapods_ei_N1100_r001", "ntetrapods_ei/tetrapods_ei_NN100_r001",
                     "ntetrapods_ei/tetrapods_ei_mamm11100_r003",
                     "fkolokotrones/kolokotrones1111_r003", "fkolokotrones/kolokotronesN111_r003", "fkolokotrones/kolokotronesNN11_r004"  
                    )
noquad_modelnames <- c("fkolokotrones/kolokotronesNN01_r003", "ntetrapods_ei/tetrapods_ei_mamm11000_r003")

## Mammals only 

## Quadratic coefficient summary
pal <- colorRampPalette(c("#1abc9c", "#3fd47d", "#3498db", "#9b59b6", "#e74c3c"))
chains <- lapply(quad_modelnames, function(x) try(readRDS(paste("../output/runs/",x,"_chain.rds", sep=""))))
noquad_chains <-lapply(noquad_modelnames, function(x) try(readRDS(paste("../output/runs/",x,"_chain.rds", sep=""))))
names(chains[[6]])[7] <- "beta_lnMass2"
postburn <- lapply(chains, function(x) floor(0.4*length(chains[[1]]$gen)):length(x$gen))
beta_lnMass2.all <- lapply(1:length(chains), function(x) chains[[x]]$beta_lnMass2[postburn[[x]]])  
cols <- sapply(1:length(beta_lnMass2.all), function(x) makeTransparent(pal(length(beta_lnMass2.all))[x],120))[sample(1:length(beta_lnMass2.all), length(beta_lnMass2.all), replace=F)]
plot(0,0, type="n", xlim=c(-0.02, 0.02), ylim=c(0, 320))
lapply(1:length(beta_lnMass2.all), function(x) polygon(density(beta_lnMass2.all[[x]]), col=cols[x]))
leg <- unname(sapply(sapply(quad_modelnames, function(x) strsplit(x, "/", fixed=TRUE)[[1]][2]), function(x) strsplit(x, "_r", fixed=TRUE)[[1]][1]))
legend(-0.02, 320, legend=leg, fill=cols)


## Curvature table
## Dataset Model Estimate HPDinterval Marginal likelihood
summary <- list()
summary$lnMass2_effSize <- c(sapply(chains, function(x) try(effectiveSize(coda::mcmc(x$beta_lnMass2)))),NA)
summary$lnMass2_coef <- c(sapply(chains, function(x) try(median(coda::mcmc(x$beta_lnMass2)))), 0)
HPDint <- rbind(t(sapply(1:length(chains), function(x) try(HPDinterval(coda::mcmc(chains[[x]]$beta_lnMass2[postburn[[x]]]))))), HPDinterval(coda::mcmc(rnorm(1000000, 0, 0.01))))
summary$lnMass2_HPDinterval <- c(apply(t(sapply(1:length(chains), function(x) try(HPDinterval(coda::mcmc(chains[[x]]$beta_lnMass2[postburn[[x]]]))))), 1, function(x) paste(round(x, 3), collapse=",")), paste(round(HPDinterval(coda::mcmc(rnorm(1000000, 0, 0.01))),3), collapse=","))
summary <- as.data.frame(summary)
rownames(summary) <- c(quad_modelnames,"prior")

require(grid)

pdf("../output/figures/curvature_runs.pdf", height=8, width=8)
par(mar=c(5,5,2,2))
plot(0,0, type="n", xlab="", ylab=expression(paste("Curvature (", beta[lnMass^2], ")",sep="")), xaxt="n", xlim=c(0.5, 8.5), ylim=c(-0.03, 0.03), cex.lab=1.5)
abline(h=0, lty=2)
points(1:8, summary$lnMass2_coef, cex=1.5, pch=21, bg=1)
lapply(1:8, function(x){
  arrows(x, HPDint[x,1], x, HPDint[x,2], code=3, angle=90)
})

dev.off()



## Marginal likelihoods
marglnL <- rep(NA,length(quad_modelnames))
names(marglnL) <- quad_modelnames
for(i in 1:length(quad_modelnames)){
  ss <- readRDS(paste("../output/runs/", quad_modelnames[i], ".ss.rds",sep=""))
  marglnL[i] <- ss$lnr
  rm(ss)
}
marglnL




##############################################################
##############################################################
################GENOME SIZE ANALYSES##########################
##############################################################
##############################################################
gs_modelnames <- c("ntetrapods_ei/tetrapods_ei_11010_r002_chain.rds", 
                   "ntetrapods_ei/tetrapods_eiN1010_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiNN010_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiS1010_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiSS010_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiP1010_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiPP010_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_ei11011_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiN1011_r003_chain.rds",
                   "ntetrapods_ei/tetrapods_eiNN011_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiS1011_r003_chain.rds",
                   "ntetrapods_ei/tetrapods_eiSS011_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiP1011_r001_chain.rds",
                   "ntetrapods_ei/tetrapods_eiPP011_r001_chain.rds"
                )

int <- c("1", "N", "N", "N", "N", "N", "N", "1", "N", "N", "N", "N", "N", "N", "1")
slope <- c("1", "1", "N", "1", "N", "1", "N", "1", "1", "N", "1", "N", "1", "N", "1")
interaction <- c("0", "0", "0", "0", "0", "0", "0", "1", "1", "1", "1", "1", "1", "1", "1")
shifts <- c("0", "N", "N", "S", "S", "P", "P", "0", "N", "N", "S", "S", "P", "P", "1")
modelspecs <- data.frame(cbind(int, slope, interaction, shifts))

chains <- lapply(gs_modelnames, function(x) readRDS(paste("../output/runs/", x, sep="")))
postburns <- lapply(chains, function(x) floor(0.4*length(x$gen)):length(x$gen))
gsHPD <- rbind(t(sapply(1:length(chains), function(x) coda::HPDinterval(coda::mcmc(chains[[x]]$beta_lnGS[postburns[[x]]])))), coda::HPDinterval(coda::mcmc(rnorm(100000, 0, 0.5))))
gsIntHPD <-sapply(1:length(chains), function(x) try(coda::HPDinterval(coda::mcmc(chains[[x]]$beta_lnMassxlnGS[postburns[[x]]]))))
gsIntHPD[sapply(gsIntHPD, class)=="try-error"] <- cbind(NA, NA)
gsIntHPD <- rbind(do.call(rbind, gsIntHPD), coda::HPDinterval(coda::mcmc(rnorm(1000000, 0, 0.25))))
gsMeds <- c(sapply(1:length(chains), function(x) median(chains[[x]]$beta_lnGS[postburns[[x]]])),0)
gsIntMeds <- c(sapply(1:length(chains), function(x) median(chains[[x]]$beta_lnMassxlnGS[postburns[[x]]])),0)
gsIntMeds[sapply(gsIntMeds, is.null)] <- NA
gsIntMeds <- do.call(c, gsIntMeds)
sapply(1:length(chains), function(x) coda::effectiveSize(coda::mcmc(chains[[x]]$beta_lnGS[postburns[[x]]])))

pdf("../output/figures/genomesize_runs.pdf", height=8, width=8)
cols <- c("#34495e", "#e74c3c", "#c0392b")
dcols <- c(1,1, "#bdc3c7", "#2ecc71", "#379adc")
par(mar=c(5,5,2,2))
plot(0,0, type="n", xlab="", ylab=expression(paste("Effect of Genome Size (", beta[lnGS], ")",sep="")), xaxt="n", xlim=c(0.5, length(gsMeds)+0.5), ylim=c(-1, 1), cex.lab=1.5)
abline(h=0, lty=3)
lapply(1:length(gsMeds), function(x){
  arrows(x-0.1, gsHPD[x,1], x-0.1, gsHPD[x,2], code=3, angle=90, length=0.1, col = cols[as.numeric(factor(apply(modelspecs[,1:2], 1, paste, collapse="")))[x]],lwd=-1.25+as.numeric(modelspecs[x,1]) + as.numeric(modelspecs[x,2]))
})
points(1:length(gsMeds)-0.1, gsMeds, cex=1.5, pch=21, bg=dcols[as.numeric(factor(modelspecs[, 4]))])
lapply(1:length(gsMeds), function(x){
  arrows(x+0.2, gsIntHPD[x,1], x+0.2, gsIntHPD[x,2], code=3, angle=90, lwd=1, lty=2, col=cols[as.numeric(factor(apply(modelspecs[,1:2], 1, paste, collapse="")))[x]], length=0.05)
})
points(1:length(gsMeds)+0.2, gsIntMeds, cex=0.75, pch=21, bg=dcols[as.numeric(factor(modelspecs[, 4]))])
#text(4.5, -0.85, labels="No interaction")
text(4.5, -0.5, "-S")
text(6.5-0.1, -0.5, "-P")
text(9, -0.5, "+Interaction")
text(11.5, -0.5, "+Int, -S")
text(13.5, -0.5, "+Int, -P")
text(15, -1.02, "Prior")

legend(0, 1.06, legend = c("Global slope & intercept", "Separate intercepts", "Separate slopes & intercepts","No shifts", "All shifts", "- Salamander shifts", "- Plethodon shift", "Interaction coefficient"), col=c(cols, 1,1,1,1,1), lwd=c(0.75, 1.75, 2.75,0, 0,0, 0, 1), pt.lwd=c(0,0,0,1,1,1,1,0), lty=c(rep(1, 3), 0,0,0,0, 2), pch=21, pt.bg=c(0,0,0,dcols[1], dcols[3], dcols[4], dcols[5],0), pt.cex=c(0,0,0,2,2,2,2,0))
dev.off()


##########################################################################
##########################################################################
############### FULL OUTPUT PARAMETER TABLE ##############################
##########################################################################
##########################################################################
require(dplyr)
require(reshape2)
require(xtable)
rj_modelnames <- c("ntetrapods_ei/tetrapods_eiRR000_r1-6", "ntetrapods_ei/kolokotronesRR01_r1-6")

quad_modelnames <- c("ntetrapods_ei/tetrapods_eiNN000_r1-4", "ntetrapods_ei/tetrapods_ei_11100_r001", "ntetrapods_ei/tetrapods_ei_N1100_r001", "ntetrapods_ei/tetrapods_ei_NN100_r001",
                     "ntetrapods_ei/tetrapods_ei_mamm11100_r003",
                     "fkolokotrones/kolokotrones1111_r003", "fkolokotrones/kolokotronesN111_r003", "fkolokotrones/kolokotronesNN11_r004"  
)
noquad_modelnames <- c("fkolokotrones/kolokotronesNN01_r003", "ntetrapods_ei/tetrapods_ei_mamm11000_r003")

gs_modelnames <- c("ntetrapods_ei/tetrapods_ei_11010_r002", 
                   "ntetrapods_ei/tetrapods_eiN1010_r001",
                   "ntetrapods_ei/tetrapods_eiNN010_r001",
                   "ntetrapods_ei/tetrapods_eiS1010_r001",
                   "ntetrapods_ei/tetrapods_eiSS010_r001",
                   "ntetrapods_ei/tetrapods_eiP1010_r001",
                   "ntetrapods_ei/tetrapods_eiPP010_r001",
                   "ntetrapods_ei/tetrapods_ei11011_r001",
                   "ntetrapods_ei/tetrapods_eiN1011_r003",
                   "ntetrapods_ei/tetrapods_eiNN011_r001",
                   "ntetrapods_ei/tetrapods_eiS1011_r003",
                   "ntetrapods_ei/tetrapods_eiSS011_r001",
                   "ntetrapods_ei/tetrapods_eiP1011_r001",
                   "ntetrapods_ei/tetrapods_eiPP011_r001"
)

allmodels <- c("ntetrapods_ei/tetrapods_eiRR000_r1-6", 
               "ntetrapods_ei/tetrapods_eiNN000_r1-4",
               "ntetrapods_ei/tetrapods_ei11000_r001",
               "ntetrapods_ei/tetrapods_eiN1000_r002",
               "ntetrapods_ei/tetrapods_ei_11100_r001",
               "ntetrapods_ei/tetrapods_ei_N1100_r001", 
               "ntetrapods_ei/tetrapods_ei_NN100_r001",
               "ntetrapods_ei/tetrapods_ei_11010_r002", 
                "ntetrapods_ei/tetrapods_eiN1010_r001",
                "ntetrapods_ei/tetrapods_eiNN010_r001",
                "ntetrapods_ei/tetrapods_eiS1010_r001",
                "ntetrapods_ei/tetrapods_eiSS010_r001",
                "ntetrapods_ei/tetrapods_eiP1010_r001",
                "ntetrapods_ei/tetrapods_eiPP010_r001",
                "ntetrapods_ei/tetrapods_ei11011_r001",
                "ntetrapods_ei/tetrapods_eiN1011_r003",
                "ntetrapods_ei/tetrapods_eiNN011_r001",
                "ntetrapods_ei/tetrapods_eiS1011_r003",
                "ntetrapods_ei/tetrapods_eiSS011_r001",
                "ntetrapods_ei/tetrapods_eiP1011_r001",
                "ntetrapods_ei/tetrapods_eiPP011_r001", 
               "ntetrapods_ei/tetrapods_ei_mamm11000_r003",
               "ntetrapods_ei/tetrapods_ei_mamm11100_r003",
               "ntetrapods_ei/kolokotronesRR01_r1-6",
               "fkolokotrones/kolokotronesNN01_r003",
               "fkolokotrones/kolokotrones1111_r003", 
               "fkolokotrones/kolokotronesN111_r003", 
               "fkolokotrones/kolokotronesNN11_r004"
)
chains <- lapply(allmodels, function(x) try(readRDS(paste("../output/runs/",x,"_chain.rds", sep=""))))
names(chains[[1]])[c(6,13)] <- c("beta_endo", "beta_lnMass")
names(chains[[24]])[c(8, 16)] <- c("beta_TempK", "beta_lnMass")
names(chains[[2]])[c(6, 13)] <- c("beta_endo", "beta_lnMass")
lapply(chains, names)

names(chains) <- unname(sapply(sapply(allmodels, function(x) strsplit(x, "/", fixed=TRUE)[[1]][[2]]), function(y) strsplit(y, "_r", fixed=TRUE)[[1]][[1]]))
printChainsTable <- function(chains, burnin=0.4){
  summarizeChain <- function(chain, burnin=0.4){
  ngen <- length(chain$gen)
  .chain <- chain[!names(chain)%in% c("gen","sb", "loc", "t2", "missing.pred", "ntheta", "pred.sig2", "pred.root")]
  postburn <- floor(burnin*ngen):ngen
  elength <- lapply(.chain, function(x) unique(sapply(x, function(y) length(y))))
  meds <- sapply(1:length(.chain), function(x) if(elength[[x]]==1) median(.chain[[x]][postburn]) else median(sapply(.chain[[x]][postburn], function(y) y[1])))
  hpd <- t(sapply(1:length(.chain), function(x) if(elength[[x]]==1) HPDinterval(coda::mcmc(.chain[[x]][postburn])) else HPDinterval(coda::mcmc(sapply(.chain[[x]][postburn], function(y) y[1])))))
  meds <- c(meds, median(log(2)/chain$alpha[postburn]))
  hpd <- rbind(hpd, HPDinterval(coda::mcmc(log(2)/chain$alpha[postburn])))
  sum <- data.frame(parameter=c(names(.chain), "halflife"), median=meds, lowerHPD=hpd[,1], upperHPD=hpd[,2])
  #sumelt <- melt(sum)
  return(sum)
  }
  
  res <- lapply(chains, summarizeChain, burnin=burnin)
  names(res) <- names(chains)
  
  tablify <- function(results, digits=3){
  modn <- names(results)
  alltab <- melt(do.call(rbind, lapply(1:length(results), function(x) data.frame("model"=modn[x], results[[x]]))), id.vars=c("model", "parameter"))
  aggrFn <- function(x){
    xx <- c(formatC(x[1]), formatC(x[2:3], digits=3))
    if(is.na(xx[2])){
        stringx <- paste(xx[1], sep="")
      } else {
      if(xx[2]!=xx[3]){
        stringx <- paste(gsub("NA", "", as.character(xx[1])), " [", xx[2], ", ", xx[3], "]", sep="")
      } else {
        stringx <- paste(xx[1], sep="")
      }
      }
    return(stringx)
  }
  dcast(alltab, model ~ parameter, fun.aggregate=aggrFn)  
  }
  tab <- tablify(res, digits=3)
  tab[,c("model", "k", "lnL", "prior", "alpha", "halflife", "sig2", "theta", sort(colnames(tab)[!colnames(tab) %in% c("model","prior", "lnL", "k", "alpha", "halflife", "sig2", "theta")]))]
}

summaryTable <- printChainsTable(chains, burnin=0.4)
summaryTable
summaryTable[summaryTable=="NA"] <- "-"
xtable::xtable(summaryTable, caption = "Median parameter estimates for all models fit in the manuscript. Numbers in brackets given the 95\\% highest posterior density estimate.")


