## Stitch together a full tree from all taxa
require(phytools)
require(geiger)
require(aRbor)
require(rotl)
setwd("~/repos/bmr/R/")

vertTree <- read.tree("../datasets/tetrapods.tre")

allwhite <- read.csv("../datasets/all_white_etal.csv")
k <- 8.617332478E-5
T2 <- 293.15

tmp <- filter_(allwhite, paste("!is.na(mass)", sep=""), paste("!is.na(temp)", sep=""), paste("!is.na(MR)", sep=""))
lnBMR <- log(tmp$MR)
lnMass <- log(tmp$mass)
temp <- 1/k*1/(tmp$temp+273.15)
mod1 <- lm(lnBMR~lnMass+temp+tmp$species)
Ei <- mod1$coef[3]
uspecies <- unique(allwhite$species)
ndat <- lapply(1:length(uspecies), function(x){tmp <- filter_(allwhite, paste("species=='", uspecies[x], "'", sep=""),
                                                      paste("!is.na(temp)", sep=""), paste("!is.na(MR)", sep=""));
                                       data.frame("species"=uspecies[x], "lnMass"=mean(log(tmp$mass)), 
                                                 "lnBMR"= mean(log(tmp$MR))+Ei/k*(mean(1/(tmp$temp+273.15))-1/T2))})
ndat <- do.call(rbind, ndat)
rownames(ndat) <- NULL
## Some data don't have raw masses, only mean mass, take the log of the means as an estimate of the mean of the logs...
ndat[which(is.na(ndat[,2])),2] <- sapply(ndat[which(is.na(ndat[,2])),1], function(x) log(mean(allwhite$mean.mass[which(allwhite$species==x)], na.rm=TRUE)))

ttol_queryOttIds <- function(ttol, nrequest=500, tree=TRUE){
  if(tree){
    ttol_names <- ttol$tip.label
  } else {
    ttol_names <- ttol
  }
  ttol_names <- gsub("_", " ", ttol_names) %>% gsub("'", "", .)
  ttol_names <- lapply(as.character(ttol_names), function(x) strsplit(x, " ", fixed=TRUE)[[1]][1:2]) %>% sapply(., function(x) paste(x, collapse=" "))
  seq1 <- c(seq(1, length(ttol_names), nrequest), length(ttol_names)+1)
  cat("Searching for exact matches....\n")
  tol_otl <- lapply(1:(length(seq1)-1), function(x) cbind(seq1[x]:(seq1[x+1]-1),tnrs_match_names(ttol_names[seq1[x]:(seq1[x+1]-1)], do_approximate_matching = FALSE)))
  tol_otl <- do.call(rbind, tol_otl)
  colnames(tol_otl) <- c("No", "search_string", "unique_name", "approximate_match", "ott_id", "number_matches", "is_synonym", "is_deprecated")
  not_found <- which(is.na(tol_otl[,2]))
  cat(paste(length(ttol_names) - length(not_found), " exact matches found in OTL database\n", sep=""))
  if(length(not_found) > 0){
    cat(paste("Performing approximate search for ", length(not_found), " unmatched taxa\n", sep=""))
    seq2 <- c(seq(1, length(not_found), nrequest), length(not_found)+1)
    missing <- lapply(1:(length(seq2)-1), function(x) cbind(not_found[seq2[x]:(seq2[x+1]-1)],tnrs_match_names(ttol_names[not_found[seq2[x]:(seq2[x+1]-1)]], do_approximate_matching = TRUE)))
    missing <- do.call(rbind, missing)
    tol_otl[not_found[which(!is.na(missing[,2]))], ] <- missing[which(!is.na(missing[,2])),]
    cat(paste("Found ", length(which(!is.na(missing[,2]))), " approximate matches\n", sep=""))
    cant_find <- not_found[which(is.na(missing[,2]))]
    if(length(cant_find>0) & tree){
      new.ttol <- drop.tip(ttol, cant_find)   
      tol_otl <- tol_otl[-cant_find,]
    } else {
      new.ttol <- gsub(" ", "_", tol_otl$unique_name)
    }
    dropped=ttol_names[cant_find]
  } else {
    dropped=NULL
    new.ttol <- gsub(" ", "_", tol_otl$unique_name)
  }
  return(list(ttol=new.ttol, tol_otl=tol_otl, dropped=dropped))
}

missing <- setdiff(ndat[,1], vertTree$tip.label)
otl_missing <- ttol_queryOttIds(missing, nrequest=500, tree=FALSE)
clean_dat <- ndat
clean_dat[,1] <- as.character(clean_dat[,1])
clean_dat[match(missing, clean_dat[,1]),1] <- otl_missing$ttol

td <- make.treedata(vertTree, clean_dat)
patches <- read.csv("../datasets/BMR_corrected.csv")
patches2 <- read.csv("../datasets/BMR_corrected2.csv")
td$dat[match(patches2$taxon, td$phy$tip.label),] <- data.frame("lnMass"=log(patches2$newMass.g), "lnBMR"=patches2$NewRate.20)
drop <- patches$taxon[which(patches$action %in% c("cut"))]
td <- treeply(td, drop.tip, as.character(drop))


saveRDS(td, "../output/data/tetrapods_ei.rds")

