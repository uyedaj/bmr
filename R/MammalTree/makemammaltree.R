## # Building a mammal tree
## Load tree and data into R
require(ape)
require(aRbor)
require(devtools)
require(bayou)
setwd("~/repos/bmr/R/MammalTree")
taxa <- "mammals"

td <- readRDS(paste("../../output/data/", taxa, ".rds", sep=""))
td <- make.treedata(td$phy, td$dat, name_column=1)
tree <- td$phy
dat <- td$dat

taxalist <- rownames(dat)

## # Querying OpenTree
require(devtools)
#install_github("rotl", user="fmichonneau")
require(rncl)
require(rotl)
require(rjson)
require(geiger)
require(foreach)
require(doParallel)
source("./otTimetreeFunctions.R")

registerDoParallel(cores=12)
tax <- foreach(i=1:length(taxalist)) %dopar% try(rotl::tnrs_match_names(taxalist[i], context="Mammals"))
tax <- do.call(rbind, tax)

ottids <- as.character(tax$ott_id) 

ottids

synth_tree <- reorder(rotl::tol_induced_subtree(ott_ids=ottids), "postorder")
pdf("mammalsynthtree.pdf", height=15)
plot(synth_tree, cex=0.25)
plot(tree, cex=0.25)
dev.off()

## Find all source trees with our focal taxa in them, filter out trees without branch lengths. 
studyTrees <- findsourcetrees(ottids)      
possibleTrees <- getStudyTrees(studyTrees)
possibleTrees$summary

possibleTrees$summary[8,]

timetrees <- possibleTrees$result[which(possibleTrees$summary$timeUnits=="Myr")]
timetrees <- lapply(timetrees, function(x) x[[1]])
timetrees_treeId <- rownames(possibleTrees$summary[which(possibleTrees$summary$timeUnits=="Myr"),])
timetrees_studyId <- unlist(studyTrees[match(timetrees_treeId, studyTrees[,2]),1])
studyMeta <- lapply(timetrees_studyId, function(x) rotl::get_study_meta(study_id = x))
citations <- lapply(studyMeta, function(x) x$nexml[[2]])
print(citations)


problemtree <- list(studyId="pg_2688", treeId="tree6240")
ptree <- rotl::get_study_tree(problemtree$studyId, problemtree$treeId, text_format="newick", file="probtree.tre")
tree_string <- ptree
tree_string <- gsub( "\\[.*?\\]", "", tree_string)
tree_string <- gsub(" ", "_", tree_string)
tree_string <- gsub("'", "", tree_string, fixed=TRUE)
ptree <- read.newick(text=tree_string)
plot(ptree)
timetrees[[1]] <- ptree

cals <- synthCalibrations(synth_tree, timetrees)
ncals <- lapply(cals, nrow)
cals <- do.call(rbind, cals)
BE_cals <- synthCalibrations(synth_tree, list(tree))[[1]]
calid <- sapply(1:nrow(cals), function(x){ paste(sort(cals[x,4:5]),collapse="")})

## Look for duplicate calibrations and span the max and min of that calibration for the age. 
calibrations <- NULL
for(i in 1:length(unique(calid))){
  tmp <- subset(cals, calid==unique(calid)[i])
  if(nrow(tmp)>1){
    minage <- min(tmp[,2])
    maxage <- max(tmp[,3])
    tmp[1,2:3] <- c(maxage, minage)
    rowi <- tmp[1,]
  } else {
    rowi <- tmp
  }
  calibrations <- rbind(calibrations, rowi)
}

## Pathd8 first with BE
target <- synth_tree
target$tip.label <- unname(sapply(synth_tree$tip.label,function(x) strsplit(x, "_ott")[[1]][1]))
target$edge.length <- rep(1,nrow(target$edge))
target <- chronopl(target, lambda=1)
tr <- replaceHashes(target, BE_cals)
phypd8 <- PATHd8.phylo(tr$scion, tr$table, base = ".tmp_PATHd8", rm = FALSE)

##Pathd8 again with better cals:
#phypd8$node.label <- target$node.label
## A set of hand-entered calibrations from Meredith 2011 supplement:
#meredithCalibrations <- list(
#                              data.frame("MRCA"=NA, "MaxAge"=238.16, "MinAge"=210.71, "taxonA"="Ornithorhynchus_anatinus", "taxonB"="Hipposideros_galeritus"),#All mammals
#                              data.frame("MRCA"=NA, "MaxAge"=116.81, "MinAge"=92.12, "taxonA"="Cercopithecus_mitis", "taxonB"="Cyclopes_didactylus"),#Placentalia
#                              data.frame("MRCA"=NA, "MaxAge"=82.64, "MinAge"=79.07, "taxonA"="Elephantulus_rufescens", "taxonB"="Procavia_capensis"),#Afrotheria
#                              data.frame("MRCA"=NA, "MaxAge"=68.42, "MinAge"=62.66, "taxonA"="Chaetophractus_vellerosus", "taxonB"="Cyclopes_didactylus"),#Xenarthra
#                              data.frame("MRCA"=NA, "MaxAge"=67.39, "MinAge"=65.22, "taxonA"="Paranyctimene_raptor", "taxonB"="Molossus_molossus"),#Chiroptera
#                              data.frame("MRCA"=NA, "MaxAge"=82.78, "MinAge"=77.46, "taxonA"="Panthera_leo", "taxonB"="Manis_javanica"),#Carnivora+Philodota
#                              data.frame("MRCA"=NA, "MaxAge"=72.95, "MinAge"=66.03, "taxonA"="Thomomys_bottae", "taxonB"="Glis_glis"),#Rodentia
#                              data.frame("MRCA"=NA, "MaxAge"=87.00, "MinAge"=75.40, "taxonA"="Cercopithecus_mitis", "taxonB"="Glis_glis"),#Glires
#                              data.frame("MRCA"=NA, "MaxAge"=90.49, "MinAge"=78.11, "taxonA"="Cercopithecus_mitis", "taxonB"="Galago_senegalensis"),#Primates+Dermoptera
#                              data.frame("MRCA"=NA, "MaxAge"=90.27, "MinAge"=71.68, "taxonA"="Didelphis_virginiana", "taxonB"="Macrotis_lagotis"),#Marsupialia
##                              data.frame("MRCA"=NA, "MaxAge"=49.57, "MinAge"=30.36, "taxonA"="Ornithorhynchus_anatinus", "taxonB"="Tachyglossus_aculeatus")#Monotremata
#                              )
meredith <- read.tree("../../datasets/meredith-mammals.tre")
meredithotts <- rotl::tnrs_match_names(meredith$tip.label, context_name = "Mammals")$ott_id
meredithtax <- gbresolve(finalTree$tip.label)

meredithCalibrations <- congruify.phylo(meredith, finalTree, meredithtax, tol=0, scale=NA)$calibrations

calibrations <- rbind(calibrations, meredithCalibrations)
tr <- replaceHashes(phypd8, calibrations)
#excluded <- c(144,146,147,148,149, 151,152,153,154,171)
rcalibrations <- reconcileCalibrations(tr$scion, tr$table)

finalTree <- PATHd8.phylo(tr$scion, rcalibrations, base = ".tmp_PATHd8", rm = FALSE)
finalTree <- reorder(finalTree, "postorder")

BETree <- reorder(tr$scion, "postorder")
#plotTree(BETree, colors='red', ftype="off")

#plot(c(0,0), c(0, 0), type="n", xlim=c(0, 170), ylim=c(0, 6), xaxt="n", yaxt="n",xlab="", ylab="", bty="n")
rsBETree <- BETree
rsBETree$edge.length <- rsBETree$edge.length/max(branching.times(rsBETree))

rsfTree <- finalTree
rsfTree$edge.length <- rsfTree$edge.length/max(branching.times(rsfTree))

pdf("trees.pdf", height=25)
#plot(BETree, cex=0.25)
#abline(v=c(166.2-seq(0, 200, 50)), col="red", lty=2)
#abline(v=166.2-65, col="blue", lty=2)
#plot(finalTree, cex=0.25)
nh2 <- max(nodeHeights(finalTree))
#abline(v=c(nh2-seq(0, 200, 50)), col="red", lty=2)
#abline(v=nh2-65, col="blue", lty=2)
plotTree(rsBETree, color="red", ftype="off")
plotTree(rsfTree, add=TRUE, ftype="off",lty=2, ftype="off")
plot(BETree,edge.color="black", cex=0.25, x.lim=c(-45,180))
abline(v=c(166.2-seq(0, 200, 50)), col="red", lty=2)
abline(v=166.2-65, col="blue", lty=2)
plot(finalTree,edge.color="black", cex=0.25, x.lim=c(0,225))
abline(v=c(nh2-seq(0, 200, 50)), col="red", lty=2)
abline(v=nh2-65, col="blue", lty=2)
dev.off()
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

## Have to match tree back with data....
ottids
tipotts <- list(sapply(synth_tree$tip.label, function(x) strsplit(x,"_ott")[[1]][2]), sapply(synth_tree$tip.label, function(x) strsplit(x,"_ott")[[1]][1]))
datnames <- unname(sapply(tax[match(tipotts[[1]][match(finalTree$tip.label, tipotts[[2]])], ottids),1], simpleCap))
newdat <- dat[datnames,]
rownames(newdat) <- finalTree$tip.label

td <- make.treedata(finalTree, newdat)
saveRDS(td, file = "../../output/data/newmammals.rds")

pdf("newtrees.pdf", height=20)
plot(finalTree, cex=0.4)
nn <- sapply(1:nrow(BE_cals), function(x) getMRCA(target, as.character(BE_cals[x,4:5])))
nodelabels(node = nn, pch=21, bg="green", cex=0.8)

nn <- sapply(1:nrow(calibrations), function(x) getMRCA(target, as.character(calibrations[x,4:5])))
nodelabels(node = nn, pch=21, bg="red", cex=0.5)
plotTree(reorder(tree,"cladewise"), lwd=2, xlim=c(0,250), fsize=0.3)
abline(v=c(166.2-seq(0, 200, 50)), col="red", lty=2)
abline(v=166.2-65, col="blue", lty=2)
plotTree(reorder(BETree, "cladewise"), lwd=1, color="red", xlim=c(0,250), fsize=0.3)
abline(v=c(166.2-seq(0, 200, 50)), col="red", lty=2)
abline(v=166.2-65, col="blue", lty=2)
plotTree(reorder(finalTree,"cladewise"), lwd=1, color="green", xlim=c(0,250), fsize=0.3)
abline(v=c(max(nodeHeights(finalTree))-seq(0, 200, 50)), col="red", lty=2)
abline(v=max(nodeHeights(finalTree))-65, col="blue", lty=2)
plotTree(reorder(meredith, "cladewise"), lwd=1, color="blue", xlim=c(0,250), fsize=0.3)
abline(v=c(max(nodeHeights(meredith))-seq(0, 200, 50)), col="red", lty=2)
abline(v=max(nodeHeights(meredith))-65, col="blue", lty=2)
dev.off()

plot(cophenetic(tree)[lower.tri(cophenetic(tree))]/2,xlim=c(0,100000), ylim=c(0,250))
points(cophenetic(finalTree)[lower.tri(cophenetic(finalTree))]/2, col="red")
tmp <- list("BE"=tree, "otBE"=BETree, "Meredith"=meredith, "final"=finalTree)
class(tmp) <- "multiPhylo"
mltt.plot(tmp)
plot(branching.times(tree), table(tree$edge[,1]), ylab="# of Polytomies", xlab="Time (mybp)", pch=21, col="green", bg="green", xlim=c(0, 250))
points(branching.times(finalTree), table(finalTree$edge[,1]), pch=21, col="red", bg="red")

