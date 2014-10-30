## # Building a mammal tree
## Load tree and data into R
require(ape)
require(aRbor)
require(devtools)
require(bayou)
setwd("~/repos/bmr/R/MammalTree")
taxa <- "mammals"

td <- readRDS(paste("../../output/data/", taxa, ".rds", sep=""))
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
tax <- foreach(i=1:length(taxalist)) %dopar% try(rotl::tnrs_match_names(taxalist[i]))
tax <- do.call(rbind, tax)

ottids <- as.character(tax$ott_id) 

ottids

synth_tree <- rotl::tol_induced_subtree(ott_ids=ottids)
pdf("mammalsynthtree.pdf", height=15)
plot(synth_tree, cex=0.25)
plot(tree, cex=0.25)
dev.off()

## Find all source trees with our focal taxa in them, filter out trees without branch lengths. 
studyTrees <- findsourcetrees(ottids)      
possibleTrees <- getStudyTrees(studyTrees)
possibleTrees$summary

timetrees <- possibleTrees$result[which(possibleTrees$summary$timeUnits=="Myr")]
timetrees <- lapply(timetrees, function(x) x[[1]])

problemtree <- list(studyId=studyTrees[1,1]$studyId, treeId=studyTrees[1,2]$treeId)
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
tr <- replaceHashes(target, BE_cals)
phypd8 <- PATHd8.phylo(tr$scion, tr$table, base = ".tmp_PATHd8", rm = FALSE)

##Pathd8 again with better cals:
#phypd8$node.label <- target$node.label
tr <- replaceHashes(phypd8, calibrations)
stids <- cumsum(unlist(ncals))
stids <- lapply(1:length(stids), function(x) c((1+c(0,stids)[x]):stids[x]))
excluded <- c(144,146,147,148,149, 151,152,153,154,171)
finalTree <- PATHd8.phylo(tr$scion, tr$table[-excluded,], base = ".tmp_PATHd8", rm = FALSE)
finalTree <- reorder(finalTree, "postorder")
BETree <- reorder(tr$scion, "postorder")
plotTree(BETree, colors='red', ftype="off")

#plot(c(0,0), c(0, 0), type="n", xlim=c(0, 170), ylim=c(0, 6), xaxt="n", yaxt="n",xlab="", ylab="", bty="n")
rsBETree <- BETree
rsBETree$edge.length <- rsBETree$edge.length/max(branching.times(rsBETree))

rsfTree <- finalTree
rsfTree$edge.length <- rsfTree$edge.length/max(branching.times(rsfTree))

pdf("trees.pdf", height=25)
plot(BETree, cex=0.25)
plot(finalTree, cex=0.25)
plotTree(rsBETree, color="red", ftype="off")
plotTree(rsfTree, add=TRUE, ftype="off",lty=2, ftype="off")
plot(branching.times(BETree), branching.times(finalTree))
dev.off()

div <- cophenetic(finalTree)
div2 <- cophenetic(BETree)
taxapair <- sample(finalTree$tip.label, 2, replace=FALSE)
#taxapair <- c("Microtus_pennsylvanicus", "Cricetus_cricetus")
finalTree$tip.label
res <- div[taxapair[1], taxapair[2]]
res2 <- div2[taxapair[1], taxapair[2]]
c(res, res2)/2
taxapair



#ttt <- lapply(1:nrow(tr$table), function(x) PATHd8.phylo(tr$scion, tr$table[x,], base = ".tmp_PATHd8", rm = FALSE))
#conccal <- calibrations[sapply(ttt, is.phylo),]
#tr <- replaceHashes(phypd8, conccal)



stcal <- sourceTreeCalibrations(synth_tree, taxalist)
stcal$summary
rotl::get_study_meta(study_id = "pg_2688")
