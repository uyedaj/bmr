## Load in command args:
args <- commandArgs(TRUE)
print(args)

## Load tree and data into R
require(ape)
require(aRbor)
require(bayou)
#require(devtools)
#load_all("~/repos/bayou/bayou_1.0")
#install_github("uyedaj/bayou", ref="dev")
td <- readRDS(paste("../output/data/", args[[1]], ".rds", sep=""))
head(td$dat)



