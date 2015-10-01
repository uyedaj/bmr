require(dplyr)
santost4 <- read.csv("../datasets/SantosSuppTable4.csv")
santost3 <- read.csv("../datasets/SantosSuppTable3.csv")

plot(log(santost4$mass.g), log(santost4$rmr.mean))

