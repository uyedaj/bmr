white <- read.csv("~/Dropbox/BodySize_Metabolism/Metabolic Rate Data_curated/Whiteetal_VertMetabolic_Appendix.csv")
head(white)
calculator <- function(DF)
{
  Q10 <- c()
  for(i in 1:dim(DF)[1])
  {
    Q10[i] <- (DF$q10smr[i]/DF$MR[i])^(10/(38-DF$temp[i]))
  }
  Q10
}


plot(DF$temp, log(q10))
q10 <- calculator(white)
white$q10 <- q10
white
head(white)
plot(white$temp, resid(lm(log(white$MR)~log(white$mean.mass))))
plot(white$temp, resid(lm(log(white$q10smr)~log(white$mean.mass))), pch=21, bg=white$color)
plot(white$temp, resid(lm(log(white$MR)~log(white$mean.mass))), pch=21, bg=white$color)

estMR38 <- function(sp){
  tmp <- filter(white, species==sp)
  y <- log(tmp$MR)
  x1 <- log(tmp$mass)
  x2 <- tmp$temp
  lmfit <- lm(y~x1+x2)
  newdat <- data.frame(x1=mean(tmp$mass, na.rm=TRUE), x2=38)
  BMR <- predict(lmfit, newdata=newdat)
  c(BMR, tmp$q10smr[!is.na(tmp$q10smr)])
}
mWhite <- names(table(white$species))[which(table(white$species)<2)]
hist(subset(white, species %in% mWhite)$q10, xlim=c(0,4), breaks=50000)
abline(v=c(1.0, 1.65, 2.21, 2.44, 2.09, 2.81), lty=2)

est38s <- sapply(mWhite,function(x) try(estMR38(x)))
est38s <- est38s[!(sapply(est38s, class)=="try-error")]
est38s <- est38s[sapply(est38s, length) ==2]
tmp <- do.call(rbind, est38s)
plot(tmp)

sp <- mWhite[1]

plot(white$temp, log(white$mean.mass))
