

PedMS <- read.table(file="~/Downloads/PedigreeAll.txt")
PedMS$Generation <- rep(1:(nrow(PedMS)/1000), each=1000)
PedMS$BV <- read.table(file="~/Downloads/EbvAll.txt")[, 2]
PedMS$BV <- PedMS$BV / sd(PedMS$BV[PedMS$Generation == 1])

PedTS <- read.table(file="~/Downloads/Pedigree_25Sires-1.txt")
PedTS$Generation <- rep(1:(nrow(PedTS)/1000), each=1000)
PedTS$BV <- read.table(file="~/Downloads/EbvAll_25Sires-1.txt")[, 2]
PedTS$BV <- PedTS$BV / sd(PedTS$BV[PedTS$Generation == 1])

library(package="doBy")
library(package="pedigreemm")

Ped2MS <- with(PedMS, pedigree(sire=PedMS[, 2], dam=PedMS[, 3], label=PedMS[, 1]))
PedMS$Inbreeding <- inbreeding(Ped2MS)

Ped2TS <- with(PedTS, pedigree(sire=PedTS[, 2], dam=PedTS[, 3], label=PedTS[, 1]))
PedTS$Inbreeding <- inbreeding(Ped2TS)

TmpMS <- summaryBy(BV ~ Generation, data=PedMS, FUN=descStat)
TmpTS <- summaryBy(BV ~ Generation, data=PedTS, FUN=descStat)

TmpFMS <- summaryBy(Inbreeding ~ Generation, data=PedMS, FUN=descStat)
TmpFTS <- summaryBy(Inbreeding ~ Generation, data=PedTS, FUN=descStat)

TmpFMS$RateOfInbreeding <- NA
TmpFMS$RateOfInbreeding[2:nrow(TmpFMS)] <- (TmpFMS$Inbreeding.mean[2:nrow(TmpFMS)] - TmpFMS$Inbreeding.mean[1:(nrow(TmpFMS)-1)]) / (1 - TmpFMS$Inbreeding.mean[1:(nrow(TmpFMS)-1)])
TmpFTS$RateOfInbreeding <- NA
TmpFTS$RateOfInbreeding[2:nrow(TmpFTS)] <- (TmpFTS$Inbreeding.mean[2:nrow(TmpFTS)] - TmpFTS$Inbreeding.mean[1:(nrow(TmpFTS)-1)]) / (1 - TmpFTS$Inbreeding.mean[1:(nrow(TmpFTS)-1)])

ColRoslinBlue   <- rgb(red= 55, green=152, blue=217, maxColorValue=255) ## blue
ColRoslinViolet <- rgb(red=193, green= 23, blue=115, maxColorValue=255) ## violet
ColRoslinGreen  <- rgb(red= 95, green=168, blue= 59, maxColorValue=255) ## green
ColRoslinOrange <- rgb(red=219, green= 79, blue= 15, maxColorValue=255) ## orange
ColRoslinGray   <- rgb(red= 83, green= 83, blue= 83, maxColorValue=255) ## orange

matplot(y=cbind(TmpTS$BV.mean, TmpMS$BV.mean), x=TmpMS$Generation,
        type="l", lwd=2, col=c(ColRoslinBlue, ColRoslinViolet), xlab="Generation", ylab="Genetic gain")
legend(x="topleft", legend=c("Truncation selection", "Mate selection"),
       lwd=2, bty="n", col=c(ColRoslinBlue, ColRoslinViolet))

matplot(y=cbind(TmpTS$BV.sd^2, TmpMS$BV.sd^2), x=TmpMS$Generation,
        type="l", lwd=2, col=c(ColRoslinBlue, ColRoslinViolet), xlab="Generation", ylab="Genetic variance")
legend(x="topright", legend=c("Truncation selection", "Mate selection"),
       lwd=2, bty="n", col=c(ColRoslinBlue, ColRoslinViolet))

matplot(y=cbind(TmpFTS$Inbreeding.mean, TmpFMS$Inbreeding.mean), x=TmpMS$Generation,
        type="l", lwd=2, col=c(ColRoslinBlue, ColRoslinViolet), xlab="Generation", ylab="Coef. of inbreeding")
legend(x="topleft", legend=c("Truncation selection", "Mate selection"),
       lwd=2, bty="n", col=c(ColRoslinBlue, ColRoslinViolet))

matplot(y=cbind(TmpFTS$RateOfInbreeding, TmpFMS$RateOfInbreeding), x=TmpMS$Generation,
        type="l", lwd=2, col=c(ColRoslinBlue, ColRoslinViolet), xlab="Generation", ylab="Rate of inbreeding")
legend(x="topright", legend=c("Truncation selection", "Mate selection"),
       lwd=2, bty="n", col=c(ColRoslinBlue, ColRoslinViolet))
