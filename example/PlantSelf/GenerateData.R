
library(package = "AlphaSimR")
FOUNDERPOP = runMacs(nInd = 370,
                     nChr = 10,
                     segSites = 2000,
                     inbred = TRUE,
                     species = "WHEAT")
SP = SimParam$new(founderPop = FOUNDERPOP)
SP$addSnpChip(nSnpPerChr = 1000)
SP$addTraitA(nQtlPerChr = 1000, mean = 0, var = 1)

nGen = 20
Pop = newPop(rawPop = FOUNDERPOP)
for (Gen in 1:nGen) {
  # Gen = 1
  Parents = selectInd(pop = Pop, nInd = 64, use = "gv")
  Pop = randCross(pop = Parents, nCrosses = 370)
}

prepareCriterionForAlphaMate = function(pop, popEdited = NULL, file) {
  tmp = data.frame(id = pop@id, crit = pop@gv)
  if (!is.null(popEdited)) {
    tmp$critEdited = popEdited@gv
  }
  write.table(x = tmp, file = file,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}

prepareNrmForAlphaMate = function(pop, file, snpChip = 1) {
  MASS::write.matrix(x = cbind(pop@id, calcGIbs(X = pullSnpGeno(pop = pop, snpChip = snpChip))),
                     file = file)
}

prepareNrmForAlphaMate(pop = Pop, file = "Nrm.txt")

prepareCriterionForAlphaMate(pop = Pop, file = "Criterion.txt")

sink(file = "AlphaMateSpec.txt")
cat("NrmMatrixFile         , Nrm.txt\n")
cat("SelCriterionFile      , Criterion.txt\n")
cat("NumberOfMatings       , 32\n")
cat("LimitContributions    , Yes\n")
cat("LimitContributionsMax , 10\n")
cat("AllowSelfing          , Yes\n")
cat("TargetDegree          , 30\n")
cat("Stop\n")
sink()