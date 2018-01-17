
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

# Generate pools
Pop@gender = sample(x = c("M", "F"), size = 370, replace = TRUE)

prepareGenderForAlphaMate = function(pop, file) {
  tmp = data.frame(id = pop@id, genderRole = pop@gender)
  tmp$genderRole = as.numeric(factor(pop@gender, levels = c("M", "F")))
  write.table(x = tmp, file = file,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
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

prepareGenderForAlphaMate(pop = Pop, file = "Pools.txt")

prepareNrmForAlphaMate(pop = Pop, file = "Nrm.txt")

prepareCriterionForAlphaMate(pop = Pop, file = "Criterion.txt")

sink(file = "AlphaMateSpec.txt")
cat("GenderFile                  , Pools.txt\n")
cat("NrmMatrixFile               , Nrm.txt\n")
cat("SelCriterionFile            , Criterion.txt\n")
cat("NumberOfMatings             , 32\n")
cat("LimitMaleContributions      , Yes\n")
cat("LimitMaleContributionsMax   ,  4\n")
cat("LimitFemaleContributions    , Yes\n")
cat("LimitFemaleContributionsMax ,  4\n")
cat("TargetDegree                , 45\n")
cat("Stop\n")
sink()