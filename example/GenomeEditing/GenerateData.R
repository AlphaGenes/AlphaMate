
library(package = "AlphaSimR")
FOUNDERPOP = runMacs(nInd = 1000,
                     nChr = 10,
                     segSites = 2000,
                     inbred = FALSE,
                     species = "CATTLE")
SP = SimParam$new(founderPop = FOUNDERPOP)
SP$setGender(gender = "yes_sys")
SP$addSnpChip(nSnpPerChr = 1000)
SP$addTraitA(nQtlPerChr = 1000, mean = 0, var = 1)

nGen = 20
Pop = newPop(rawPop = FOUNDERPOP)
for (Gen in 1:nGen) {
  # Gen = 1
  Fathers = selectInd(pop = Pop[Pop@gender == "M"], nInd =  25,                         use = "gv")
  Mothers = selectInd(pop = Pop[Pop@gender == "F"], nInd = Pop[Pop@gender == "F"]@nInd, use = "gv")
  Pop = randCross2(females = Mothers, males = Fathers, nCrosses = 1000)
}

PopEdited1  = editGenomeTopQtl(pop = Pop, ind = 1:Pop@nInd, nQtl =  1)
PopEdited5  = editGenomeTopQtl(pop = Pop, ind = 1:Pop@nInd, nQtl =  5)
PopEdited20 = editGenomeTopQtl(pop = Pop, ind = 1:Pop@nInd, nQtl = 20)

plot(PopEdited1@gv ~ Pop@gv);  abline(a = 0, b = 1)
plot(PopEdited5@gv ~ Pop@gv);  abline(a = 0, b = 1)
plot(PopEdited20@gv ~ Pop@gv); abline(a = 0, b = 1)

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

prepareGenderForAlphaMate(pop = Pop, file = "Gender.txt")

prepareNrmForAlphaMate(pop = Pop, file = "Nrm.txt")

prepareCriterionForAlphaMate(pop = Pop,                          file = "Criterion.txt")
prepareCriterionForAlphaMate(pop = Pop, popEdited = PopEdited1,  file = "CriterionEdit1.txt")
prepareCriterionForAlphaMate(pop = Pop, popEdited = PopEdited5,  file = "CriterionEdit5.txt")
prepareCriterionForAlphaMate(pop = Pop, popEdited = PopEdited20, file = "CriterionEdit20.txt")

sink(file = "AlphaMateSpec.txt")
cat("GenderFile                  , Gender.txt\n")
cat("NrmMatrixFile               , Nrm.txt\n")
cat("SelCriterionFile            , Criterion.txt\n")
cat("OutputBasename              , Edit0_\n")
cat("NumberOfMatings             , 500\n")
cat("NumberOfMaleParents         ,  25\n")
cat("EqualizeMaleContributions   , Yes\n")
cat("NumberOfFemaleParents       , 500\n")
cat("EqualizeFemaleContributions , Yes\n")
cat("EvaluateFrontier            , Yes\n")
cat("Stop\n")
sink()

sink(file = "AlphaMateSpecEdit1.txt")
cat("GenderFile                  , Gender.txt\n")
cat("NrmMatrixFile               , Nrm.txt\n")
cat("SelCriterionFile            , CriterionEdit1.txt\n")
cat("OutputBasename              , Edit1_\n")
cat("NumberOfMatings             , 500\n")
cat("NumberOfMaleParents         ,  25\n")
cat("EqualizeMaleContributions   , Yes\n")
cat("NumberOfFemaleParents       , 500\n")
cat("EqualizeFemaleContributions , Yes\n")
cat("EvaluateFrontier            , Yes\n")
cat("PAGEMales                   , Yes\n")
cat("PAGEMalesMax                ,   5\n")
cat("Stop\n")
sink()

sink(file = "AlphaMateSpecEdit5.txt")
cat("GenderFile                  , Gender.txt\n")
cat("NrmMatrixFile               , Nrm.txt\n")
cat("SelCriterionFile            , CriterionEdit5.txt\n")
cat("OutputBasename              , Edit5_\n")
cat("NumberOfMatings             , 500\n")
cat("NumberOfMaleParents         ,  25\n")
cat("EqualizeMaleContributions   , Yes\n")
cat("NumberOfFemaleParents       , 500\n")
cat("EqualizeFemaleContributions , Yes\n")
cat("EvaluateFrontier            , Yes\n")
cat("PAGEMales                   , Yes\n")
cat("PAGEMalesMax                ,   5\n")
cat("Stop\n")
sink()

sink(file = "AlphaMateSpecEdit20.txt")
cat("GenderFile                  , Gender.txt\n")
cat("NrmMatrixFile               , Nrm.txt\n")
cat("SelCriterionFile            , CriterionEdit20.txt\n")
cat("OutputBasename              , Edit20_\n")
cat("NumberOfMatings             , 500\n")
cat("NumberOfMaleParents         ,  25\n")
cat("EqualizeMaleContributions   , Yes\n")
cat("NumberOfFemaleParents       , 500\n")
cat("EqualizeFemaleContributions , Yes\n")
cat("EvaluateFrontier            , Yes\n")
cat("PAGEMales                   , Yes\n")
cat("PAGEMalesMax                ,   5\n")
cat("Stop\n")
sink()