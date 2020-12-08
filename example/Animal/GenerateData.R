
# install.packages(pkg = "AlphaSimR")
# devtools::install_bitbucket("hickeyjohnteam/AlphaMME")
library(package = "AlphaSimR") # for simulation
# install.packages(pkg = "pedigreemm")
library(package = "pedigreemm") # for pedigree NRM

# Small simulation
nInd = 20
nFathers = 5

# Larger simulation
nInd = 1000
nFathers = 25

FOUNDERPOP = runMacs(nInd = nInd,
                     nChr = 10,
                     segSites = 2000,
                     inbred = FALSE,
                     species = "CATTLE")
SP = SimParam$new(founderPop = FOUNDERPOP)
SP$setSexes(sexes = "yes_sys")
SP$addSnpChip(nSnpPerChr = 1000)
SP$addTraitA(nQtlPerChr = 1000, mean = 0, var = 1)
SP$setTrackPed(isTrackPed = TRUE)

nGen = 20
Pop = newPop(rawPop = FOUNDERPOP)
for (Gen in 1:nGen) {
  # Gen = 1
  Fathers = selectInd(pop = Pop[Pop@sex == "M"], nInd =  nFathers,                   use = "gv")
  Mothers = selectInd(pop = Pop[Pop@sex == "F"], nInd = Pop[Pop@sex == "F"]@nInd, use = "gv")
  Pop = randCross2(females = Mothers, males = Fathers, nCrosses = nInd)
}

prepareSexForAlphaMate = function(pop, file) {
  tmp = data.frame(id = pop@id, sex = pop@sex)
  tmp$sex = as.numeric(factor(pop@sex, levels = c("M", "F")))
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

preparePedigreeNrmForAlphaMate = function(SimParam, pop = NULL, file = NULL, fileInverse = NULL, generations = NULL) {
  ped = cbind(individual = 1:nrow(SimParam$pedigree), father = SimParam$pedigree[, "father"], mother = SimParam$pedigree[, "mother"])
  if (!is.null(pop)) {
    cat("TODO: expand code to use only parts of the pedigree that is relevant to a pop and given number of generations\n")
  }
  ped2 = pedigree(sire = ped[, "father"], dam = ped[, "mother"], label = ped[, "individual"])
  if (!is.null(file)) {
    # U = relfactor(ped = ped2)
    # A = crossprod(U)
    # rownames(A) = ped2@label
    A = getA(ped = ped2)
    A = as.matrix(A)
    colnames(A) = NULL
    MASS::write.matrix(x = cbind(rownames(A), A), file = paste(file))
  }
  if (!is.null(fileInverse)) {
    # TInv = as(object = ped2, Class = "sparseMatrix")
    # DInv = Diagonal(x = 1 / Dmat(ped = ped2))
    # AInv = crossprod(sqrt(DInv) %*% TInv)
    # rownames(AInv) = ped2@label
    AInv = getAInv(ped = ped2)
    AInv = as.matrix(AInv)
    colnames(AInv) = NULL
    MASS::write.matrix(x = cbind(rownames(AInv), AInv), file = fileInverse)
  }
}

prepareGenomicNrmForAlphaMate = function(pop, file, snpChip = 1) {
  MASS::write.matrix(x = cbind(pop@id, AlphaMME::calcGIbs(X = pullSnpGeno(pop = pop, snpChip = snpChip))),
                     file = file)
}

prepareSexForAlphaMate(pop = Pop, file = "Sex.txt")

preparePedigreeNrmForAlphaMate(SimParam = SP, file = "PedigreeNrm.txt", fileInverse = "PedigreeNrmInverse.txt")

prepareGenomicNrmForAlphaMate(pop = Pop, file = "GenomicNrm.txt")

prepareCriterionForAlphaMate(pop = Pop, file = "Criterion.txt")

sink(file = "AlphaMateSpec.txt")
cat("GenderFile                  , Sex.txt\n")
cat("NrmMatrixFile               , Nrm.txt\n")
cat("SelCriterionFile            , Criterion.txt\n")
cat("NumberOfMatings             , 500\n")
cat("NumberOfMaleParents         ,  25\n")
cat("EqualizeMaleContributions   , Yes\n")
cat("NumberOfFemaleParents       , 500\n")
cat("EqualizeFemaleContributions , Yes\n")
cat("TargetDegree                ,  30\n")
cat("TargetMinInbreedingPct      ,  75\n")
cat("Stop\n")
sink()
