
# ---- Code requirements ----

# install.packages(pkg = c("tidyverse", "optiSel", mipfp"))
library(package = "tidyverse") # for tidy data handling
library(package = "optiSel") # for Optimum Contribution Selection via quadratic programming
library(package = "mipfp") # for Iterative proportional fitting (used in MateAtRandom)
source(file = "Functions.R")

# ---- Data ----

# Breeding values
Data = read.table(file = "../example/Animal/Criterion.txt")[, -3]
colnames(Data) = c("Indiv", "Ebv1")
# ... standardize them
Tmp = scale(Data$Ebv1)
CurrentSelCriterionMean = attr(Tmp, "scaled:center")
CurrentSelCriterionSd   = attr(Tmp, "scaled:scale")
Data$Ebv1S = Tmp[, 1]

# Add gender
Data$Sex = c("male", "female")[read.table(file = "../example/Animal/Gender.txt")[, 2]]

# Relationship matrix
Nrm = as.matrix(read.table(file = "../example/Animal/Nrm.txt")[, -1])
nInd = nrow(Data)
dimnames(Nrm) = list(Data$Indiv, Data$Indiv)
# ... make it coancestry
Kin = Nrm / 2

str(Data); head(Data)
str(Nrm); Nrm[1:10, 1:10]
str(Kin); Kin[1:10, 1:10]

# Define maximum possible contributions
MaxUse = rep(4 / (2 * nInd), times = nInd)
names(MaxUse) = Data$Indiv

# Construct candidates object
Candidates = candes(phen = Data, sKin = Kin, N = nInd)

# ---- Optimisation for maximal selection criterion ----

ConstraintsMaxSelCriterion = list()
ConstraintsMaxSelCriterion$ub = MaxUse
ContributionsMaxSelCriterion = opticont(method = "max.Ebv1", cand = Candidates, con = ConstraintsMaxSelCriterion)
ContributionsMaxSelCriterion$info
ContributionsMaxSelCriterion$mean
ContributionsMaxSelCriterion$parent$nOff = noffspring(ContributionsMaxSelCriterion$parent, N = nInd)$nOff
ContributionsMaxSelCriterion$parent %>%
  filter(nOff > 0) %>%
  arrange(desc(nOff))

# ---- Optimisation for minimum coancestry ----

ConstraintsMinCoancestry = list()
ConstraintsMinCoancestry$ub = MaxUse
ContributionsMinCoancestry = opticont(method = "min.sKin", cand = Candidates, con = ConstraintsMinCoancestry)
ContributionsMinCoancestry$info
ContributionsMinCoancestry$mean
ContributionsMinCoancestry$parent$nOff = noffspring(ContributionsMinCoancestry$parent, N = nInd)$nOff
ContributionsMinCoancestry$parent %>%
  filter(nOff > 0) %>%
  arrange(desc(nOff))

# ---- Penalty degrees ----

CurrentCoancestry = Candidates$current["sKin", "Val"]

(SelCriterionAtMaxSelCriterion    = ContributionsMaxSelCriterion$mean$Ebv1)
(SelCriterionStdAtMaxSelCriterion = ContributionsMaxSelCriterion$mean$Ebv1S)
(CoancestryAtMaxSelCriterion      = ContributionsMaxSelCriterion$mean$sKin)
(CoancestryRateAtMaxSelCriterion  = Coancestry2CoancestryRate(CurrentCoancestry = CurrentCoancestry,
                                                               FutureCoancestry = CoancestryAtMaxSelCriterion))

(SelCriterionAtMinCoancestry    = ContributionsMinCoancestry$mean$Ebv1)
(SelCriterionStdAtMinCoancestry = ContributionsMinCoancestry$mean$Ebv1S)
(CoancestryAtMinCoancestry      = ContributionsMinCoancestry$mean$sKin)
(CoancestryRateAtMinCoancestry  = Coancestry2CoancestryRate(CurrentCoancestry = CurrentCoancestry,
                                                             FutureCoancestry = CoancestryAtMinCoancestry))

if (FALSE) {
  CurrentCoancestry = 0.17307
  CurrentSelCriterionMean = 7.86059
  CurrentSelCriterionSd = 1.56991

  # AlphaMate
  SelCriterionAtMaxSelCriterion = 10.2517
  CoancestryAtMaxSelCriterion   = 0.27191
  SelCriterionAtMinCoancestry   = 9.03691
  CoancestryAtMinCoancestry     = 0.17547

  # optiSel
  SelCriterionAtMaxSelCriterion = 10.2517
  CoancestryAtMaxSelCriterion   = 0.2746
  SelCriterionAtMinCoancestry   = 9.0713
  CoancestryAtMinCoancestry     = 0.1745

  (SelCriterionStdAtMaxSelCriterion =  SelCriterion2SelCriterionStd(SelCriterion = SelCriterionAtMaxSelCriterion,
                                                                    Mean = CurrentSelCriterionMean,
                                                                    Sd = CurrentSelCriterionSd))
  (CoancestryRateAtMaxSelCriterion  = Coancestry2CoancestryRate(CurrentCoancestry = CurrentCoancestry,
                                                                FutureCoancestry = CoancestryAtMaxSelCriterion))
  (SelCriterionStdAtMinCoancestry = SelCriterion2SelCriterionStd(SelCriterion = SelCriterionAtMinCoancestry,
                                                                 Mean = CurrentSelCriterionMean,
                                                                 Sd = CurrentSelCriterionSd))
  (CoancestryRateAtMinCoancestry  = Coancestry2CoancestryRate(CurrentCoancestry = CurrentCoancestry,
                                                              FutureCoancestry = CoancestryAtMinCoancestry))
}

TargetDegree = 30

(MaxCriterionAtTargetDegree     = Degree2MaxCriterionPct(Degree = TargetDegree))
(SelCriterionStdAtTargetDegree  = Degree2SelCriterionStd(Degree = TargetDegree,
                                                         MinSelCriterionStd = SelCriterionStdAtMinCoancestry,
                                                         MaxSelCriterionStd = SelCriterionStdAtMaxSelCriterion))
(SelCriterionAtTargetDegree     = SelCriterionStd2SelCriterion(SelCriterionStd = SelCriterionStdAtTargetDegree,
                                                               Mean = CurrentSelCriterionMean,
                                                               Sd = CurrentSelCriterionSd))
(MinCoancestryPctAtTargetDegree = Degree2MinCoancestryPct(Degree = TargetDegree))
(CoancestryRateAtTargetDegree   = Degree2CoancestryRate(Degree = TargetDegree,
                                                        MinCoancestryRate = CoancestryRateAtMinCoancestry,
                                                        MaxCoancestryRate = CoancestryRateAtMaxSelCriterion))
(CoancestryAtTargetDegree       = CoancestryRate2Coancestry(CoancestryRate = CoancestryRateAtTargetDegree,
                                                            CurrentCoancestry = CurrentCoancestry))

if (FALSE) {
  # AlphaMate
  SelCriterionAtOptTarget = 10.17
  CoancestryAtOptTarget = 0.220

  # optiSel
  SelCriterionAtOptTarget = 10.19
  CoancestryAtOptTarget = 0.225

  (SelCriterionStdAtOptTarget = SelCriterion2SelCriterionStd(SelCriterion = SelCriterionAtOptTarget,
                                                             Mean = CurrentSelCriterionMean,
                                                             Sd = CurrentSelCriterionSd))
  (MaxCriterionAtOptTarget = SelCriterionStd2MaxCriterionPct(SelCriterionStd = SelCriterionStdAtOptTarget,
                                                             MinSelCriterionStd = SelCriterionStdAtMinCoancestry,
                                                             MaxSelCriterionStd = SelCriterionStdAtMaxSelCriterion))
  (CoancestryRateAtOptTarget = Coancestry2CoancestryRate(CurrentCoancestry = CurrentCoancestry,
                                                         FutureCoancestry = CoancestryAtOptTarget))
  (MinCoancestryAtOptTarget = CoancestryRate2MinCoancestryPct(CoancestryRate = CoancestryRateAtOptTarget,
                                                              MinCoancestryRate = CoancestryRateAtMinCoancestry,
                                                              MaxCoancestryRate = CoancestryRateAtMaxSelCriterion))
  (DegreeAtOptTarget = MaxCriterionPct2Degree(MaxCriterionPct = MaxCriterionAtOptTarget))
}

# ---- Optimisation for optimal objective ----

ConstraintsOpt = list()
ConstraintsOpt$ub = MaxUse
ConstraintsOpt$ub.sKin = CoancestryAtTargetDegree
ContributionsOpt = opticont(method = "max.Ebv1", cand = Candidates, con = ConstraintsOpt)
ContributionsOpt$info
ContributionsOpt$mean
ContributionsOpt$parent$nOff = noffspring(ContributionsOpt$parent, N = nInd)$nOff
ContributionsOpt$parent %>%
  filter(nOff > 0) %>%
  arrange(desc(nOff))

# ---- Mating plan ---

# Random mating
MatingPlan = MateAtRandom(nMatings = nInd / 10,
                          Ind = ContributionsOpt$parent$Indiv,
                          nContributions = ContributionsOpt$parent$nOff,
                          Gender = (ContributionsOpt$parent$Sex == "female") + 1,
                          AllowSelfing = FALSE,
                          AllowRepeated = FALSE)

MatingRandomPoissonFamilySize = Mating2Progeny(x = MatingPlan, nProgenyPerMating = 10)
head(MatingRandomPoissonFamilySize); nrow(MatingRandomPoissonFamilySize)
MatingRandomFixedFamilySize   = Mating2Progeny(x = MatingPlan, nProgenyPerMating = 10, Poisson = FALSE)
head(MatingRandomFixedFamilySize); nrow(MatingRandomFixedFamilySize)

# Minimize expected inbreeding of progeny
MatingMinInb <- matings(ContributionsOpt$parent, Kin = Kin, N = nInd)
head(MatingMinInb)

