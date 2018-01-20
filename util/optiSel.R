
# ---- Code requirements ----

# install.packages(pkg=c("tidyverse", "optiSel", mipfp"))
library(package = "tidyverse") # for tidy data handling
library(package = "optiSel") # for Optimum Contribution Selection via quadratic programming
library(package = "mipfp") # for Iterative proportional fitting (used in MateAtRandom)
source(file = "Functions.R")

# ---- Data ----

# Breeding values
Data = read.table(file="../example/test/Ebv.txt")[, -3]
colnames(Data) = c("Indiv", "Ebv1")
# ... standardize them
Tmp = scale(Data$Ebv1)
EbvMean = attr(Tmp, "scaled:center")
EbvSd   = attr(Tmp, "scaled:scale")
Data$Ebv1S = Tmp[, 1]

# Add gender
Data$Sex = c("male", "female")[read.table(file="../example/test/Gender.txt")[, 2]]

# Relationship matrix
Nrm = as.matrix(read.table(file="../example/test/RelMat.txt")[, -1])
nInd = nrow(Data)
dimnames(Nrm) = list(1:nInd, 1:nInd)
# ... make it coancestry
Kin = Nrm / 2

str(Data); head(Data)
str(Nrm); Nrm[1:10, 1:10]

# Define maximum possible contributions
MaxUse = rep(4 / (2 * nInd), times=nInd)
names(MaxUse) = 1:nInd

# Construct candidates object
Candidates = candes(phen=Data, sKin=Kin, N=nInd)

# ---- Optimisation for maximal selection criterion ----

ConstraintsMaxSelCriterion = list()
ConstraintsMaxSelCriterion$ub = MaxUse
ContributionsMaxSelCriterion = opticont(method="max.Ebv1", cand=Candidates, con=ConstraintsMaxSelCriterion)
ContributionsMaxSelCriterion$info
ContributionsMaxSelCriterion$mean
ContributionsMaxSelCriterion$parent$nOff = noffspring(ContributionsMaxSelCriterion$parent, N=nInd)$nOff
ContributionsMaxSelCriterion$parent %>%
  filter(nOff > 0) %>%
  arrange(desc(nOff))

# ---- Optimisation for minimum coancestry ----

ConstraintsMinCoancestry = list()
ConstraintsMinCoancestry$ub = MaxUse
ContributionsMinCoancestry = opticont(method="min.sKin", cand=Candidates, con=ConstraintsMinCoancestry)
ContributionsMinCoancestry$info
ContributionsMinCoancestry$mean
ContributionsMinCoancestry$parent$nOff = noffspring(ContributionsMinCoancestry$parent, N=nInd)$nOff
ContributionsMinCoancestry$parent %>%
  filter(nOff > 0) %>%
  arrange(desc(nOff))

# ---- Penalty degrees ----

CurrentCoancestry = Candidates$current["sKin", "Val"]

(SelCriterionAtMaxSelCriterion     = ContributionsMaxSelCriterion$mean$Ebv1)
(SelCriterionStdAtMaxSelCriterion = ContributionsMaxSelCriterion$mean$Ebv1S)
(CoancestryAtMaxSelCriterion       = ContributionsMaxSelCriterion$mean$sKin)
(CoancestryRateAtMaxSelCriterion   = Coancestry2CoancestryRate(CurrentCoancestry=CurrentCoancestry,
                                                               FutureCoancestry=CoancestryAtMaxSelCriterion))

(SelCriterionAtMinCoancestry     = ContributionsMinCoancestry$mean$Ebv1)
(SelCriterionStdAtMinCoancestry = ContributionsMinCoancestry$mean$Ebv1S)
(CoancestryAtMinCoancestry       = ContributionsMinCoancestry$mean$sKin)
(CoancestryRateAtMinCoancestry   = Coancestry2CoancestryRate(CurrentCoancestry=CurrentCoancestry,
                                                             FutureCoancestry=CoancestryAtMinCoancestry))

TargetDegree = 30

(MaxCriterionAtTargetDegree     = Degree2MaxCriterionPct(Degree=TargetDegree))
(SelCriterionStdAtTargetDegree = Degree2SelCriterionStd(Degree=TargetDegree,
                                                        MinSelCriterionStd=SelCriterionStdAtMinCoancestry,
                                                        MaxSelCriterionStd=SelCriterionStdAtMaxSelCriterion))
(SelCriterionAtTargetDegree     = SelCriterionStd2SelCriterion(SelCriterionStd=SelCriterionStdAtTargetDegree,
                                                               Mean=EbvMean, Sd=EbvSd))
(MinCoancestryPctAtTargetDegree = Degree2MinCoancestryPct(Degree=TargetDegree))
(CoancestryRateAtTargetDegree   = Degree2CoancestryRate(Degree=TargetDegree,
                                                        MinCoancestryRate=CoancestryRateAtMinCoancestry,
                                                        MaxCoancestryRate=CoancestryRateAtMaxSelCriterion))
(CoancestryAtTargetDegree       = CoancestryRate2Coancestry(CoancestryRate=CoancestryRateAtTargetDegree,
                                                            CurrentCoancestry=CurrentCoancestry))

# ---- Optimisation for optimal objective ----

ConstraintsOpt = list()
ConstraintsOpt$ub = MaxUse
ConstraintsOpt$ub.sKin = CoancestryAtTargetDegree
ContributionsOpt = opticont(method="max.Ebv1", cand=Candidates, con=ConstraintsOpt)
ContributionsOpt$info
ContributionsOpt$mean
ContributionsOpt$parent$nOff = noffspring(ContributionsOpt$parent, N=nInd)$nOff
ContributionsOpt$parent %>%
  filter(nOff > 0) %>%
  arrange(desc(nOff))

# ---- Mating plan ---

# Random mating
MatingPlan = MateAtRandom(nMatings=nInd / 10,
                          Ind=ContributionsOpt$parent$Indiv,
                          nContributions = ContributionsOpt$parent$nOff,
                          Gender = (ContributionsOpt$parent$Sex == "female") + 1,
                          AllowSelfing = FALSE,
                          AllowRepeated = FALSE)

MatingRandomPoissonFamilySize = Mating2Progeny(x=MatingPlan, nProgenyPerMating=10)
head(MatingRandomPoissonFamilySize); nrow(MatingRandomPoissonFamilySize)
MatingRandomFixedFamilySize   = Mating2Progeny(x=MatingPlan, nProgenyPerMating=10, Poisson=FALSE)
head(MatingRandomFixedFamilySize); nrow(MatingRandomFixedFamilySize)

# Minimize expected inbreeding of progeny
MatingMinInb <- matings(ContributionsOpt$parent, Kin=Kin, N=nInd)
head(MatingMinInb)

