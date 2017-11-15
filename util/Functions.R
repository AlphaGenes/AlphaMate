
# ---- SelCriterion2SelIntensity --------------------------------------------------------------

#' Selection criterion to selection intensity
#'
#' `SelCriterion2SelIntensity` converts selection criterion to selection intensity
#'
#' @param SelCriterion numeric, selection criterion
#' @param Mean numeric, mean of selection criterion; if `NULL` computed from `SelCriterion`
#' @param Sd numeric, standard deviation of selection criterion; if `NULL` computed from `SelCriterion`
#'
#' @return numeric, selection intensity
#'
#' @export
#'
#' @examples
#' x <- rnorm(n=10, mean=10, sd=1)
#' SelCriterion2SelIntensity(x, Mean=10, Sd=1)
SelCriterion2SelIntensity = function(SelCriterion, Mean=NULL, Sd=NULL) {
  if (is.null(Mean)) {
    Mean = mean(SelCriterion)
  }
  if (is.null(Sd)) {
    Sd = sd(SelCriterion)
  }
  (SelCriterion - Mean) / Sd
}

# ---- SelIntensity2SelCriterion --------------------------------------------------------------

#' Selection intensity to selection criterion
#'
#' `SelIntensity2SelCriterion` converts selection intensity to selection criterion
#'
#' @param SelIntensity numeric, selection criterion
#' @param Mean numeric, mean of selection criterion
#' @param Sd numeric, standard deviation of selection criterion
#'
#' @return numeric, selection criterion
#'
#' @export
#'
#' @examples
#' x <- rnorm(n=10, mean=10, sd=1)
#' SelIntensity2SelCriterion(SelCriterion2SelIntensity(x, Mean=10, Sd=1))
SelIntensity2SelCriterion = function(SelIntensity, Mean, Sd) {
  SelIntensity * Sd + Mean
}

# ---- Coancestry2CoancestryRate --------------------------------------------------------------

#' Coancestries to coancestry rate
#'
#' `Coancestry2CoancestryRate` computes coancestry rate from current and future coancestries
#'
#' @param CurrentCoancestry numeric, current coancestry
#' @param FutureCoancestry numeric, future coancestry
#'
#' @return numeric, coancestry rate
#'
#' @export
#'
#' @examples
#' Coancestry2CoancestryRate(CurrentCoancestry=0.1, FutureCoancestry=0.2)
Coancestry2CoancestryRate = function(CurrentCoancestry, FutureCoancestry) {
  Diff    = FutureCoancestry - CurrentCoancestry
  MaxDiff =              1.0 - CurrentCoancestry
  if (MaxDiff == 0) {
    if (Diff >= 0) {
      CoancestryRate = 1.0
    } else if (Diff == 0) {
      CoancestryRate = 0.0
    } else {
      CoancestryRate = -1.0
    }
  } else {
    CoancestryRate = Diff / MaxDiff
  }
  CoancestryRate
}

# ---- CoancestryRate2Coancestry --------------------------------------------------------------

#' Future coancestry from current coancestry and coancestry rate
#'
#' `CoancestryRate2Coancestry` converts coancestry rate and current coancestry to future coancestry
#'
#' @param CoancestryRate numeric, coancestry rate
#' @param CurrentCoancestry numeric, current coancestry
#'
#' @return numeric, future coancestry
#'
#' @export
#'
#' @examples
#' CoancestryRate2Coancestry(CoancestryRate=Coancestry2CoancestryRate(CurrentCoancestry=0.1, FutureCoancestry=0.2),
#'                           CurrentCoancestry=0.1)
CoancestryRate2Coancestry = function(CoancestryRate, CurrentCoancestry) {
  CoancestryRate + (1.0 - CoancestryRate) * CurrentCoancestry
}

# ---- MinCoancestryPct2Degree --------------------------------------------------------------

#' Frontier degree from MinCoancestryPct
#'
#' `MinCoancestryPct2Degree` converts MinCoancestryPct to frontier degree
#'
#' @param MinCoancestryPct numeric, percentage of minimum coancestry achieved (100 means we achieved the minimum possible coancestry)
#'
#' @return numeric, frontier degree
#'
#' @export
#'
#' @details Assuming unit circular selection/coancestry frontier and given the
#' percentage of minimum inbreeding achieved (x-axis) we can evaluate the angle
#' between the maximum selection line and the solution line recognising that
#' sin(angle) = opposite/hypothenuse, where opposite = MinCoancestryPct/100 and
#' hypothenuse = 1 (unit circle). Then angle = asin(MinCoancestryPct/100).
#'
#' @examples
#' MinCoancestryPct2Degree(  0)
#' MinCoancestryPct2Degree( 25)
#' MinCoancestryPct2Degree( 50)
#' MinCoancestryPct2Degree( 75)
#' MinCoancestryPct2Degree(100)
MinCoancestryPct2Degree = function(MinCoancestryPct) {
  asin(MinCoancestryPct / 100.0) * 180.0 / pi
}

# ---- Degree2MinCoancestryPct --------------------------------------------------------------

#' MinCoancestryPct from Frontier degree
#'
#' `Degree2MinCoancestryPct` convverts frontier degree to MinCoancestryPct
#'
#' @param Degree numeric, frontier degree
#'
#' @return numeric, MinCoancestryPct (percentage of minimum coancestry achieved (100 means we achieved the minimum possible coancestry))
#'
#' @export
#'
#' @details Assuming unit circular selection/coancestry frontier and given the
#' angle between the maximum selection line and the solution line we can evaluate
#' the percentage of maximum selection intensity achieved (y-axis) by recognising
#' that sin(angle) = opposite/hypothenuse, where opposite = MinCoancestryPct/100
#' and hypothenuse = 1 (unit circle). Then MinCoancestryPct = sin(angle) * 100.
#'
#' @examples
#' Degree2MinCoancestryPct( 0)
#' Degree2MinCoancestryPct(45)
#' Degree2MinCoancestryPct(90)
Degree2MinCoancestryPct = function(Degree) {
  sin(Degree * pi / 180.0) * 100.0
}

# ---- MaxCriterionPct2Degree --------------------------------------------------------------

#' Frontier degree from MaxCriterionPct
#'
#' `MaxCriterionPct2Degree` converts MaxCriterionPct to frontier degree
#'
#' @param MaxCriterionPct numeric, percentage of maximum criterion achieved (100 means we achieved the maximum possible selection intensity)
#'
#' @return numeric, frontier degree
#'
#' @export
#'
#' @details Assuming unit circular selection/coancestry frontier and given the
#' percentage of maximum selection intensity achieved (y-axis) we can evaluate
#' the angle between the maximum selection line and the solution line recognising
#' that cos(angle) = adjacent/hypothenuse, where adjacent = MaxCriterionPct/100
#' and hypothenuse = 1 (unit circle). Then angle = acos(MaxCriterionPct/100).
#'
#' @examples
#' MaxCriterionPct2Degree(  0)
#' MaxCriterionPct2Degree( 25)
#' MaxCriterionPct2Degree( 50)
#' MaxCriterionPct2Degree( 75)
#' MaxCriterionPct2Degree(100)
MaxCriterionPct2Degree = function(MaxCriterionPct) {
  acos(MaxCriterionPct / 100.0) * 180.0 / pi
}

# ---- Degree2MaxCriterionPct --------------------------------------------------------------

#' MaxCriterionPct from Frontier degree
#'
#' `Degree2MaxCriterionPct` convverts frontier degree to MaxCriterionPct
#'
#' @param Degree numeric, frontier degree
#'
#' @return numeric, MaxCriterionPct (percentage of maximum criterion achieved (100 means we achieved the maximum possible selection intensity))
#'
#' @export
#'
#' @details Assuming unit circular selection/coancestry frontier and given the
#' angle between the maximum selection line and the solution line we can evaluate
#' the percentage of maximum selection intensity achieved (y-axis) by recognising
#' that cos(angle) = adjacent/hypothenuse, where adjacent = MaxCriterionPct/100
#' and hypothenuse = 1 (unit circle). Then MaxCriterionPct = cos(angle) * 100.
#'
#' @examples
#' Degree2MaxCriterionPct( 0)
#' Degree2MaxCriterionPct(45)
#' Degree2MaxCriterionPct(90)
Degree2MaxCriterionPct = function(Degree) {
  cos(Degree * pi / 180.0) * 100.0
}

# ---- MinCoancestryPct2CoancestryRate --------------------------------------------------------------

#' Coancestry rate from MinCoancestryPct
#'
#' `MinCoancestryPct2CoancestryRate` converts MinCoancestryPct to coancestry rate
#'
#' @param MinCoancestryPct numeric, MinCoancestryPct of a solution
#' @param MinCoancestryRate numeric, Minimum possible coancestry rate
#' @param MaxCoancestryRate numeric, Maximum possible coancestry rate
#'
#' @return numeric, coancestry rate
#'
#' @export
#'
#' @details Compute difference between Max and Min and take (1 - MinCoancestryPct)
#' of the difference as we express MinCoancestryPct as 100% when solution achieves
#' Min coancestry
#'
#' @examples
#' MinCoancestryPct2CoancestryRate(MinCoancestryPct=100, MinCoancestryRate=0.1, MaxCoancestryRate=0.2)
MinCoancestryPct2CoancestryRate = function(MinCoancestryPct, MinCoancestryRate, MaxCoancestryRate) {
  MinCoancestryRate + (100.0 - MinCoancestryPct) / 100.0 * (MaxCoancestryRate - MinCoancestryRate)
}

# ---- CoancestryRate2MinCoancestryPct --------------------------------------------------------------

#' MinCoancestryPct from coancestry rate
#'
#' `CoancestryRate2MinCoancestryPct` converts coancestry rate to MinCoancestryPct
#'
#' @param CoancestryRate numeric, coancestry rate at a given MinCoancestryPct
#' @param MinCoancestryRate numeric, Minimum possible coancestry rate
#' @param MaxCoancestryRate numeric, Maximum possible coancestry rate
#'
#' @return numeric, MinCoancestryPct
#'
#' @export
#'
#' @details Compute difference between Max and Min and take (1 - MinCoancestryPct)
#' of the difference as we express MinCoancestryPct as 100% when solution achieves
#' Min coancestry
#'
#' @examples
#' CoancestryRate2MinCoancestryPct(CoancestryRate=0.1, MinCoancestryRate=0.1, MaxCoancestryRate=0.2)
CoancestryRate2MinCoancestryPct = function(CoancestryRate, MinCoancestryRate, MaxCoancestryRate) {
  Diff    = MaxCoancestryRate - CoancestryRate
  MaxDiff = MaxCoancestryRate - MinCoancestryRate
  if (MaxDiff == 0) {
    if (Diff >= 0) {
      MinCoancestryPct = 100.0
    } else {
      MinCoancestryPct =   0.0
    }
  } else {
    MinCoancestryPct = Diff / MaxDiff * 100.0
  }
  MinCoancestryPct
}

# ---- MaxCriterionPct2SelIntensity --------------------------------------------------------------

#' Selection intensity from MaxCriterionPct
#'
#' `MaxCriterionPct2SelIntensity` convert MaxCriterionPct to selection intensity
#'
#' @param MaxCriterionPct numeric, MaxCriterionPct of a solution
#' @param MinSelIntensity numeric, Minimum possible selection intensity
#' @param MaxSelIntensity numeric, Maximum possible selection intensity
#'
#' @return numeric, selection intensity
#'
#' @export
#'
#' @details Compute difference between Max and Min and take MaxCriterionPct of
#' the difference as we express MaxCriterionPct as 100% when solution achieves
#' Max criterionCompute difference between Max and Min and take (1 - MinCoancestryPct)
#' of the difference as we express MinCoancestryPct as 100% when solution achieves
#' Min coancestry
#'
#' @examples
#' MaxCriterionPct2SelIntensity(MaxCriterionPct=100, MinSelIntensity=0.5, MaxSelIntensity=2.0)
MaxCriterionPct2SelIntensity = function(MaxCriterionPct, MinSelIntensity, MaxSelIntensity) {
  MinSelIntensity + MaxCriterionPct / 100.0 * (MaxSelIntensity - MinSelIntensity)
}

# ---- SelIntensity2MaxCriterionPct --------------------------------------------------------------

#' MaxCriterionPct from selection intensity
#'
#' `SelIntensity2MaxCriterionPct` convert selection intensity to MaxCriterionPct
#'
#' @param SelIntensity numeric, Selection intensity of a solution
#' @param MinSelIntensity numeric, Minimum possible selection intensity
#' @param MaxSelIntensity numeric, Maximum possible selection intensity
#'
#' @return numeric, MaxCriterionPct
#'
#' @export
#'
#' @examples
#' SelIntensity2MaxCriterionPct(SelIntensity=2, MinSelIntensity=0.5, MaxSelIntensity=2)
SelIntensity2MaxCriterionPct = function(SelIntensity, MinSelIntensity, MaxSelIntensity) {
  Diff = SelIntensity - MinSelIntensity
  MaxDiff = MaxSelIntensity - MinSelIntensity
  if (MaxDiff == 0) {
    if (Diff >= 0) {
      MaxCriterionPct = 100.0
    } else {
      MaxCriterionPct =   0.0
    }
  } else {
    MaxCriterionPct = Diff / MaxDiff * 100.0
  }
  MaxCriterionPct
}

# ---- Degree2SelIntensity --------------------------------------------------------------

#' Selection intensity from frontier degree
#'
#' `Degree2SelIntensity` convert degree to selection intensity
#'
#' @param Degree numeric, Degree
#' @param MinSelIntensity numeric, Minimum possible selection intensity
#' @param MaxSelIntensity numeric, Maximum possible selection intensity
#'
#' @return numeric, selection intensity
#'
#' @export
#'
#' @examples
#' Degree2SelIntensity(Degree=45, MinSelIntensity=0.5, MaxSelIntensity=2)
Degree2SelIntensity = function(Degree, MinSelIntensity, MaxSelIntensity) {
  MaxCriterionPct2SelIntensity(MaxCriterionPct=Degree2MaxCriterionPct(Degree=Degree),
                               MinSelIntensity=MinSelIntensity,
                               MaxSelIntensity=MaxSelIntensity)
}

# ---- SelIntensity2Degree --------------------------------------------------------------

#' Frontier degree from selection intensity
#'
#' `SelIntensity2Degree` convert selection intensity to degree
#'
#' @param SelIntensity numeric, selection intensity
#' @param MinSelIntensity numeric, Minimum possible selection intensity
#' @param MaxSelIntensity numeric, Maximum possible selection intensity
#'
#' @return numeric, frontier degree
#'
#' @export
#'
#' @examples
#' SelIntensity2Degree(SelIntensity=1.56, MinSelIntensity=0.5, MaxSelIntensity=2)
SelIntensity2Degree = function(SelIntensity, MinSelIntensity, MaxSelIntensity) {
  MaxCriterionPct2Degree(MaxCriterionPct=SelIntensity2MaxCriterionPct(SelIntensity=SelIntensity,
                                                                      MinSelIntensity=MinSelIntensity,
                                                                      MaxSelIntensity=MaxSelIntensity))
}

# ---- Degree2CoancestryRate --------------------------------------------------------------

#' Coancestry rate from frontier degree
#'
#' `Degree2CoancestryRate` convert degree to coancestry rate
#'
#' @param Degree numeric, Degree
#' @param MinCoancestryRate numeric, Minimum possible coancestry rate
#' @param MaxCoancestryRate numeric, Maximum possible coancestry rate
#'
#' @return numeric, coancestry rate
#'
#' @export
#'
#' @examples
#' Degree2CoancestryRate(Degree=45, MinCoancestryRate=0.1, MaxCoancestryRate=1)
Degree2CoancestryRate = function(Degree, MinCoancestryRate, MaxCoancestryRate) {
  MinCoancestryPct2CoancestryRate(MinCoancestryPct=Degree2MinCoancestryPct(Degree=Degree),
                                  MinCoancestryRate=MinCoancestryRate,
                                  MaxCoancestryRate=MaxCoancestryRate)
}

# ---- CoancestryRate2Degree --------------------------------------------------------------

#' Frontier degree from coancestry rate
#'
#' `CoancestryRate2Degree` convert coancestry rate to degree
#'
#' @param CoancestryRate numeric, coancestry rate
#' @param MinCoancestryRate numeric, Minimum possible coancestry rate
#' @param MaxCoancestryRate numeric, Maximum possible coancestry rate
#'
#' @return numeric, frontier degree
#'
#' @export
#'
#' @examples
#' CoancestryRate2Degree(CoancestryRate=0.36, MinCoancestryRate=0.1, MaxCoancestryRate=1)
CoancestryRate2Degree = function(CoancestryRate, MinCoancestryRate, MaxCoancestryRate) {
  MinCoancestryPct2Degree(MinCoancestryPct=CoancestryRate2MinCoancestryPct(CoancestryRate=CoancestryRate,
                                                                           MinCoancestryRate=MinCoancestryRate,
                                                                           MaxCoancestryRate=MaxCoancestryRate))
}

# ---- MateAtRandom ----

#' Mate/cross at random
#'
#' `MateAtRandom` mates/crosses individuals at random
#'
#' @param nMatings numeric, number of matings/crosses
#' @param Ind numeric/integer/character, a vector with individual names/ids
#' @param nContributions numeric, a vector with individual integer! contributions
#' @param Gender numeric, a vector with individual gender (1=males, 2=females),
#'        can be `NULL` when individuals are bi-sexual
#' @param AllowSelfing logical, should selfing be allowed (applies only when `Gender=NULL`)
#' @param AllowRepeated logical, allow repeated matings/crosses between two individuals
#'
#' @return data frame, a mating/crossing plan
#'
#' @export
#'
#' @examples
#' MateAtRandom(nMatings=3, Ind=1:4, nContributions=c(2, 1, 2, 2))
MateAtRandom = function(nMatings, Ind, nContributions, Gender=NULL, AllowSelfing=FALSE, AllowRepeated=FALSE) {
  # First remove non-contributers
  Sel = nContributions > 0
  nContributions = nContributions[Sel]
  Ind = Ind[Sel]
  if (!is.null(Gender)) {
    Gender = Gender[Sel]
  }
  # Specify the two groups of parents
  if (!is.null(Gender)) {
    Sel = Gender == 1
    Parent1 = Ind[Sel]
    nContributionsParent1 = nContributions[Sel]
    Parent2 = Ind[!Sel]
    nContributionsParent2 = nContributions[!Sel]
  } else {
    Parent1 = Ind
    nContributionsParent1 = nContributions
    Parent2 = Ind
    nContributionsParent2 = nContributions
  }
  # Setup the initial (seed) matrix
  nParent1 = length(Parent1)
  nParent2 = length(Parent2)
  Init = array(data=1, dim=c(nParent1, nParent2), dimnames=list(Parent1, Parent2))
  if (is.null(Gender) & !AllowSelfing) {
    diag(Init) = 0
  }
  # Run iterative proportional fitting to find matings
  # (could equally use a sampling approach based on nContributions as well; we sample
  #  the permissible matings later on anyhow)
  Result = Ipfp(seed=Init, target.list=list(1, 2), target.data=list(nContributionsParent1, nContributionsParent2))
  if (!Result$conv) {
    stop("Ipfp did not converge!")
  }
  # Sample matings
  Mating = Array2Vector(outer(Parent1, Parent2, paste, sep="-"))
  Prob = Array2Vector(Result$x.hat)
  Tmp = sample(x=Mating, size=nMatings, prob=Prob, replace=AllowRepeated)
  # Produce an easy to use data frame
  Tmp = strsplit(x=Tmp, split="-")
  Return = data.frame(Mating=1:nMatings, Parent1=rep(NA, times=nMatings), Parent2=rep(NA, times=nMatings))
  for (Mating in 1:nMatings) {
    Return[["Parent1"]][Mating] = Tmp[[Mating]][1]
    Return[["Parent2"]][Mating] = Tmp[[Mating]][2]
  }
  Return
}

# ---- Mating2Progeny ----

#' Generate progeny pedigree from a mating plan
#'
#' `Mating2Progeny` generates progeny pedigree from a mating plan
#'
#' @param x data frame, mating/crossing plan with columns Mating (1:n), Parent1, and Parent2
#' @param nProgenyPerMating numeric/integer, average number of progeny per mating
#' @param Poisson logical, should number of progeny per mating be sampled from a
#'        Poisson distribution with lambda=nProgenyPerMating or should it be fixed
#'        to nProgenyPerMating for all matings; Other distributions could be considered
#'        in the future
#'
#' @return data frame, a pedigree
#'
#' @export
#'
#' @examples
#' MatingPlan = MateAtRandom(nMatings=3, Ind=1:4, nContributions=c(2, 1, 2, 2))
#' Mating2Progeny(x=MatingPlan, nProgenyPerMating=10)
#' Mating2Progeny(x=MatingPlan, nProgenyPerMating=10, Poisson=FALSE)
Mating2Progeny = function(x, nProgenyPerMating, Poisson=TRUE) {
  nMatings = nrow(x)
  if (Poisson) {
    nPerMating = rpois(n=nMatings, lambda=nProgenyPerMating)
  } else {
    nPerMating = rep(nProgenyPerMating, times=nMatings)
  }
  n = sum(nPerMating)
  Return = data.frame(Mating=NA, Progeny=1:n, Parent1=NA, Parent2=NA)
  Row = 0
  for (Mating in 1:nMatings) {
    for (Progeny in 1:nPerMating[Mating]) {
      Row = Row + 1
      Return[["Mating"]][[Row]] = Mating
      Return[["Parent1"]][[Row]] = x[["Parent1"]][Mating]
      Return[["Parent2"]][[Row]] = x[["Parent2"]][Mating]
    }
  }
  Return
}
