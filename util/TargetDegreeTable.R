
source(file = "Functions.R")
Degree = seq(from = 0, to = 90, by = 5)
x = cbind(Degree, MaxCriterionPct = Degree2MaxCriterionPct(Degree), MinCoancestryPct = Degree2MinCoancestryPct(Degree))
round(x, digits = 1)
