#!/usr/bin/env Rscript

# ---- Required packages and functions ----

RequiredPackages = c("optparse", "tidyverse")
for (Package in RequiredPackages) {
  Test = require(package = Package, quietly = TRUE, character.only = TRUE)
  if (!Test) {
    install.packages(pkgs = Package)
  }
}

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

# ---- Arguments ----

ArgumentsList = list(
  make_option(opt_str = c("-f", "--frontier"), type = "character", default = NULL,
              help = "frontier file name [default is %default]"),
  make_option(opt_str = "--frontier2", type = "character", default = NULL,
              help = "frontier file name 2 [default is %default]"),
  make_option(opt_str = "--frontier3", type = "character", default = NULL,
              help = "frontier file name 3 [default is %default]"),
  make_option(opt_str = "--frontier4", type = "character", default = NULL,
              help = "frontier file name 4 [default is %default]"),
  make_option(opt_str = "--frontier5", type = "character", default = NULL,
              help = "frontier file name 5 [default is %default]"),

  make_option(opt_str = c("-t", "--targets"), type = "character", default = NULL,
              help = "targets file name [default is %default]"),

  make_option(opt_str = c("-l", "--log"), type = "character", default = NULL,
              help = "optimisation log file name [default is %default]"),
  make_option(opt_str = "--log2", type = "character", default = NULL,
              help = "optimisation log file name 2 [default is %default]"),
  make_option(opt_str = "--log3", type = "character", default = NULL,
              help = "optimisation log file name 3 [default is %default]"),
  make_option(opt_str = "--log4", type = "character", default = NULL,
              help = "optimisation log file name 4 [default is %default]"),
  make_option(opt_str = "--log5", type = "character", default = NULL,
              help = "optimisation log file name 5 [default is %default]"),
  make_option(opt_str = "--log6", type = "character", default = NULL,
              help = "optimisation log file name 6 [default is %default]"),
  make_option(opt_str = "--log7", type = "character", default = NULL,
              help = "optimisation log file name 7 [default is %default]"),
  make_option(opt_str = "--log8", type = "character", default = NULL,
              help = "optimisation log file name 8 [default is %default]"),
  make_option(opt_str = "--log9", type = "character", default = NULL,
              help = "optimisation log file name 9 [default is %default]"),
  make_option(opt_str = "--log10", type = "character", default = NULL,
              help = "optimisation log file name 10 [default is %default]"),
  make_option(opt_str = "--log11", type = "character", default = NULL,
              help = "optimisation log file name 11 [default is %default]"),

  make_option(opt_str = c("-p", "--pop"), type = "character", default = NULL,
              help = "optimisation log (whole population) file name [default is %default]"),

  make_option(opt_str = c("-o", "--out"), type = "character", default = "Frontier.png",
              help = "output file name [default is %default]"),
  make_option(opt_str = c("-bs", "--base_size"), type = "numeric", default = 11,
              help = "base font size [default is %default]"),
  make_option(opt_str = c("-us", "--out_size"), type = "numeric", default = 10,
              help = "output size [default is %default]")
);

ArgumentsParser = OptionParser(option_list = ArgumentsList)
Arguments = parse_args(ArgumentsParser);

# Arguments$frontier = "Editing_demo_Edit0_Frontier.txt"
# Arguments$targets = "Editing_demo_Edit0_Targets.txt"
# Arguments$log = "Editing_demo_Edit0_OptimisationLogModeMaxCriterion.txt"

# ---- Read in the data & processing ----

# Baseline frontier
if (!is.null(Arguments$frontier)) {
  Frontier <- read.table(file = Arguments$frontier, header = TRUE)
  Frontier <- arrange(Frontier, FrontierDegree)
  Frontier$TargetDegree     <- round(Frontier$FrontierDegree)
  Frontier$TargetDegreeChar <- as.factor(paste(Frontier$TargetDegree, "°", sep = ""))
}

# Targets
if (!is.null(Arguments$targets)) {
  Targets <- read.table(file = Arguments$targets, header = TRUE)
  Targets <- arrange(Targets, FrontierDegree)
  Targets$TargetDegree     <- round(Targets$FrontierDegree)
  Targets$TargetDegreeChar <- as.factor(paste(Targets$TargetDegree, "°", sep = ""))

  if (is.null(Arguments$frontier)) {
    Frontier = Targets
  } else {
    colnames(Targets) = colnames(Frontier)
    Frontier = rbind(Frontier, Targets)
    Frontier = arrange(Frontier, FrontierDegree) %>%
      distinct(ModeOrPoint, .keep_all = TRUE)
  }
}

# Limits
MaxSelCriterPctLimits = c(0, 100)
MinCoancestryPctLimits = c(0, 100)

# Mapping formulae
# plot(Frontier$SelIntensity   ~ Frontier$MaxSelCriterPct)
fit <- coef(lm(Frontier$SelIntensity   ~ Frontier$MaxSelCriterPct))
SelIntFormula <- as.formula(paste("~", fit[1], "+ . *", fit[2]))

# plot(Frontier$CoancestryRate ~ Frontier$MinCoancestryPct)
fit <- coef(lm(Frontier$CoancestryRate ~ Frontier$MinCoancestryPct))
RatCoaFormula <- as.formula(paste("~", fit[1], "+ . *", fit[2]))

# Additional frontiers (must appear after setting the initial limits!!!)
FrontierAdd = vector(mode = "list", length = 4)
for (Tmp in 2:5) {
  Arg = paste("frontier", Tmp, sep = "")
  if (!is.null(Arguments[[Arg]])) {
    FrontierAdd[[Tmp - 1]] <- read.table(file = Arguments[[Arg]], header = TRUE)
    FrontierAdd[[Tmp - 1]] <- arrange(FrontierAdd[[Tmp - 1]], FrontierDegree)
    FrontierAdd[[Tmp - 1]]$MaxSelCriterPct = SelIntensity2MaxCriterionPct(SelIntensity = FrontierAdd[[Tmp - 1]]$SelIntensity,
                                                                          MinSelIntensity = filter(Frontier, ModeOrPoint == "ModeMinCoancestry")$SelIntensity,
                                                                          MaxSelIntensity = filter(Frontier, ModeOrPoint == "ModeMaxCriterion")$SelIntensity)
    MaxSelCriterPctLimits = range(c(MaxSelCriterPctLimits, range(FrontierAdd[[Tmp - 1]]$MaxSelCriterPct)))
    FrontierAdd[[Tmp - 1]]$MinCoancestryPct = CoancestryRate2MinCoancestryPct(CoancestryRate = FrontierAdd[[Tmp - 1]]$CoancestryRate,
                                                                              MinCoancestryRate = filter(Frontier, ModeOrPoint == "ModeMinCoancestry")$CoancestryRate,
                                                                              MaxCoancestryRate = filter(Frontier, ModeOrPoint == "ModeMaxCriterion")$CoancestryRate)
    MinCoancestryPctLimits = range(c(MinCoancestryPctLimits, range(FrontierAdd[[Tmp - 1]]$MinCoancestryPct)))
  }
}

# Optimisation path
if (!is.null(Arguments$log)) {
  Path <- read.table(file = Arguments$log, header = TRUE)
  Path$MaxSelCriterPct = SelIntensity2MaxCriterionPct(SelIntensity = Path$SelIntensity,
                                                      MinSelIntensity = filter(Frontier, ModeOrPoint == "ModeMinCoancestry")$SelIntensity,
                                                      MaxSelIntensity = filter(Frontier, ModeOrPoint == "ModeMaxCriterion")$SelIntensity)
  Path$MinCoancestryPct = CoancestryRate2MinCoancestryPct(CoancestryRate = Path$CoancestryRate,
                                                          MinCoancestryRate = filter(Frontier, ModeOrPoint == "ModeMinCoancestry")$CoancestryRate,
                                                          MaxCoancestryRate = filter(Frontier, ModeOrPoint == "ModeMaxCriterion")$CoancestryRate)
}

PathAdd = vector(mode = "list", length = 10)
for (Tmp in 2:11) {
  Arg = paste("log", Tmp, sep = "")
  if (!is.null(Arguments[[Arg]])) {
    PathAdd[[Tmp - 1]] = read.table(file = Arguments[[Arg]], header = TRUE)
    PathAdd[[Tmp - 1]]$MaxSelCriterPct = SelIntensity2MaxCriterionPct(SelIntensity = PathAdd[[Tmp - 1]]$SelIntensity,
                                                                      MinSelIntensity = filter(Frontier, ModeOrPoint == "ModeMinCoancestry")$SelIntensity,
                                                                      MaxSelIntensity = filter(Frontier, ModeOrPoint == "ModeMaxCriterion")$SelIntensity)
    PathAdd[[Tmp - 1]]$MinCoancestryPct = CoancestryRate2MinCoancestryPct(CoancestryRate = PathAdd[[Tmp - 1]]$CoancestryRate,
                                                                          MinCoancestryRate = filter(Frontier, ModeOrPoint == "ModeMinCoancestry")$CoancestryRate,
                                                                          MaxCoancestryRate = filter(Frontier, ModeOrPoint == "ModeMaxCriterion")$CoancestryRate)
  }
}

# Optimisation log (whole population)
if (!is.null(Arguments$logpop)) {
  PopPath <- read.table(file = Arguments$logpop, header = TRUE)
}

# ---- Basic frontier/targets plot ----

BaseSize <- Arguments$base_size
theme  <- theme_bw(base_size = BaseSize)
GeomTextSize <- BaseSize / (14 / 5)

# Roslin colors
ColRoslinViolet <- rgb(red = 193, green =  23, blue = 115, maxColorValue = 255) ## violet
ColRoslinBlue   <- rgb(red =  55, green = 152, blue = 217, maxColorValue = 255) ## blue
ColRoslinGreen  <- rgb(red =  95, green = 168, blue =  59, maxColorValue = 255) ## green
ColRoslinOrange <- rgb(red = 219, green =  79, blue =  15, maxColorValue = 255) ## orange

p <- ggplot(data = Frontier,
            mapping = aes(x = MinCoancestryPct, y = MaxSelCriterPct)) +
  geom_vline(xintercept = 0,   linetype = 2, color = "grey20") +
  geom_vline(xintercept = 100, linetype = 2, color = "grey20") +
  geom_hline(yintercept = 0,   linetype = 2, color = "grey20") +
  geom_hline(yintercept = 100, linetype = 2, color = "grey20") +
  geom_segment(mapping = aes(x = 0, y = 0, xend = MinCoancestryPct, yend = MaxSelCriterPct),
               data = filter(Frontier, FrontierDegree > 0, FrontierDegree < 90),
               color = "grey92", linetype = 2) +
  scale_y_continuous(limits = MaxSelCriterPctLimits, name = "Maximum gain (%)",
                     sec.axis = sec_axis(trans = SelIntFormula,
                                         name = "Selection intensity")) +
  scale_x_continuous(limits = MinCoancestryPctLimits, name = "Minimum coancestry (%)",
                     sec.axis = sec_axis(trans = RatCoaFormula,
                                         name = "Rate of coancestry")) +
  geom_line() +
  theme

# Additional frontiers
for (Tmp in 2:5) {
  Arg = paste("frontier", Tmp, sep = "")
  if (!is.null(Arguments[[Arg]])) {
    p = p + geom_line(data = FrontierAdd[[Tmp - 1]],
                      mapping = aes(x = MinCoancestryPct, y = MaxSelCriterPct))
  }
}

# ---- Add degree labels ----

# ... slightly away from the frontier
for (d in Frontier$TargetDegree) {
  # d <- 0
  x <- filter(Frontier, round(FrontierDegree) == d)
  if (nrow(x) > 0) {
    x <- rbind(x, x)
    x[2, "MaxSelCriterPct"]  <- 0
    x[2, "MinCoancestryPct"] <- 0
    # plot(x$MaxSelCriterPct ~ x$MinCoancestryPct)
    fit <- coef(lm(x$MaxSelCriterPct ~ x$MinCoancestryPct))
    if (d == 0) {
      x[2, "MaxSelCriterPct"] <- 0.95 * x[1, "MaxSelCriterPct"]
    } else {
      x[2, "MinCoancestryPct"] <- 0.95 * x[1, "MinCoancestryPct"]
      x[2, "MaxSelCriterPct"] <- fit[1] + fit[2] * x[2, "MinCoancestryPct"]
    }
    x <- x[2, ]
    p <- p + geom_label(mapping = aes(label = TargetDegreeChar), fill = "white", color = "grey20",
                        size = GeomTextSize, label.padding = unit(0, "lines"), label.size = NA,
                        data = x, show.legend = FALSE)
  }
}

# ---- Optimisation path (whole population) ----

if (!is.null(Arguments$logpop)) {
p <- p + geom_point(data = PopPath, show.legend = FALSE, color = ColRoslinViolet, size = .1, alpha = .1,
                    mapping = aes(x = MinCoancestryPct, y = MaxSelCriterPct))
}

# ---- Optimisation path ----

if (!is.null(Arguments$log)) {
  print(Path)
  p <- p + geom_path(data = Path, show.legend = FALSE, color = ColRoslinBlue,
                     mapping = aes(x = MinCoancestryPct, y = MaxSelCriterPct))
}
for (Tmp in 2:11) {
  Arg = paste("log", Tmp, sep = "")
  if (!is.null(Arguments[[Arg]])) {
    print(PathAdd[[Tmp - 1]])
    p <- p + geom_path(data = PathAdd[[Tmp - 1]], show.legend = FALSE, color = ColRoslinBlue,
                       mapping = aes(x = MinCoancestryPct, y = MaxSelCriterPct))
  }
}

# ---- Save plot ----

ggsave(plot = p, filename = Arguments$out, width = Arguments$out_size, height = Arguments$out_size, unit = "cm")

q()

