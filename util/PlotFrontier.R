#!/usr/bin/env Rscript

# TODO: add some options
#       - plots the paths, just the min and opt one or all of them
#       - handling only x limits and only y limits
#       - handle several frontiers
#       - plot the population "swarm"?

library(package="methods")
library(package="utils")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  FindYLim <- TRUE
  FindXLim <- TRUE
  ylim <- xlim <- NULL
  PlotPaths <- TRUE
} else {
  FindYLim <- FALSE
  FindXLim <- FALSE
  ylim <- as.numeric(c(args[1], args[2]))
  xlim <- as.numeric(c(args[3], args[4]))
  PlotPaths <- as.logical(args[5])
}

ColRoslinBlue   <- rgb(red= 55, green=152, blue=217, maxColorValue=255) ## blue
ColRoslinViolet <- rgb(red=193, green= 23, blue=115, maxColorValue=255) ## violet
ColRoslinGreen  <- rgb(red= 95, green=168, blue= 59, maxColorValue=255) ## green
ColRoslinOrange <- rgb(red=219, green= 79, blue= 15, maxColorValue=255) ## orange
ColRoslinGray   <- rgb(red= 83, green= 83, blue= 83, maxColorValue=255) ## orange

LogFiles <- dir(pattern="OptimisationLog")
if (length(LogFiles) < 1 & (!file.exists("Frontier.txt") | !file.exists("Targets.txt"))) {
  stop("ERROR: No optimisation and frontier or target log files to plot!")
}

pdf(file="Frontier.pdf", width=10, height=10*2/3, pointsize=14)

if (length(LogFiles) > 0) {
  ## Get the data
  LogFileCount <- 0
  Dat <- vector(mode="list", length=length(LogFiles))
  for (LogFile in LogFiles) {
    ## LogFile <- "OptimisationLogModeMinCoancestry.txt"
    ## LogFile <- "OptimisationLogModeMaxCriterion.txt"
    LogFileCount <- LogFileCount + 1
    Dat[[LogFileCount]] <- read.table(file=LogFile, header=TRUE)
    #Dat[[LogFileCount]] <- Dat[[LogFileCount]][order(Dat[[LogFileCount]]$CoancestryRate), ]
    if (FindYLim) {
      ylim <- range(c(ylim, Dat[[LogFileCount]]$SelIntensity))
    }
    if (FindXLim) {
      xlim <- range(c(xlim, Dat[[LogFileCount]]$CoancestryRate))
    }
  }

  ## Plot the data
  LogFileCount <- 0
  for (LogFile in LogFiles) {
    ## LogFile <- "OptimisationLogModeMinCoancestry.txt"
    ## LogFile <- "OptimisationLogOptimumGain.txt"
    LogFileCount <- LogFileCount + 1

    if (LogFile == "OptimisationLogModeMinCoancestry.txt") {
      Col <- ColRoslinBlue
    } else if (grepl(pattern="OptimisationLogModeOptTarget",x=LogFile)) {
      Col <- ColRoslinGreen
    } else if (LogFile == "OptimisationLogModeMaxCriterion.txt") {
      Col <- ColRoslinOrange
    } else if (LogFile == "OptimisationLogModeRan.txt") {
      Col <- ColRoslinViolet
    } else {
      Col <- ColRoslinGray
    }
    if (LogFile == "OptimisationLogModeRan.txt") {
      Sel <- which.max(Dat[[LogFileCount]]$Step)
      Ratio <- 2
      Type <- "o"
    } else {
      # Option to plot less points
      Sel <- rep(TRUE, times=nrow(Dat[[LogFileCount]]))
      Test <- which.max(Dat[[LogFileCount]]$Iteration)
      Ratio <- Dat[[LogFileCount]][Sel, "Objective"] / Dat[[LogFileCount]][Sel, "Objective"][Test]
      Ratio <- Ratio + min(Ratio)
      Ratio <- Ratio / max(Ratio)
      Ratio <- 1
      Type <- "b"
    }
    if (LogFileCount == 1) {
      plot(y=Dat[[LogFileCount]][Sel, "SelIntensity"], x=Dat[[LogFileCount]][Sel, "CoancestryRate"], type=Type,
           pch=21, lwd=0.5, ylim=ylim, xlim=xlim, cex=0.5*Ratio,
           col=ifelse(PlotPaths, Col, "white"), bg=ifelse(PlotPaths, Col, "white"),
           xlab="Rate of coancestry", ylab="Selection intensity")
    } else {
      if (PlotPaths) {
        points(y=Dat[[LogFileCount]][Sel, "SelIntensity"], x=Dat[[LogFileCount]][Sel, "CoancestryRate"], type=Type,
               pch=21, lwd=0.5, ylim=ylim, xlim=xlim, col=Col, bg=Col, cex=0.5*Ratio)
      }
    }
    if (LogFile %in% c("OptimisationLogModeMinCoancestry.txt", "OptimisationLogModeMaxCriterion.txt")) {
      abline(h=Dat[[LogFileCount]][Test, ]$SelIntensity,  lwd=1, lty=2, col=Col)
      abline(v=Dat[[LogFileCount]][Test, ]$CoancestryRate, lwd=1, lty=2, col=Col)
    }
    #readline(prompt=cat("Pause:", LogFile))
  }
}

if (file.exists("Frontier.txt") | file.exists("Targets.txt")) {
  if (file.exists("Frontier.txt")) {
    Frontier <- read.table(file="Frontier.txt", header=TRUE)
    Frontier <- Frontier[order(Frontier$CoancestryRate), ]
  }
  if (file.exists("Targets.txt")) {
    Targets <- read.table(file="Targets.txt", header=TRUE)
    Targets <- Targets[order(Targets$CoancestryRate), ]
    if (file.exists("Frontier.txt")) {
      names(Targets)[1] <- names(Frontier)[1]
      Frontier <- rbind(Frontier, Targets)
      Frontier <- Frontier[order(Frontier$CoancestryRate), ]
    } else {
      Frontier <- Targets
    }
  }
  # Tmp <- installed.packages()
  # if (!("cobs" %in% Tmp[, "Package"])) {
  #   install.packages(pkg="cobs")
  # }
  # if ("cobs" %in% Tmp[, "Package"]) {
  #   library(package="cobs")
  #   Mat <- matrix(nrow=nrow(Frontier),ncol=3)
  #   Mat[,1] <- 1 # fitted value will be >= observed value (0 for =)
  #   Mat[,2] <- Frontier$CoancestryRate
  #   Mat[,3] <- Frontier$SelIntensity
  #   Tmp <- cobs(x=Frontier$CoancestryRate, y=Frontier$SelIntensity, pointwise=Mat, ic="BIC")#, nknots=10)
  #   x <- seq(from=min(Frontier$CoancestryRate),
  #            to=max(Frontier$CoancestryRate),
  #            by=0.0001)
  #   y <- predict(Tmp, z=x)
  #   lines(y=y[, 2], x=x, pch=21, cex=.5, lwd=2, col=ColRoslinGray)
  #   #lines(y=yPAGE[, 2], x=x[-length(x)], pch=21, cex=.5, lwd=2, col=ColRoslinGray, lty=2)
  # } else {
    lines(y=Frontier$SelIntensity, x=Frontier$CoancestryRate,
          pch=21, cex=.5, lwd=2, col=ColRoslinGray)
  # }
}

dev.off()
