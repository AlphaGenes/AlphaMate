#!/usr/bin/env Rscript

# TODO: add some options
#       - plots the paths, just the min and opt one or all of them
#       - handling only x limits and only y limits
#       - handle several frontiers
#       - plot the population "swarm"?

library(package="methods")
library(package="utils")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
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

LogFileCount <- 0
LogFiles <- dir(pattern="OptimisationLog")
if (length(LogFiles) < 1 & !file.exists("Frontier.txt")) {
  stop("ERROR: No optimisation and frontier log files to plot!")
}

pdf(file="Frontier.pdf", width=10, height=10*2/3, pointsize=14)

if (length(LogFiles) > 0) {
  ## Get the data
  Dat <- vector(mode="list", length=length(LogFiles))
  for (LogFile in LogFiles) {
    ## LogFile <- "OptimisationLogMinimumInbreeding.txt"
    ## LogFile <- "OptimisationLogOptimumGain.txt"
    LogFileCount <- LogFileCount+1
    Dat[[LogFileCount]] <- read.table(file=LogFile, header=TRUE)
    #Dat[[LogFileCount]] <- Dat[[LogFileCount]][order(Dat[[LogFileCount]]$RatePopInb), ]
    if (FindYLim) {
      ylim <- range(c(ylim, Dat[[LogFileCount]]$GainStand))
    }
    if (FindXLim) {
      xlim <- range(c(xlim, Dat[[LogFileCount]]$RatePopInb))
    }
  }

  ## Plot the data
  LogFileCount <- 0
  for (LogFile in LogFiles) {
    ## LogFile <- "OptimisationLogMinimumInbreeding.txt"
    ## LogFile <- "OptimisationLogOptimumGain.txt"
    LogFileCount <- LogFileCount + 1

    if (LogFile == "OptimisationLogMinimumInbreeding.txt") {
      Col <- ColRoslinBlue
    } else if (LogFile == "OptimisationLogRandomMating.txt") {
      Col <- ColRoslinOrange
    } else if (LogFile == "OptimisationLogOptimumGain.txt") {
      Col <- ColRoslinViolet
    } else {
      Col <- ColRoslinGray
    }
    if (LogFile == "OptimisationLogRandomMating.txt") {
      Sel <- which.max(Dat[[LogFileCount]]$Step)
      Ratio <- 2
      Type <- "o"
    } else {
      # Option to plot less points
      Sel <- rep(TRUE, times=nrow(Dat[[LogFileCount]]))
      Test <- which.max(Dat[[LogFileCount]]$Step)
      Ratio <- Dat[[LogFileCount]][Sel, "Criterion"] / Dat[[LogFileCount]][Sel, "Criterion"][Test]
      Ratio <- Ratio + min(Ratio)
      Ratio <- Ratio / max(Ratio)
      Type <- "l"
    }
    if (LogFileCount == 1) {
      plot(y=Dat[[LogFileCount]][Sel, "GainStand"], x=Dat[[LogFileCount]][Sel, "RatePopInb"], type=Type,
           pch=21, lwd=0.5, ylim=ylim, xlim=xlim, cex=0.5*Ratio,
           col=ifelse(PlotPaths, Col, "white"), bg=ifelse(PlotPaths, Col, "white"),
           xlab="Rate of inbreeding", ylab="Genetic gain (standardized)")
    } else {
      if (PlotPaths) {
        points(y=Dat[[LogFileCount]][Sel, "GainStand"], x=Dat[[LogFileCount]][Sel, "RatePopInb"], type=Type,
               pch=21, lwd=0.5, ylim=ylim, xlim=xlim, col=Col, bg=Col, cex=0.5*Ratio)
      }
    }
    if (LogFile %in% c("OptimisationLogMinimumInbreeding.txt", "OptimisationLogOptimumGain.txt")) {
      abline(h=Dat[[LogFileCount]][Test, ]$GainStand,  lwd=1, lty=2, col=Col)
      abline(v=Dat[[LogFileCount]][Test, ]$RatePopInb, lwd=1, lty=2, col=Col)
    }
    #readline(prompt=cat("Pause:", LogFile))
  }
  DeltaF <- axis(side=1)
  Tmp <- read.csv(file="ConstraintPopulationInbreeding.txt", header=FALSE)
  CoefF <- Tmp[Tmp$V1 == "Old_coancestry", "V2"]
  CoefF <- DeltaF*(1-CoefF) + CoefF
  axis(side=3, at=DeltaF, labels=round(CoefF,digits=3))
  mtext(side=3, text="Coef. of inbreeding", line=2.5)
}

if (file.exists("Frontier.txt")) {
  Frontier <- read.table(file="Frontier.txt", header=TRUE)
  Frontier <- Frontier[order(Frontier$RatePopInb),]
  Tmp <- installed.packages()
  if (!("cobs" %in% Tmp[, "Package"])) {
    install.packages(pkg="cobs")
  }
  if ("cobs" %in% Tmp[, "Package"]) {
    library(package="cobs")
    Mat <- matrix(nrow=nrow(Frontier),ncol=3)
    Mat[,1] <- 1 # fitted value will be >= observed value (0 for =)
    Mat[,2] <- Frontier$RatePopInb
    Mat[,3] <- Frontier$GainStand
    Tmp <- cobs(x=Frontier$RatePopInb, y=Frontier$GainStand, pointwise=Mat, ic="BIC")#, nknots=10)
    x <- seq(from=min(Frontier$RatePopInb),
             to=max(Frontier$RatePopInb),
             by=0.0001)
    y <- predict(Tmp, z=x)
    lines(y=y[, 2], x=x, pch=21, cex=.5, lwd=2, col=ColRoslinGray)
    lines(y=yPAGE[, 2], x=x[-length(x)], pch=21, cex=.5, lwd=2, col=ColRoslinGray, lty=2)
  } else {
    lines(y=Frontier$GainStand, x=Frontier$RatePopInb,
          pch=21, cex=.5, lwd=2, col=ColRoslinGray)
  }
}

dev.off()
