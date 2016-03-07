#!/usr/bin/env Rscript

# TODO: add some options
#       - plots the paths, just the min and opt one or all of them
#       - handling only x limits and only y limits

library(package="methods")
library(package="utils")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  FindYLim <- TRUE
  FindXLim <- TRUE
  ylim <- xlim <- NULL
} else {
  FindYLim <- FALSE
  FindXLim <- FALSE
  ylim <- as.numeric(c(args[1], args[2]))
  xlim <- as.numeric(c(args[3], args[4]))
}

ColRoslinBlue   <- rgb(red= 55, green=152, blue=217, maxColorValue=255) ## blue
ColRoslinViolet <- rgb(red=193, green= 23, blue=115, maxColorValue=255) ## violet
ColRoslinGreen  <- rgb(red= 95, green=168, blue= 59, maxColorValue=255) ## green
ColRoslinOrange <- rgb(red=219, green= 79, blue= 15, maxColorValue=255) ## orange
ColRoslinGray   <- rgb(red= 83, green= 83, blue= 83, maxColorValue=255) ## orange

LogFileCount <- 0
LogFiles <- dir(pattern="OptimisationLog")
LogFiles <- LogFiles[!grepl(pattern="Initial.txt", x=LogFiles)]
if (length(LogFiles) < 1 & !file.exists("Frontier.txt")) {
  stop("ERROR: No optimisation and frontier log files to plot!")
}

pdf(file="Frontier.pdf")

if (length(LogFiles) > 0) {
  ## Get the data
  Dat <- vector(mode="list", length=length(LogFiles))
  for (LogFile in LogFiles) {
    ## LogFile <- "OptimisationLog1MinimumInbreeding.txt"
    ## LogFile <- "OptimisationLog2OptimumGain.txt"
    LogFileCount <- LogFileCount+1
    Dat[[LogFileCount]] <- read.table(file=LogFile, header=TRUE)
    colnames(Dat[[LogFileCount]]) <- c("Step", "AcceptRate", "Criterion", "Penalties", "Gain", "GainStand", "PopInbreed", "RatePopInb", "PopInbree2", "IndInbreed")
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
    ## LogFile <- "OptimisationLog1MinimumInbreeding.txt"
    ## LogFile <- "OptimisationLog2OptimumGain.txt"
    LogFileCount <- LogFileCount + 1
    
    if (LogFile == "OptimisationLog1MinimumInbreeding.txt") {
      Col <- ColRoslinBlue
    } else if (LogFile == "OptimisationLog2OptimumGain.txt") {
      Col <- ColRoslinViolet
    } else {
      Col <- ColRoslinGray
    }
    Test <- which.max(Dat[[LogFileCount]]$Step)
    if (LogFileCount == 1) {
      plot(y=Dat[[LogFileCount]]$GainStand, x=Dat[[LogFileCount]]$RatePopInb, type="o",
           pch=21, lwd=0.5, ylim=ylim, xlim=xlim, col=Col, bg=Col,
           cex=0.5*Dat[[LogFileCount]]$Criterion/Dat[[LogFileCount]][Test, "Criterion"],
           xlab="Rate of inbreeding", ylab="Genetic gain (standardized)")
    } else {
      points(y=Dat[[LogFileCount]]$GainStand, x=Dat[[LogFileCount]]$RatePopInb, type="o",
             pch=21, lwd=0.5, ylim=ylim, xlim=xlim, col=Col, bg=Col,
             cex=0.5*Dat[[LogFileCount]]$Criterion/Dat[[LogFileCount]][Test, "Criterion"])
    }
    if (LogFile == "OptimisationLog1MinimumInbreeding.txt") {
      DeltaF <- axis(side=1)
      CoefF <- Dat[[LogFileCount]][nrow(Dat[[LogFileCount]]), "PopInbreed"]
      CoefF <- DeltaF*(1-CoefF) + CoefF
      axis(side=3, at=DeltaF, labels=round(CoefF,digits=3))
      mtext(side=3, text="Coef. of inbreeding", line=2.5)
    }
    if (LogFile %in% c("OptimisationLog1MinimumInbreeding.txt", "OptimisationLog2OptimumGain.txt")) {
      abline(h=Dat[[LogFileCount]][Test, ]$GainStand,  lwd=1, lty=2, col=Col)
      abline(v=Dat[[LogFileCount]][Test, ]$RatePopInb, lwd=1, lty=2, col=Col)
    }
  }
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
    Tmp <- cobs(x=Frontier$RatePopInb, y=Frontier$GainStand, pointwise=Mat,
                ic="BIC")#, nknots=10)
    x <- seq(from=min(Frontier$RatePopInb),
             to=max(Frontier$RatePopInb),
             by=0.0001)
    y <- predict(Tmp, z=x)
    lines(y=y[, 2], x=x, pch=21, cex=.5, lwd=2, col=ColRoslinGray)
  } else {
    lines(y=Frontier$GainStand, x=Frontier$RatePopInb,
          pch=21, cex=.5, lwd=2, col=ColRoslinGray)
  }
}

dev.off()
