setwd(dir="~/Gitbox/AlphaSuite/AlphaMate/test")
Ebv <- read.table(file="Ebv.txt")
Gender <- read.table(file="Gender.txt")

n <- nrow(Ebv)
Mat <- matrix(nrow=n, ncol=n)
for (i in 1:n) {
  for (j in 1:n) {
    Mat[i,j] <- runif(n=1)*(Ebv[i,2]+Ebv[j,2])/2
  }
}
Mat[1:10,1:10]
head(Ebv)

## !GenderMatters, SelfingAllowed
k <- 0
for (j in 1:n) {
  for (i in j:n) {
    if (k == 0) {
      Append <- FALSE
    } else {
      Append <- TRUE
    }
    cat(i, j, Mat[i,j], Mat[i,j]/2, "\n", file="MatingVal.txt", append=Append)
    k <- k + 1
  }
}

## !GenderMatters, !SelfingAllowed
k <- 0
for (j in 1:n) {
  for (i in j:n) {
    if (i > j) {
      if (k == 0) {
        Append <- FALSE
      } else {
        Append <- TRUE
      }
      cat(i, j, Mat[i,j], Mat[i,j]/2, "\n", file="MatingValNoSelfing.txt", append=Append)
      k <- k + 1
    }
  }
}

## GenderMatters, !SelfingAllowed
k <- 0
for (j in 1:n) {
  for (i in 1:n) {
    if (i != j) {
      if (i %in% Gender[Gender[, 2] == 1, 1]) {
        if (j %in% Gender[Gender[, 2] == 2, 1]) {
          if (k == 0) {
            Append <- FALSE
          } else {
            Append <- TRUE
          }
          cat(i, j, Mat[i,j], Mat[i,j]/2, "\n", file="MatingValGender.txt", append=Append)
          k <- k + 1
        }
      }
    }
  }
}
