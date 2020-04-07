# notes lecture on: False Discovery Rate and Benjamini-Hochberg Procedure

library(downloader)
library(genefilter)

set.seed(1)


url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)

population = unlist( read.csv("femaleControlsPopulation.csv") )
set.seed(1)
alpha <- 0.05
N <- 12
m <- 10000
p0 <- 0.90 ##10% of diets work, 90% don't
m0 <- m*p0
m1 <- m-m0
nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1))
delta <- 3

controls <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
treatments <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)

treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),] + delta

# combine to form one matrix
dat <- cbind(controls,treatments)

# markers: first 12 elements in each row are group1 (controls) and next 12 elements in each row are group2 (treatments)
g <- factor( c(rep(0,N), rep(1,N)))

pvals <- rowttests(dat, g)$p.value

# count false discovery when using bonferroni correction
sum(pvals <= alpha/m)
# we see that out of the 10% (1000) diets that work only 2 get called significant -> means very high false negative rate
## --> this is one of the issues with the bonferroni correction: very high false negative rate 
##     because the correction is too conservative/strict

# FDR

calls <- pvals <= alpha
R <- sum(calls)
V <- sum(nullHypothesis & calls) # find falsely identified significant genes -> false positives
Q <- ifelse(R>0, V/R, 0) # Q is the fraction of tests which we call significant but actually are false positive

# Q is a random variable, it is not the FDR, but the expectation of Q is the FDR which we can get by monte carlo simulation
## => E[Q] = FDR

B <- 1000
Qs <- replicate(B, {
  controls <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
  treatments <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
  
  treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),] + delta
  
  # combine to form one matrix
  dat <- cbind(controls,treatments)
  
  # markers: first 12 elements in each row are group1 (controls) and next 12 elements in each row are group2 (treatments)
  g <- factor( c(rep(0,N), rep(1,N)))
  
  pvals <- rowttests(dat, g)$p.value
  
  calls <- pvals <= alpha
  R <- sum(calls)
  V <- sum(nullHypothesis & calls) # find falsely identified significant genes -> false positives
  Q <- ifelse(R>0, V/R, 0) # Q is the fraction of tests which we call significant but actually are false positive
  return(Q)
})

FDR <- mean(Qs)
print(FDR)

# controlling FDR

library(rafalib)
mypar(1,1)
hist(Qs)

fdr <- p.adjust(pvals, method = "fdr")
mypar(1,1)
plot(x = pvals, fdr, log="xy")

alpha <- 0.05
B <- 100 ## number of simulatoins. we should increase for more precision.
res <- replicate(B, {
  controls <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
  treatments <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
  
  treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),] + delta
  
  # combine to form one matrix
  dat <- cbind(controls,treatments)
  
  # markers: first 12 elements in each row are group1 (controls) and next 12 elements in each row are group2 (treatments)
  g <- factor( c(rep(0,N), rep(1,N)))
  
  pvals <- rowttests(dat, g)$p.value
  adjusted_pvals <- p.adjust(pvals, method = "fdr")
  
  calls <- adjusted_pvals <= alpha
  R <- sum(calls)
  V <- sum(nullHypothesis & calls) # find falsely identified significant genes -> false positives
  Q <- ifelse(R>0, V/R, 0) # Q is the fraction of tests which we call significant but actually are false positive
  return(c(R,Q))
})

Qs <- res[2,]
mypar(1,1)
hist(Qs)
FDR <- mean(Qs)
print(FDR)

#### the benjamini-hochberg procedure allows you to control the false discovery rate (FDR)
#### the bonferroni prodecure allows you to control the familywise error rate (FWER)
