# exercises week 2: Bonferroni Correction

### EXERCISE #1
alphas <- seq(0,0.25,0.01)
ms <- c(2,10,100,1000)

# setings multi panel plot
par(mfrow=c(2,2)) #first argument (#rows) and second argument (#columns)

for (m in ms) {
  bonf <- alphas/m
  sidak <- 1 - ( 1 - alphas)^(1/m)
  # plot bonf in RED
  plot(alphas, bonf, col = "red", type = "l", lty = 1,
       xlab = "alpha", ylab = "p-value cutoff", main = paste("m =",m))
  # plot sidak in blue
  lines(alphas, sidak, col = "blue", type = "l", lty = 2)
  legend("topleft", legend = c("bonferonni", "sidak"), col = c("red", "blue"), lty = 1:2, inset = 0.05)
}

### EXERCISE #2

# NOTE: FWER is the probability of rejecting at least one null hypothesis given that all null hypotheses are true

# we know that p-values have a uniform distribution when the null (H0) is true, so let's generate 8793 p-values given H0 is true

set.seed(1)
n_simulation <- 10000
n_test <- 8793
bonf_alpha <- 0.05/n_test
sidak_alpha <- 1 - (1- 0.05)^(1/n_test)

# 1) bonferroni
out_bonf <- replicate(n_simulation, {
  # get p values
  pvals <- runif(n_test, 0, 1)
  # determine if at least one p value is below adjusted significance level
  out <- any(pvals < bonf_alpha)
})

# calculate FWER for both adjustment methods
FWER_bonf <- mean(out_bonf)
# print results
FWER_bonf

### EXERCISE #3

# 2) sidak
out_sidak <- replicate(n_simulation, {
  # get p values
  pvals <- runif(n_test, 0, 1)
  # determine if at least one p value is below adjusted significance level
  out <- any(pvals < sidak_alpha)
})


# calculate FWER for both adjustment methods
FWER_sidak <- mean(out_sidak)

# print results
FWER_sidak

