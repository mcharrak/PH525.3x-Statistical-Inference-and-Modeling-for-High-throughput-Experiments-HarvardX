# exercises week 2 FDR Exercises

library(devtools)
library(rafalib)

# run the following two lines of code if you do not have the genomics data or the packages genefilter and qvalue
# install_github("genomicsclass/GSE5859Subset")
# BiocManager::install(c("genefilter", "qvalue"))

library(GSE5859Subset)
data(GSE5859Subset)

### EXERCISE 1

groups <- factor(sampleInfo$group)
res <- rowttests(geneExpression, groups)
pvals <- res$p.value
count_ex1 <- sum(pvals < 0.05)
print(count_ex1)

### EXERCISE 2

# define bonf alpha
n_test <- dim(geneExpression)[1]
bonf_alpha <- 0.05/n_test
count_ex2 <- sum(pvals < bonf_alpha)
print(count_ex2)

### EXERCISE 3

fdrs <- p.adjust(p = pvals, method = "fdr")
count_ex3 <- sum(fdrs < 0.05)
print(count_ex3)

### EXERCISE 4

res_ex4 <- qvalue(pvals)
qvals <- res_ex4$qvalues
count_ex4 <- sum(qvals < 0.05)
print(count_ex4)

### EXERCISE 5

pi0 <- res_ex4$pi0
print(pi0)

### EXERCISE 6

ratios <- qvals/fdrs
plot(ratios)
abline(h=res_ex4$pi0, col=2)

# extras/bonus
hist(pvals, breaks = seq(0,1,length.out = 21))
expectedfreq <- length(pvals)/20 #per bin
abline(h=expectedfreq*pi0,col=2,lty=2) #red line (col=2) and dotted (lty=2)

### EXERCISES 7 to 12

n <- 24 #number of samples (e.g. participants)
m <- 8793 #number of features (e.g. genes)
mat <- matrix(rnorm(n*m), nrow = m, ncol = n)

# now for 500 features, set a difference of 2 between group 1 (1-12) and group 2 (13-24)
delta <- 2
positives <- 500
mat[1:positives, 1:(n/2)] <- mat[1:positives, 1:(n/2)] + delta

# set ground truth vector indicating which features have nullhypothesis true
m1 <- 500
m0 <- m - m1
nullHypothesis <- c( rep(FALSE, m1), rep(TRUE, m0))

# now the null hypothesis is only true for 8793-500 genes. So we have m0=8293 and m1=500 such that m=m0+m1
set.seed(1)
# bonf correction to achieve FWER of 0.05
alpha <- 0.05
bonf_alpha <- alpha/m

B <- 1000
res_ex7 <- replicate(B, {
  # 1) create data
  mat <- matrix(rnorm(n*m), nrow = m, ncol = n)
  # now for 500 features, set a difference of 2 between group 1 (1-12) and group 2 (13-24)
  delta <- 2
  positives <- 500
  mat[1:positives, 1:(n/2)] <- mat[1:positives, 1:(n/2)] + delta
  # 2) calculate p-values
  # markers: first 12 elements in each row are group1 (controls) and next 12 elements in each row are group2 (treatments)
  g <- factor( c(rep(0,(n/2)), rep(1,(n/2))))
  pvals <- rowttests(mat, g)$p.value
  
  #1
  calls_bonf <-  (pvals <= bonf_alpha)
  #2
  adjusted_pvals <- p.adjust(pvals, method = "fdr")
  calls_padjust <- (adjusted_pvals <= alpha)
  #3
  qvalue_output <- qvalue(pvals)
  qvals <- qvalue_output$qvalues
  calls_qvals <- (qvals <= alpha)
  
  # return 3 tuple
  return(c(calls_bonf, calls_padjust, calls_qvals))
})

#### EXERCISE 7

# bonferonni => WE ONLY CONSIDER GENES WHERE H0 IS TRUE
# => (we can see this from our denominator being m0 (#count of genes where H0 (nullHypothesis is true)))
start <- 0*m+1+m1
end <- 1*m
calls_bonf <- res_ex7[start : end, ]
# false positives are the cases where H0=TRUE & calls=TRUE
bonf_FPs <- colSums(calls_bonf & nullHypothesis[-(1:positives)]) # FPs is here Vs, which represents false positives
bonf_FPs_div_m0 <- bonf_FPs/m0
bonf_FPR <- mean(bonf_FPs_div_m0)
print(bonf_FPR)


#### EXERCISE 8

start <- 0*m+1
end <- start+positives-1
calls_bonf <- res_ex7[start : end, ]
# false positivies are the cases where pval>0.05, therefore, the negation to calls_bonf which checked for pval <= 0.05
bonf_FNs <- colSums(!calls_bonf)
bonf_FNs_div_m1 <- bonf_FNs/m1
bonf_FNR <- mean(bonf_FNs_div_m1)
print(bonf_FNR)

#### EXERCISE 9

start <- 1*m+1+m1
end <- 2*m
calls_padjusted <- res_ex7[start : end, ]
# false positives are the cases where H0=TRUE & calls=TRUE
padjusted_FPs <- colSums(calls_padjusted & nullHypothesis[-(1:positives)]) # FPs is here Vs, which represents false positives
padjusted_FPs_div_m0 <- padjusted_FPs/m0
padjusted_FPR <- mean(padjusted_FPs_div_m0)
print(padjusted_FPR)
# => the larger m1 (cases where nullHypothesis is not true, the more conservative is this FDR approximation (i.e. the higher will the FDR approximation be))

#### EXERCISE 10

start <- 1*m+1
end <- start+positives-1
calls_padjusted <- res_ex7[start : end, ]
# false positivies are the cases where pval>0.05, therefore, the negation to calls_padjusted which checked for pval <= 0.05
padjusted_FNs <- colSums(!calls_padjusted)
padjusted_FNs_div_m1 <- padjusted_FNs/m1
padjusted_FNR <- mean(padjusted_FNs_div_m1)
print(padjusted_FNR)

#### EXERCISE 11

start <- 2*m+1+m1
end <- 3*m
calls_qvals <- res_ex7[start : end, ]
# false positives are the cases where H0=TRUE & calls=TRUE
qvals_FPs <- colSums(calls_qvals & nullHypothesis[-(1:positives)]) # FPs is here Vs, which represents false positives
qvals_FPs_div_m0 <- qvals_FPs/m0
qvals_FPR <- mean(qvals_FPs_div_m0)
print(qvals_FPR)
# => the larger m1 (cases where nullHypothesis is not true, the more conservative is this FDR approximation (i.e. the higher will the FDR approximation be))

#### EXERCISE 12

start <- 2*m+1
end <- start+positives-1
calls_qvals <- res_ex7[start : end, ]
# false positivies are the cases where pval>0.05, therefore, the negation to calls_qvals which checked for pval <= 0.05
qvals_FNs <- colSums(!calls_qvals)
qvals_FNs_div_m1 <- qvals_FNs/m1
qvals_FNR <- mean(qvals_FNs_div_m1)
print(qvals_FNR)
