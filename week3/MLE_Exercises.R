# week 3: MLE exercises

library(devtools)
# make sure you have the latest version of dagdata 
#install_github("genomicsclass/dagdata")

library(dagdata)
# load the palindrome data from the Human cytomegalovirus genome
data(hcmv)

library(rafalib)
mypar()
plot(locations,rep(1,length(locations)), ylab="", yaxt="n") #yaxt="n" removes ticks and numbers on y-axis

breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))

hist(counts) # should follow a Poisson distribution => The distribution seems about right

# likelihood "L" of observing the data with lambda "lam"
# L(lam) = Pr(X_1 = x_1 AND X_2 = x_2 AND ... X_n = x_n; lam)
# now we assume that all X are independent, thus:
# L(lam) = Pr(X_1 = x_1; lam) * Pr(X_2 = x_2; lam) * ... * Pr(X_n = x_n; lam)

# example for lambda "lam" = 4
probs <- dpois(x = counts, lambda = 4)
likelihood <- prod(probs)
print(likelihood)
# observation: note that this is a tiny number. It is usually more convenient to compute log-likelihoods
logprobs <- dpois(x = counts, lambda = 4, log = T)
loglikelihood <- sum(logprobs)
print(loglikelihood)

## exercise 1
lambdas = seq(0,15,len=300)

# function to compute LLs given a list of lambdas and our data (i.e. counts)
compute_LL <- function(lambda, counts) {
  logprobs <- dpois(x = counts, lambda = lambda, log = T)
  LL <- sum(logprobs)
  return(LL)
}

LLs <- sapply(X = lambdas, FUN = compute_LL, counts = counts) #sapply return vector, lapply return list
plot(lambdas, LLs, type = "l")
max_idx <- which( max(LLs) == LLs)
opt_lambda <- lambdas[max_idx]
print(opt_lambda)
# make vertical line at optimal lambda
abline(v = opt_lambda, col="red")

## exercise 2

breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
binLocation=(breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation,counts,type="l",xlab="")
max_idx_count <- which.max(counts)
center_of_bin_with_highest_count <- binLocation[max_idx_count]
print(center_of_bin_with_highest_count)

## exercise 3

max_count <- max(counts)
print(max_count)

## exercise 4

# calculate probability of seeing value as extreme as our max count of 14
# compute lambda using the approximatoin mean(data)
lambda_ = mean(counts[-max_idx_count]) #exclude the max value of 14

# Pr(count >= 14) = 1 - Pr(count <= 13)
p_ex4 <- 1 - ppois(q = max_count-1, lambda = lambda_)
# hint: this probability can be interpreted as a "p-value"
print(p_ex4)

## exercise 5

# answer: We selected the highest region out of 57 and need to adjust for multiple testing.

## exercise 6

# for bonferroni we have to divide pvalue cut-off by number of comparisons
alpha <- 0.05
n_comparisons = length(counts)
pval_cutoff <- alpha/n_comparisons
print(pval_cutoff)

## exercise 7

ps <- (seq(along=counts) - 0.5)/length(counts) # define the quantiles
lambda <- mean( counts[-which.max(counts)])
poisq <- qpois(p = ps, lambda = lambda)
qqplot(poisq, counts)
abline(0,1)
