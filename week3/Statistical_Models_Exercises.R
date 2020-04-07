# week 3: Statistical Models Exercises

## exercise 1

# Pr(X_i = 1) = p = Pr(ball = red)
# Pr(X_i = 0) = 1- p = Pr(ball = blue)

# binomial distribution: Pr(S_N = k) = (N over k) p^k (1-p)^(N-k) = Pr(k red balls after N bernoulli trials)
# in R we can use the function "pbinom() to calculate Pr(S_N <= k)

N_kids=4
p_girl=0.49

#Pr(2 girls AND 2 boys) = Pr(S_N = 2)
p_ex1 <- pbinom(2, size = 4, prob = p_girl) - pbinom(1, size = 4, prob = p_girl)
p_ex1_alternative <- dbinom(x = 2, size = 4, prob = p_girl)
print(p_ex1)
print(p_ex1_alternative)

## exercise 2

#Pr(4 girls AND 6 boys) = Pr(S_N = 4)
p_ex2 <- dbinom(x = 4, size = 10, prob = p_girl)
print(p_ex2)

## exercise 3

# Pr(base == C) = Pr(base == G) = 0.2
p_C_or_G=0.4
N_bases=20
k=10
p_ex3 <- 1 - pbinom(q = k, size = N_bases, prob = p_C_or_G)
print(p_ex3)

## exercise 4
p_win=1/175223510
N_tickets=189000000 
## in the question it is mentioned that tickets are generated with replacement, this means p_win does not change
## (e.g. in real world, once a winning ticket was change, the probability of generating a second winning ticket gets smaller
## because we would have exactly 1 winning ticket less)

#Pr(S_N >= 1) = 1 - Pr(S_N = 0)
p_ex4 <- 1 - pbinom(q = 0, size = N_tickets, prob = p_win)
# alternatively we can use dbinom because for k=0 pbinom==dbinom see:
p_ex4_alternative <- 1 - dbinom(x = 0, size = N_tickets, prob = p_win)
print(p_ex4)
print(p_ex4_alternative)

## exercise 5

#Pr(S_N >= 2) = 1 - Pr(S_N <= 1)
p_ex5 <- 1 - pbinom(q = 1, size = N_tickets, prob = p_win)
print(p_ex5)

## exercise 6

lower_bound <- 0.35*20 #7 (excluded)
upper_bound <- 0.45*20 #9 (included)

# Pr(lower_bound < S_N <= upper_bound) = Pr(7 < S_N <= 9) = Pr(S_N = 8) + Pr(S_N = 9) = Pr(S_N <= 9) - Pr(S_N <= 7)
p_ex6 <- dbinom(x = 8, size = N_bases, prob = p_C_or_G) + dbinom(x = 9, size = N_bases, prob = p_C_or_G)
print(p_ex6)
# altertnative method using pbinom
p_ex6_alternative <- pbinom(q = upper_bound, size = N_bases, prob = p_C_or_G) - pbinom(q = lower_bound, size = N_bases, prob = p_C_or_G)
print(p_ex6_alternative)

## exercise 7

# E[S_N] = N*p
mu_S_N = N_bases * p_C_or_G
# Var[S_N] = N*p*(1-p)  
var_S_N = N_bases * p_C_or_G * (1-p_C_or_G)
sd_S_N = sqrt(var_S_N)

# compute Z-score for lower and upper bound
z_lower = (lower_bound - mu_S_N)/sd_S_N
z_upper = (upper_bound - mu_S_N)/sd_S_N

p_ex7 <- pnorm(q = z_upper) - pnorm(q = z_lower)
print(p_ex7)


## exercise 8 (repeat exercise 3 experiments but with N_bases = 1000)

N_bases_new=1000

lowerBound=0.35*N_bases_new
upperBound=0.45*N_bases_new


# 1) normal approximation probability
mu_ = N_bases_new * p_C_or_G
var_ = N_bases_new * p_C_or_G * (1 - p_C_or_G)
sd_ = sqrt(var_)
zUpper = (upperBound - mu_)/sd_
zLower = (lowerBound - mu_)/sd_
p_approx <-  pnorm(q = zUpper) - pnorm(q = zLower)

# 2) binomial exact probability
p_exact <- pbinom(q = upperBound, size = N_bases_new, prob = p_C_or_G) - pbinom(q = lowerBound, size = N_bases_new, prob = p_C_or_G)

# finally calculate difference
p_diff <- abs(p_exact - p_approx)
print(p_diff)
# observation from result p_diff: we see that as N becomes larger (20 -> 1000), the approximation gets close to the exact probability
# CONCLUSION: THE BINOMIAL DIST IS APPROX. NORMAL WITH LARGE N AND p NOT TOO CLOSE TO p=0 OR p=1

## exercise 9


Ns <- c(5,10,30,100) # number of pieces that contain C
ps <- c(0.01,0.10,0.5,0.9,0.99) # proportion of pieces that have methylated C
library(rafalib)
mypar(4,5)

for (N in Ns) {
  ks <- 1:(N-1)
  for (p in ps) {
    p_exact = dbinom(x = ks, size = N, prob = p)
    
    # because k is an integer
    upperBounds <- ks+0.5
    lowerBounds <- ks-0.5
    
    mu_ = N*p
    sd_ = sqrt(N*p*(1-p))
      
    upperZs <- (upperBounds - mu_)/sd_
    lowerZs <- (lowerBounds - mu_)/sd_
    p_approx = pnorm(upperZs) - pnorm(lowerZs)
    
    LIM <- range(c(p_approx, p_exact))
    
    plot(p_exact, p_approx, main=paste("N=",N," p =",p), xlim=LIM, ylim=LIM, col=1, pch=16)
    # add line to better see deviation beteween exact and approx probabilities for different k-values
    abline(0,1)
  }
}

# observation: when p is very small, the normal approximation for the poisson distribution is bad. 
# however, if N is large, the Poission approximation becomes good

## exercise 10

N <- 189000000
p <- 1/175223510
# 1) exact: prob of 2 people winning the lottery
p_2win_exact <- dbinom(2,N,p)
print(p_2win_exact) 
# 2) normal approx: prob of 2 people winning the lottery
a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
p_2win_approx <- pnorm(a) - pnorm(b)
print(p_2win_approx)
# observation: from p_2win_approx we can see that using the normal approximation we highly overestimate the probability of 2 people winning the lottery
# therefore: let us use the Poission approximation here in order to approximate the binomial distributed probability

# define lambda == #tickets per 189,000,000 that win the lottery, which is simply N_tickets*p_win
lambda_ = N_tickets * p_win
# 3) poisson approx: prob of 2 people winning the lottery
p_2win_poisson <- dpois(x = 2, lambda = lambda_)
print(p_2win_poisson)
# observation: possion makes a good approximation when N is large, here N is very very large,
# so the exact binomial prob is practically the same as the approximation using poisson

# Pr(S_N >= 2) = 1 - Pr(S_N <= 1)
p_ex10 <- 1 - ppois(q = 1, lambda = lambda_)
print(p_ex10)