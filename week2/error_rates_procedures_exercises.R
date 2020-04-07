# exercises week 2 Error Rates and Procedures Exercises


#Q2
B <- 1000
# we reject H0 iff p-value smaller alpha (0.05)
B<-1000
minpval <- replicate(B, min(runif(8793,0,1))>0.05) #fetch smallest p-value and check if it is larger than alpha 
p <- 1- mean(minpval) # p(at least 1 rejection) =  1 - P(no rejection) ## we use the average of no rejections we make in these 1000 tests
p

#Q3: 0.05 must equal p(at least 1 rejection) = 1 - p(no rejection) = 1 - (1-alpha)^8793 --> solve for new alpha
## thus: alpha_new = 1 - 0.95^(1/8793)

B <- 10000
candidates <- 10^seq(-7,-4,0.1)

probs = sapply(candidates, function(candidate) {
  minpval = replicate(B, min(runif(8793,0,1))>candidate) #fetch smallest p-value and check if it is larger than alpha  
  1- mean(minpval) # p(at least 1 rejection) =  1 - P(no rejection) ## we use the average of no rejections we make in these 1000 tests
})
candidates[which.min(abs(probs-0.05))]

