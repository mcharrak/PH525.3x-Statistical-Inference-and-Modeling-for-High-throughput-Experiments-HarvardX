# notes lecture week 2

# load libraries
library(downloader)
library(genefilter)

### Error Rates and Procedures Examples

url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)

### 1) completely null model
set.seed(1)
population = unlist(read.csv("femaleControlsPopulation.csv"))

alpha <- 0.05
N <- 12
m <- 10000 # m is the number of experiments/tests

pvals <- replicate(m, {
  control <- sample(x = population, size = N)
  treatment <- sample(x = population, size = N)
  t.test(treatment,control)$p.value
})
sum(pvals < 0.05)

### 2) alternative hypothesis being true for 10%
N <- 12
m <- 10000
p0 <- 0.9 # 10% of the diets work, 90% don't work
m0 <- m*p0
m1 <- m - m0

nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1))
delta <- 3

# now we are ready to simulate 10,000 tests
calls <- sapply(1:m, function(i){
  control <- sample(x = population, size = N)
  treatment <- sample(x = population, size = N)
  # change 10% of the treatments s.t. alternative is true
  if( !nullHypothesis[i]) {treatment <- treatment + delta}
  ifelse( t.test(treatment, control)$p.value < alpha,
          "called significant",
          "not called significant")
})

# turn the results in calls into a factor for better readability
null_hypothesis <- factor(nullHypothesis,levels = c("TRUE","FALSE"))
table(null_hypothesis,calls)

# 3) repeat above simulation multiple times (i.e. B times)
B <- 10
system.time(
VandS <- replicate(B, {
  calls <- sapply(1:m, function(i){
    control <- sample(x = population, size = N)
    treatment <- sample(x = population, size = N)
    # change 10% of the treatments s.t. alternative is true
    if( !nullHypothesis[i]) {treatment <- treatment + delta}
    t.test(treatment, control)$p.value < alpha
  })
  cat("V =", sum(nullHypothesis & calls), "S =", sum(!nullHypothesis & calls),"\n") #note: sum(TRUE & FALSE) equals 1
  c(sum(nullHypothesis & calls), sum(!nullHypothesis & calls))
})
)

##### vectorizing the code from above because in R, 
##### operations based on matrices are typically much faster than operations performed within loops or sapply()

alpha <- 0.05
N <- 12
m <- 10000
p0 <- 0.90 ##10% of diets work, 90% don't
m0 <- m*p0
m1 <- m-m0
nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1))
delta <- 3

library(genefilter) ## rowttests is here
set.seed(1)

### define groups to be used with rowttests
g <- factor( c(rep(0,N), rep(1,N)))
B <- 10 ## number of simulations

system.time(
VandS <- replicate(B, {
  ## matrix with control data, (rows are tests, columns are mice (mäuse))
  controls <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
  
  ## matrix with control data, (rows are tests, columns are mice)
  treatments <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
  
  ## add effect to 10% of the rows in treatments matrix
  treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),] + delta
  
  ## combine to form one matrix
  dat <- cbind(controls, treatments)
  
  ## perform ttest for each row, using g, the factor which code the grouping to be tested
  calls <- rowttests(dat, g)$p.value < 0.05
  
  # count false positive (type-1-error) and false negatives (type-2-error)
  c(sum(nullHypothesis & calls), sum(!nullHypothesis & calls))
})
)

#### video lecture: False Discovery Rate and Benjamini–Hochberg procedure ####
controls <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)
treatments <- matrix(sample(population, size = N*m, replace = TRUE), nrow = m)

treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),] + delta

# combine to form one matrix
dat <- cbind(controls,treatments)

# markers: first 12 elements in each row are group1 (controls) and next 12 elements in each row are group2 (treatments)
g <- factor( c(rep(0,N), rep(1,N)))

pvals <- rowttests(dat, g)$p.value