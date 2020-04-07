# week1 exercise: Inference in Practice Exercises

set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)

# simulate 1000 pvalues
pvals <- replicate(n = 1000, {
  N = 12
  control = sample(population[,1], N)
  treatment = sample(population[,1], N)
  t.test(treatment,control)$p.value
})

head(pvals)
hist(pvals)

# Q1
total_freq <- 1000
below_alpha1_freq <-  sum(pvals <= 0.05)

prop1 <- below_alpha1_freq/total_freq
prop1

#Q2
below_alpha2_freq <-  sum(pvals <= 0.01)

prop2 <- below_alpha2_freq/total_freq
prop2

#Q3
set.seed(100)
n_diets <- 20
pvals_mice_diets <- replicate(n = n_diets,{
  cases = rnorm(10,30,2)
  controls = rnorm(10,30,2)
  t.test(cases,controls)$p.value
})
count_below_alpha3 <- sum(pvals_mice_diets <= 0.05)
count_below_alpha3

#Q4
set.seed(100)
numbers <- replicate(n = 1000, {
  curr_pvals <- replicate(n = n_diets,{
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.value
  })
  sum(curr_pvals <= 0.05)
})
numbers_avg <- mean(numbers)
numbers_avg

#Q5
set.seed(100)
binaries <- replicate(n = 1000, {
  curr_pvals <- replicate(n = n_diets,{
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.value
  })
  curr_count <- sum(curr_pvals <= 0.05)
  if (curr_count == 0) {
    return(0)
  } else {
    return(1)
  }
})
binaries
binaries_avg <- mean(binaries) # which is the proportion of rejecting the null hypothesis even though it is true
print(binaries_avg)
# alternatively:
#mean(numbers > 0)