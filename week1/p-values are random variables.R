population = read.csv("/Users/Amine/Desktop/ML_classes/Stats_R_harvard/femaleControlsPopulation.csv")
population <- population[,1]
N <- 12
B <- 10^4
pvals <- replicate(B,{
  control = sample(population,N)
  treatment = sample(population,N)
  t.test(treatment,control)$p.val
})
# lines in curly brackets is the function that replicate re-uses
hist(pvals)
# the histogram shows that our experimental p-values (for the null distribution) follow a uniform distribution and the p-values are indeed a random variable (RV)
# mathematically it can be shown that under H0, p-values follow a uniform distribution
