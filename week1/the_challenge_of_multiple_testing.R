library(GSE5859Subset)
data(GSE5859Subset)

g <- sampleInfo$group
# choose 25th gene
e <- geneExpression[25,]

library(rafalib)

# qqplot to check if data is normal - checking assumption of a t-test
mypar(1,2)
qqnorm(e[g==1])
qqline(e[g==1])

qqnorm(e[g==0])
qqline(e[g==1])

group1 <- e[g==1]
group2 <- e[g==0]
t.test(group2,group1)

# repeat for all features/genes
mytest <- function(x){
  t.test(x[g==1],x[g==0],var.equal = TRUE)$p.value
}

mytest(geneExpression[25,])

# apply mytest function (t-test) on every gene
pvals <- apply(geneExpression, MARGIN = 1, mytest)
length(pvals)

# count pvals which are significant
sum(pvals <= 0.05)

# let us test what happens if the data is totally randomly distributed with all values from the same population
# then the null hypothesis needs to be true for every row of our new randomData matrix
m <- nrow(geneExpression)
n <- ncol(geneExpression)
randomData <- matrix(rnorm(n*m), nrow = m, ncol = n)

nullpvals <- apply(randomData, MARGIN = 1, mytest)
# count # of pvalues <= 0.05
sum(nullpvals<=0.05)
# we see that we get 400+ p-values which are significant even though there should be none#
# the problem is that p-values are random variables (RV) and by chance, we expect to see 5% of them to be lower than 0.05
# -> this multiple testing using 2 sample t-test will generate many false positives (FP) even when the null hypothesis is true for all features