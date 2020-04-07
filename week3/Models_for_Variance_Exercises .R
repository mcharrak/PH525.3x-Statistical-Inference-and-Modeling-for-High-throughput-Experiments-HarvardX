# week 3: Models for Variance Exercises

# Install and load the following data library:
  
library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)

# Now load this data and select the columns related to endometrium:

data("tissuesGeneExpression")
library(genefilter)
y = e[,which(tissue=="endometrium")]

## exercise 1

# compute var across all samples for each feature
vars <- rowVars(y)
sds <- sqrt(vars)

mypar(1,2)
# make normal QQ plot
qqnorm(y = vars, main = "Variance")
# plot a line into QQ plot to evaluate whether you see a clear deviation from normality or not
qqline(y = vars)

qqnorm(y = sds, main = "Standard deviation")
# plot a line into QQ plot to evaluate whether you see a clear deviation from normality or not
qqline(y = sds)

# answer: 
# The normal distribution is not a useful approximation here: the left tail is over estimated and the right tail is underestimated.

## exercise 2

# install limma package from Bioconductor if you do not have it yet
#BiocManager::install("limma")
library(limma)

# now use F-dist to fit data
n_df1 <-  14
F_fit <- fitFDist(x = vars, df1 = n_df1)
print(F_fit$scale)

## exercise 3

# theoretical F distribution values
n_df2 <-  F_fit$df2
s0 <-  F_fit$scale
ps <- (seq(along=vars)-0.5) /length(vars)

theo <- qf(p = ps, df1 = n_df1, df2 = n_df2)*s0
obs <- vars

#LIM <- range(c( sqrt(theo), sqrt(obs) ))
LIM <- range(c( sqrt(theo), sqrt(obs) ))
mypar(1,2)
qqplot( sqrt(theo), sqrt(vars), ylim = LIM, xlim = LIM)
abline(0,1)

# remove upper 5% quantiles to get a "zoomed" view for the QQ plot
# get value of the 95% percentile from sqrt(vars) and use this as max value on both axes

axis_max_val <- sqrt( quantile(sqrt(vars), 0.95) )
new_LIM <- c(0,axis_max_val)
qqplot( sqrt(theo), sqrt(vars), ylim = new_LIM, xlim = new_LIM)
abline(0,1)

# observation and answer: If we exclude the genes with the highest variances (top 5%), then the F-distribution provides a good fit.