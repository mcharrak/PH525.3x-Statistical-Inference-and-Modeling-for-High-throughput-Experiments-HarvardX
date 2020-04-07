# week 4 lectures notes:

# exploratory data analysis techniques
# plots; Volcano plots and p-value histograms, boxplots, and MAplots

library(genefilter)
library(GSE5859Subset)
data(GSE5859Subset)

# groups (cases and controls)
g <- factor(sampleInfo$group)
# perform t-test for each gene expression
results <- rowttests(x = geneExpression, fac = g)

# create another dataset where we know that for all genes/rows/features, the null hypothesis (H0) is true
m <- nrow(geneExpression)
n <- ncol(geneExpression)
randomData <- matrix(rnorm(n*m), nrow = m, ncol = n)
nullresults <- rowttests(x = randomData, fac = g)

### p-value histograms

library(rafalib)
mypar(1,2)

hist(nullresults$p.value, ylim = c(0,1400), main = "")
axis(side = 1, at = seq(0,1,0.1))
# observation nullresults: as expected, the p-values look uniformly distributed
hist(results$p.value, ylim = c(0,1400), main = "")
axis(side = 1, at = seq(0,1,0.1))

### volcano plots
mypar(1,1)
plot(results$dm, -log10(results$p.value), xlab = "estimated effect size", ylab = "-log(base 10) p-values")
# we use -log10() because we want significant genes (i.e. small associated p-value) to be high on the volcano plot

### data boxplots and histograms

library(Biobase)
#devtools::install_github("genomicsclass/GSE5859") #install GSE5859 library if you have not done yet
library(GSE5859)
data(GSE5859)
ge <- exprs(e) #ge for gene expressions
#dim(ge) #the data has 208 samples for 8793 gene expressions

# 1) lets plot boxplot of data distribution for all 208 samples
library(rafalib)
mypar(1,1)
boxplot(x = ge, 
        range = 0, 
        names = 1:ncol(e), 
        col = ifelse(1:ncol(ge) == 49, yes = "green", no = "red"),
        main = "Data with no error") # future outlier nr.49 will be color code 1 (i.e. green)


# now, let us manipulate the data from sample number 49 by introducing an error
ge[,49] <- ge[,49]/log2(exp(1)) ## error:data curators use log(base 2) [i.e. log2()] instead of log(base e) [i.e. log()] to transform data

library(rafalib)
boxplot(x = ge, 
        range = 0, 
        names = 1:ncol(e), 
        col = ifelse(1:ncol(ge) == 49, yes = "green", no = "red"),
        main = "Data with error in sample #49")

# to better visualize this outlier we can plot the quantiles instead
# compute the 5%, 25%, 50%, 75%, 95% percentiles
qs <- t(apply(X = ge, MARGIN = 2, FUN = quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
# now plot the graphs for each quantile (i.e. each column of qs)
mypar(1,1)
matplot(x = qs, type = "l", lty = 1, xlab = "sample IDs", yaxt = "n", main = "quantiles (qs) vs. samples (sample IDs)")
at_ = as.numeric(apply(X = qs, MARGIN = 2, FUN = median))
labels_ = c("5%", "25%" , "50%" , "75%" , "95%")
axis(side = 2, at = at_, labels = labels_, las = 1) # argument las to rotate labels of tick marks)
# NOTE: we refer to this plot as a "KABOXPLOT" because Karl Broman was the first we saw to use it as an alternative to boxplots

# alternative: plot all "smooth" histograms for each sample ID (total:208)
mypar(1,1)
shist(ge, unit = 0.5) #smooth histograms (shist) (i.e. empirical density)

### MA plot

x <- ge[ ,1]
y <- ge[ ,2]
mypar(1,2)
plot(x,y)
diff_btw_x_y <- x-y
mean_btw_x_y <- (x+y)/2
plot((x+y)/2, x-y)