# week 4: Hierarchical Models in Practice Exercises

# install SpikeInSubset data from Bioconductor if you haven't yet
# BiocManager::install("SpikeInSubset")
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95) #' data of 6 replicate samples
# RNA from 16 genes were artificially added in different quantities to each of the six samples

cat("The following gene names were manipulated:\n")
print( colnames(pData(rma95)) ) # here, these quantities and gene IDs are stored

g <- factor(rep(0:1,each=3))

# create an index of which rows are associated with the artificially added genes
spike <- rownames(y) %in% colnames(pData(rma95)) # spike is boolean vector

## exercise 1
library(genefilter)

ttests <- rowttests(x = y, fac = g)
pvals <- ttests$p.value
alpha <- 0.01
calls <- (pvals < alpha)

true_positives <- sum (calls & spike)
positives <- sum(calls)
false_positives <- positives - true_positives

prop_ex1 <- false_positives/positives
print(prop_ex1)

# volcano plot: x-axis [difference of the group means] & y-axis [-log10(p value)]

# create mask for pvalue smaller than significance level alpha of 0.01 && small effect-size (i.e. difference in means)
mask <- with(ttests, abs(dm) < 0.2 & p.value < 0.01)
cols <- ifelse(mask, "red", ifelse(spike, "dodgerblue", "black"))

# red: abs(dm) < 0.2 & pvalue < 0.01
# dodgerblue: if gene is truly significant
# black: if gene is not significant

with(ttests, 
     plot(-dm, -log10(p.value), cex = 0.8, pch = 16, 
          xlim = c(-1,1), ylim = c(0,5), xlab = "Difference in means / effect size",
          col = cols, main = "Volcano plot")
     )
abline(h = 2, v = c(-0.2, 0.2))
# vertical line indicates border between genes with "large" difference in means and genes with "small" difference in means
# horizontal line indicates genes with very small p-value (b/c small pvalue => large log(pvalue))

## exercise 2

# get indices for true positives, false positives, true negatives, and false negatives
TPs <- (calls==T) & (spike==T)
FPs <- (calls==T) & (spike==F)
TNs <- (calls==F) & (spike==F)
FNs <- (calls==F) & (spike==T)

# use only data from group 1
y_group2 <- y[,1:3]

# take subset of data from group 1 for each of the 4 cases (TP,FP,TN,FN)
y_TP <- y_group2[TPs,]
sd_TP <- apply(X = y_TP, MARGIN = 1, FUN = sd)

y_FP <- y_group2[FPs,]
sd_FP <- apply(X = y_FP, MARGIN = 1, FUN = sd)

y_FN <- y_group2[FNs,]
sd_FN <- apply(X = y_FN, MARGIN = 1, FUN = sd)

y_TN <- y_group2[TNs,]
sd_TN <- apply(X = y_TN, MARGIN = 1, FUN = sd)

mypar(1,4)
LIM <- range(c(sd_FN, sd_FP, sd_TP, sd_TN))
# boxplot TP
boxplot(x = sd_TP, main = "TP within group SD", xlab = "TP", ylab = "gene expression", ylim = LIM)
# boxplot FP
boxplot(x = sd_FP, main = "FP within group SD", xlab = "FP", ylab = "gene expression", ylim = LIM)
# boxplot TN
boxplot(x = sd_TN, main = "TN within group SD", xlab = "TN", ylab = "gene expression", ylim = LIM)
# boxplot FN
boxplot(x = sd_FN, main = "FN within group SD", xlab = "FN", ylab = "gene expression", ylim = LIM)

#### alternative solution
mypar(1,1)
sds <- rowSds(y[,g==0])
index <- paste0( as.numeric(spike), as.numeric(ttests$p.value<0.01))
index <- factor(index,levels=c("11","01","00","10"),labels=c("TP","FP","TN","FN"))
boxplot(split(sds,index))

# OBSERVATION FROM EXERCISES 1&2: 
# random variability associated with the sample standard deviation (denominator of test statistic) leads to t-statistics 
# that are large by chance. Note that the sample standard deviation we use in the t-test is an estimate and
# that with just 3 samples per gene, the variability associated with the denominator in the t-test can be large.

## exercise 3

# install limma package from Bioconductor if you haven't yet
# BiocManager::install("limma")
library(limma)

fit <- lmFit(object = y_group2)
colnames(coef(fit))
fit <- eBayes(fit)

# extract sample SD
sampleSD <- fit$sigma
posteriorSD <- sqrt(fit$s2.post)

# first plot sampleSD in red
y1 <- rep(x = 1, times = length(sampleSD))
y2 <- rep(x = 2, times = length(posteriorSD))
plot(sampleSD, y1, ylim = c(1, 2), pch = 19, yaxt = "n",
     col = "red", ylab = "Before (red) vs. after (blue)", xlab = "SD estimates", main = "Standard deviation estimates")
# add points from posteriorSD
points(posteriorSD, y2,
       col = "dodgerblue", pch = 19)
axis(side = 2, at = 1:2, labels = c("before","after"), las = 1) # argument las to rotate labels of tick marks

cat("Mean of SD from samples:\n")
cat(mean(sampleSD))

cat("Mean of SD using 'limma' analysis:\n")
cat(mean(posteriorSD))

# answer & observation: the small SD values are pushed toward higher SD values and larger SD values are pushed toward smaller SD values
# therfore final answer: Moves all the estimates of SD closer to 0.12 (here move from 0.11917 to 0.11954)

## exercise 4

# now we use the new sample standard deviation estimates for our t-statistics in order to re-calcualte p values
##second coefficient relates to diffences between group
fit = lmFit(y, design=model.matrix(~ g))
fit = eBayes(fit)
##second coefficient relates to diffences between group
corrected_pvals = fit$p.value[,2]

alpha <- 0.01
corrected_calls <- (corrected_pvals < alpha)

corrected_true_positives <- sum (corrected_calls & spike)
corrected_positives <- sum(corrected_calls)
corrected_false_positives <- corrected_positives - corrected_true_positives

prop_ex4 <- corrected_false_positives/corrected_positives
print(prop_ex4)

# corrected volcano plot: x-axis [difference of the group means] & y-axis [-log10(p value)]

# create mask for pvalue smaller than significance level alpha of 0.01 && small effect-size (i.e. difference in means)
corrected_dms <- fit$coefficients[,2]
corrected_mask <- abs(corrected_dms) < 0.2 & corrected_pvals < 0.01
corrected_cols <- ifelse(corrected_mask, "red", ifelse(spike, "dodgerblue", "black"))

# red: abs(dm) < 0.2 & pvalue < 0.01
# dodgerblue: if gene is truly significant
# black: if gene is not significant

##second coefficient relates to diffences between group
plot(corrected_dms, -log10(corrected_pvals), cex = 0.8, pch = 16,
     xlim = c(-1,1), ylim = c(0,5), xlab = "Difference in means / effect size", 
     col = cols, main = "Corrected volcano plot")
abline(h = 2, v = c(-0.2, 0.2))
# vertical line indicates border between genes with "large" difference in means and genes with "small" difference in means
# horizontal line indicates genes with very small p-value (b/c small pvalue => large log(pvalue))

# OBSERVATION: compared to previous volcano plot we can see that we no longer have small p-values for genes with small effect sizes

