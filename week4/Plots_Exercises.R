# week 4: plots exercises

# download and install the Bioconductor package SpikeInSubset and then load the library and the mas133 data

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("SpikeInSubset") #install SpikeInSubset if you do not have it yet
library(SpikeInSubset)
data(mas133)


## exercise 1

e <- exprs(mas133)
# plot first two samples
sample1 <- e[,1]
sample2 <- e[,2]
cor_samples <- cor(sample1, sample2)
mypar(1,1)
plot(sample1, sample2, main = paste("corr =", signif(cor_samples, digits = 3)), cex = 0.5)

k <- 3000
b <- 1000 #a buffer
polygon(c(-b, k, k, -b), c(-b, -b, k, k), border = "red")

# we see that all values are strictly positive
max_val <- k
count_points_inside <- sum(sample1 <= max_val & sample2 <= max_val)
count_points_total <- length(sample1)

prop_ex1 <- count_points_inside/count_points_total
print(prop_ex1)

## exercise 2

log_sample1 <- log2(sample1)
log_sample2 <- log2(sample2)
cor_log_samples <- cor(log_sample1, log_sample2)
mypar(1,1)
plot(log_sample1, log_sample2, main = paste("log-values & corr =", signif(cor_log_samples, digits = 3)), cex = 0.5)

k_new <- log2(k)
b_new <- log2(0.5) #log2(1/2) = log2(1) - log2(2) = 0 - 1 = -1
polygon(c(-b_new, k_new, k_new, -b_new), c(-b_new, -b_new, k_new, k_new), border = "red")
# OBSERVATION:
# compared to the previous plot, when we now take the log2() transformation of our data, 
# then 95% of our data (range from ~ 5000 to ~30000) is no longer in a tiny section of plot

# answer: The tails (here actually the left tail) do not donimate the plot

## exercise 3

# convert expression data into log2 scale
e <- log2(exprs(mas133))

# create and MA plot (is equiv. to Bland-Altman plot) for two samples s1 and s2
# MA plot x-axis: average of 2 samples: (s1+s2)/2
# MA plot y-axis: difference btw 2 samples: s1-s2
s1 <- e[,1]
s2 <- e[,2]
plot(x = (s1+s2)/2, y = s2-s1, cex = 0.5, main = "MA plot: Difference vs. average (aka Bland-Altman plot)")


# because the values are in log scale, we need to consider that the value we are looking for, here log ratio: log2(s2/s1)
# equals to log2(s2) -log2(s1). Now because we already applied the log2() transformation in line 54, we can simply compute the log ratio as:
# s2 - s1. Thus, we need to compute the SD (standard deviation) of log-ratio= s2 - s1
log_ratio <- s2 - s1
SD_log_ratio <- sd(log_ratio)
print(SD_log_ratio)
cat("alternative answer:\n")
print(sqrt(mean( (e[,2]-e[,1])^2)))

## exercise 4

diffs <- s2 - s1
abs_diffs <- abs(diffs)
# question: how many cases do we have, where the difference is larger than a factor of 2. 
# We know that for log2() transformed values, a increase of factor 2 in the original data corresponds to increase of 1 unit in the log2() transformed values
counts_ex4 <- sum(abs_diffs > 1)
print(counts_ex4)