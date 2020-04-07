# week 3 lectures notes: Poisson Example for RNA-seq

### simulated gene expression data

N=10000 #number of genes
lambdas=2^seq(1,16,len=N)

#create two independent poisson data replicates
y=rpois(n = N, lambda = lambdas)
x=rpois(n = N, lambda = lambdas)
#find indices where both x and y are positive, because we want to use log-ratio as summary statistic to summarize the differences 
ind=which(y>0 & x > 0)

library(rafalib)
mypar()
splot(log2(lambdas), log2(y/x), subset = ind)
# observation from plot: when lambda small, difference between x and y is much higher than at the high end of the lamda spectrum

### ### real gene expression data (now with real RNA-seq data)

library(parathyroidSE) # if missing run: "BiocManager::install("parathyroidSE")"
data(parathyroidGenesSE)
se <- parathyroidGenesSE

# take two random replicates from the data
x <- assay(se)[,23]
y <- assay(se)[,24]
# again only consider cases where values are positive for both replicates
ind=which(y>0 & x>0) # make sure no 0s due to ratio and log
# use ovaerage of log2(x) and log2(y) as an estimate of lambda
splot( (log2(y)+log2(x))/2 , log2(y/x), subset = ind)
# observation: very similar behaviour with the actual data i.e. plots similar trend in both case: simulation and real RNA data