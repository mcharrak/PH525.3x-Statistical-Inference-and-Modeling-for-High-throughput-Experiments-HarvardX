# lecture notes week 1

library(devtools)
install_github("genomicsclass/GSE5859Subset") #dataset created for this course

library(GSE5859Subset)
data(GSE5859Subset) ##this loads the three tables

head(geneExpression)
dim(geneExpression) # we have 8793 features (rows) for 24 individuals (columns)

head(sampleInfo)
dim(sampleInfo) # each row represents a column in the geneExpression matrix

identical(colnames(geneExpression), sampleInfo$filename)
## colnames(geneExpression) == sampleInfo$filename

g <- sampleInfo$group
e <- geneExpression[25,]

library(rafalib)
mypar(1,2)
# 1 => case group
qqnorm(e[g==1])
qqline(e[g==1])

# 0 => control group
qqnorm(e[g==0])
qqline(e[g==0])
### we can see that both data is approximately normal distributed

t.test(e[g==1],e[g==0])
### lets repeat this ttest on every feature of the geneexpression dataset
mytest <- function(x){
  t.test(x[g==1],x[g==0],var.equal = TRUE)$p.value
}

mytest(geneExpression[25,])

pvals <- apply(geneExpression, 1, mytest)
length(pvals)