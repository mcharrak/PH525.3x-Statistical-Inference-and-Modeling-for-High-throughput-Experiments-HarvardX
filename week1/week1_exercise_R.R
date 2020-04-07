# week1 R exercises

library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset) ##this loads the three tables

#ex1
sum(sampleInfo$date == "2005-06-27")
#ex2
bools_chrY <-geneAnnotation$CHR == "chrY"
bools_chrY <- na.omit(bools_chrY) #drop NAs
sum(bools_chrY)
sum(geneAnnotation$CHR=="chrY",na.rm = TRUE)
#ex3
expression_col_name <- sampleInfo$filename[sampleInfo$date=="2005-06-10"]
geneAnnotation_row_idx <- which(geneAnnotation$SYMBOL=="ARPC1A")
geneAnnotation_row <- geneAnnotation[geneAnnotation_row_idx,]
expression_row_name <- geneAnnotation_row$PROBEID
subject_expression <- geneExpression[expression_row_name, expression_col_name]
subject_expression
# alternative solution
i <- which(geneAnnotation$SYMBOL == "ARPC1A")
j <- which(sampleInfo$date == "2005-06-10")
geneExpression[i,j]
#ex4
medians <- apply(geneExpression,2,median)
median <- median(medians)
median
#ex5
function_ex5 <- function(e,group){
  pval <- t.test( e[group==1], e[group==0])$p.value
  return(pval)
}

g <- factor(sampleInfo$group)
pvals <- apply(X = geneExpression, MARGIN = 1, FUN = function_ex5, g)
min_pval <- min(pvals)
min_pval