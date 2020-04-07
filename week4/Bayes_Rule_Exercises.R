# week 3: Bayes Rule Exercises

## exercise 1

# Pr(D|+) = [Pr(+|D)*Pr(D)]/Pr(+)
# Pr(+) = Pr(+,D) + Pr(+,no D) = p(D) * p(+|D) + p(no D) * p(+|no D)

p_pos_D <- 0.99
p_neg_D <- 1 - p_pos_D

p_neg_noD <- 0.99
p_pos_noD <- 1 - p_neg_noD

p_D <- 1/4000
p_noD <- 1 - p_D

p_pos <- p_D*p_pos_D + p_noD*p_pos_noD
p_D_pos <- p_pos_D*p_D/p_pos
print(p_D_pos)
# for visual see here https://simplystatistics.org/2014/10/17/bayes-rule-in-a-gif/