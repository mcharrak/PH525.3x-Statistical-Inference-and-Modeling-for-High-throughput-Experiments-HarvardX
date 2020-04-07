# week 4: Bayes Rule in Practice Exercises

library(dplyr)

# First download some baseball statistics.
tmpfile <- tempfile()
tmpdir <- tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
##this shows us files
filenames <- unzip(tmpfile,list=TRUE)
players <- read.csv(unzip(tmpfile,files="Batting.csv",exdir=tmpdir),as.is=TRUE)
unlink(tmpdir)
file.remove(tmpfile)

## exercise 1

# select batting averages for players with more than 500 at bats (500) in year 2012
# here we have: H... hits, AB... at bats (i.e. hit trials)
players_ss_ex1 <- filter(.data = players, yearID == 2012) %>% mutate(AVG=H/AB) %>% filter(AB>=500)
AVGs_2012 <- players_ss_ex1 %>% select(AVG)

## exercise 2

players_ss_ex2 <- filter(.data = players, yearID %in% c(2010,2011,2012)) %>% filter(AB >= 500)  %>% mutate(AVG = H/AB)
AVGs_ex2 <- players_ss_ex2 %>% select(AVG) %>% pull()
AVG_ex2 <- mean(AVGs_ex2)
print(AVG_ex2)

# alternative solution from class:
#dat <- filter(players,yearID>=2010, yearID <=2012) %>% mutate(AVG=H/AB) %>% filter(AB>500) %>% select(AVG)
#print(mean(dat$AVG))

## exercise 3

sd_ex3 <- sd(AVGs_ex2)
print(sd_ex3)

## exercise 4

mypar(1,2)

# normal dist qqplot
qqnorm(y = AVGs_ex2, main = "Normal QQ plot")
qqline(y = AVGs_ex2)

# histogram
hist(x = AVGs_ex2, main = "Histogram of AVG across players")
# answer and observation: Normal

## exercise 5

p <- 0.450
q <- 1 - p
N <- 20 # sample size
n <- 20 # trials

# var(1/N sum_i=1^N[X_i)] = 1/(N^2) * sum_i=1^N[var(X_i)] = N*Var(X_i)/(N^2) = Var(X_i)/N = n*p*q/N = var(...)
# we know that for a binomial dist. we have variance of: n*p*(1-p)
# now we can consider each X_i as a binomial dist with n (#trials) n=1
# therfore we get var(...) = 1*p*q/N

# for further explanation see stackexchange explanation:
# https://stats.stackexchange.com/questions/29641/standard-error-for-the-mean-of-a-sample-of-binomial-random-variables

var_sample_mean <- 1*p*q/N
sd_sample_mean <- sqrt(var_sample_mean)
print(sd_sample_mean)

## exercise 6

theta_prior <- 0.45 # Jose Iglesias batting score average at 20 bats
Y <- theta_prior
sigma <- 0.11 # Jose Iglesias battings score standard deviation at 20 bats

mu <- 0.275 # across players batting score average
tau <- 0.027 # across players battings score standard deviation
# note that tau >> sd

B <- sigma^2 / (sigma^2 + tau^2)

# Y... observed average - here Y = 0.45
# E[theta|Y] = B*mu + (1-B)*Y is the posterior mean prediction; posterior because after we see the data Y from Jose Iglesias hits
theta_posterior <- B*mu + (1-B)*Y
print(theta_posterior)