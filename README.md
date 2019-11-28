# Notes

## Hypothesis Testing

For tests, use z-dist for n>100 and t-dist for less.
Further, pay attention to critical region given one sided vs two 

### Generate a t-distribution

```R
x <- seq(-5,5,0.1)
y <- dt(x, df=10)
plot(x, y, type='l')
```

### Given some samples generate a CI

Remember qt gives left hand side by default, while pt gives right hand side.

```R
samples <- c(25, 31, 30)

sample.mean <- mean(samples)
sample.stdev <- sd(samples)
degrees.of.freedom = length(samples) - 1

tvalues <- qt(c(0.025, 0.975), df=degrees.of.freedom)
ci.upper <- sample.mean + tvalues[2] * sample.stdev/sqrt(degrees.of.freedom + 1)
ci.lower <- sample.mean + tvalues[1] * sample.stdev/sqrt(degrees.of.freedom + 1) 
```

### Using your own t-test

This is assuming equal variance. Unequal variance formula on pg 15 of Lecture 3. 
```R
t = (mean1 - mean2) / (pooled_stdev * sqrt((1/n1) + (1/n2)))
```

### Using t-test function

Remember, by default the t.test function uses a Welsh correction so variances are NOT assumed equal

```R
# Two sample, two sided
t.test(samples, other.samples)

# One sample, two sided)
t.test(samples, mu=mu0)

# One sample, one sided
t.test(sample, mu=mu0, alternative="greater")
```

### Test variance

```R
var.test(sample, other.samples)
```

## Regression

### Linear models

```R
model <- lm(log(lp) ~ X, data = dataframe)
summary(model)
slope <- lm.r$coefficients["X"]
```

### Assess linear model

This shows four diagnostic plots. Look to lab4 at the end for guidance on what they mean. Further consult Dyer page 115 (or 97 in the doc).

```R
layout(matrix(1:4, 2, 2))
plot(model)
```

## Growth Models

### Discrete growth rates between generations

```R
Nt1 <-data$Population[-1]
Nt<-data$Population[1:length(data$Population)-1]
pgrow<-Nt1/Nt
```

### Arithmetic mean model

```R
arith_mean <- mean(pgrow)

time <- 1:length(data$Population)
N0 <- data$Population[1]
badger_model.arith <-sapply(arith_mean, function(lb) N0*lb^time)

plot(x=time, y=Nt.s, log="y")
lines(x=data$Year, y=badger_model.arith)
```

### Geo mean model

```R
r <- mean(log(pgrow))

geo_lambda <- exp(r)
time <- 1:length(data$Population)
N0 <- data$Population[1]
badger_model.geo <-sapply(geo_lambda, function(lb) N0*lb^time)

plot(data, log="y")
lines(x=data$Year, y=badger_model.geo)
```

### Logistic model

Based on building a linear model on growth rates and density

```R
nt <- data$pop.density[1:length(data$pop.density)-1]
nt1 <- data$pop.density[-1]

growth_rates <- nt1 / nt

growth <- data.frame(nt, growth_rates)
model <- lm("growth_rates ~ nt", growth)

yintercept <- model$coefficients[1] # r value
slope <- model$coefficients[2] # 
xintercept <- -yintercept / slope # K-value

dlogistic<-function(K, lambda, N0=2, t=15){ 
  N<-c(N0, numeric(t))
  for (i in 1:t) {
    N[i+1]<- N[i] +  lambda * (1 - (N[i] / K)) * N[i]
  } 
  return(N) 
}

t <- 20
Nts <- dlogistic(K=xintercept, lambda=yintercept, t=t)
```
### Stochastic Growth Models

```R
threshold=24         
project=30
runs=1000  

stoch.pop=matrix(NA,project,runs)
stoch.pop.dens=matrix(NA,project,runs)
stoch.pop[1,]=N0

# two nested loops to create stochastic population sizes

for (i in 1:runs){			
	# looping over 1000 runs of the stochastic model
  for (t in 2:project){		
  	# and looping over 30 years of projection within each of 1000 runs
    lambda=exp(rnorm(1,geomean,sqrt(var.r)))
    # draw a value of lambda from a lognormal distribution
    
    #lambda=exp(sample(x=r,size=1,replace=T)) 
    # or sample from within the existing observed growth rates
    
    stoch.pop[t,i]=stoch.pop[(t-1),i]*lambda
         
    if(stoch.pop[t,i]<=threshold) break  
    # leave the loop if pop <= threshold
  }
}
```

### Stochastic Mean & CI

```R
pop.sd=apply(stoch.pop, 1, function(x) sd(x, na.rm=TRUE))
stoch.pop.mean=apply(stoch.pop, 1, function(x) mean(x, na.rm=TRUE))
ucl =((stoch.pop.dens.mean)+1.96*pop.sd/sqrt(nrow(stoch.pop)))
```

## Matrix Models

### Projecting Population

```R
N.proj<-matrix(0,nrow=nrow(L), ncol=years+1)
N.proj[,1] <- N0
years = 6

for(i in 1:years){ 
  N.proj[,i+1]<-L%*%N.proj[,i]
}
```

### Dominant Eigenvalue and Left Eigenvector
This vector represents the stage proportions at the stable age distribution

```R
eigen.vals <- eigen(L)
q <- which.max((eigen.vals$values))
dom.eig=eigen.vals$values[q]

sad1<-(eigen.vals$vectors[,q])
sad2<-sad1/sum(sad1)
```

### Right Eigenvector
This vector represents the reproductive value of each stage.
```R
tel<-eigen(t(L))
repv<-Re(tel$vectors[,which.max(tel$values)])
repv<-repv/repv[1]
```

### Calculating Sensitivities and Elasticities

```R
S.dem<-sum(repv*sad2)
s=L
n=dim(L)[1]
for(i in 1:n) {
  for(j in 1:n) {
    s[i,j]=repv[i]*sad2[j]/S.dem;
  }
}

elas<-Re(L/dom.eig*s)
```

## Plotting

### Bar Plots With Segments
Here you pass the plot object (b) and give it the upper and lower bounds of the intervals

```R 
b <- barplot(c(sample.mean, other.sample.mean), names=c("Our Data", "Other Team's Data"), ylab="Yellow Beans per 30 ml", xlab='Groups')
segments(b, c(myu, thu), b, c(myl, thl), lwd=1.5)
```

### Scatter Plot Options

Where pch is the type of point and col is color
```R
plot(x,y, ylab="Endangered species recognized by COSEWIC", xlab="Population of Canada", pch=20, col="red")
```

### Plotting regression

```R
abline(model)
```

### Plotting least squared lines

This code snippet plots distances between your model and the observed values. It takes coordinates where the first two terms are the vectors of x1, y1 with y1 being the predicted point and the second two terms are x2, y2 where y2 is the osberved point. (Obviously they have the same x value).

```R
segments(popdata$X, fitted(model), popdata$X, log(popdata$y))
```
### Matrices

```R
matplot(x=(year1+1):(year1+project), y=stoch.pop, log="y", type="l", ylab='Ln Population Size', xlab=paste('Years'), col=cm.colors(ncol(stoch.pop)))
x=(year1):(year1+project)
lines(x,y=N[20]*exp(geomean)^(1:length(x)), col='blue',type="l")

lines(x[-1],ucl,'l', col  = 'red', lwd=2, lty=2)
abline(h=threshold, col='red', lwd=2, lty=2)

```
### Extinction Probabilities
```R
time=apply(stoch.pop,2, function(x) max(which(x>0)))
#plot histogram of extinction times
nbin<-hist(time,breaks=seq(0,project,by=1), main="distribution of extinction times")

#plot cumulative probabiliy of extinction
td<-cumsum(nbin$counts/runs)
plot(1:project-1,td, xlab="years", ylab="cumulative probability of extinction", "l")
mdext<-median(time)
abline(v=mdext,lwd=2, col="red")
mlabel<-paste("median time to quasi-extinction = ", mdext, sep=" ")
text(x=14, y=.2,mlabel)
```


## Sampling

Rememeber to keep everything in the same units.

### Calculating mean and variation of mark recapture sampling protocol

```R
mean <- (m1 + 1) * (n2 + 1) / (m2 + 1)
var <- (mr_xbar^2) * (n2 - m2 + 1) / ((n2 + 1) * (m2 + 2))
stdev <- sqrt(mr_var)
df <- n2 - 1
```

## General R tips

### Apply functions

Useful link: http://www.datasciencemadesimple.com/apply-function-r/

tapply is by index of one column apply function in another

sapply applies an anonymous function against a vector returning a vector of the same length

### Create a dataframe from vectors

```R
growth <- data.frame(nt, growth_rates)
```

### Importing functions

```R
source("dlogistic.R")
```

### Matrix Notation

When doing matrix selection it's [row, column]:
* x[1, 2] selects the element in the first row at the second column
* x[1, ] selects all the data in the first row
* x[, 1] selects all the data in the first column

### Matrix Multiplication

```R
L%*%N0
```
