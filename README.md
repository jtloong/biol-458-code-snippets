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
