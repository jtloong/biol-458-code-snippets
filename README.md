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
