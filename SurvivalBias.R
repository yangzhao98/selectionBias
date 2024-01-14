## ---
## Eliminating Survival Bias in Two-stage Instrumental Variable Estimators
## Ref. Vansteelandt et al. (2018). Epidemiology; 29; 536-541.
## ---
## install.packages("timereg")
library(timereg)

## Generate artificial test dataset
n <- 500                             ## Total sample size
u <- rnorm(n)                        ## Unconfounders 
z <- rbinom(n,1,0.5)                 ## The IVs
a <- rnorm(n,3+z-u)                  ## The exposure
ti<- rexp(n, 1/abs(0.5*a + 0.25*u))  ## The observed survival (or censoring) times
d <- ifelse(ti < 10, 1, 0)           ## The survival indicator
t0 <- runif(n, 0, 0.5)               ## The observed entry times
id <- 1:length(a)                    ## The IDs

## Vector of observed event times
times <- sort(unique(c(t0,ti)[t0<ti]))
n <- length(times)

## Construct longitudinal dataset which expresses all observed at risk periods
dataset <- data.frame(
    #y = rep(y, each = n),
    ti = rep(ti, each = n),
    t0 = rep(t0, each = n),
    d = rep(d, each = n),
    a = rep(a, each = n),
    z = rep(z, each = n),
    id = rep(id, each = n),
    start = rep(c(0, times[-n]), length = n),
    stop = rep(times, length = n)
)

dataset <- dataset[dataset$t0 < dataset$ti,]
dataset <- dataset[dataset$t0 <= dataset$start,]
dataset <- dataset[dataset$ti>=dataset$stop,]
dataset$d[dataset$ti > dataset$stop] <- 0


## Construct predictions M(t) - denoted m

dataset$m <- NULL
for (i in 1:(length(times) - 1)) {
    ty <- times[i]
    s <- (dataset$ti >= ty) & (dataset$t0 <= ty)
    dataset$m[dataset$start == ty] <- predict(
        lm(a~z, data = dataset[s,]),
        newdata = data.frame(z = dataset$z[dataset$start == ty]))
}

aalen(Surv(start,stop,d)~const(m),data = dataset)


res <- matrix(0,ncol=3,nrow=100)
delta <- resid(lm(a~z))
for (i in 1:100){
    # choose a range of values theta up to the maximum accumulated exposure effect
    # here an exposure effect of 0.08 which accumulates for at most 30 years
    res[i,1] <- (30*i/100)*0.08
    # calculate the correlation
    res[i,2] <- cor(z,exp(-res[i,1]*delta))
    # test for evidence of a correlation
    res[i,3] <- cor.test(z,exp(-res[i,1]*delta))$p.value
}
par(mfrow=c(1,2))
plot(res[,1],res[,2],xlab=expression(theta),ylab="Correlation")
plot(res[,1],log(res[,3]),xlab=expression(theta),ylab="log p-value",
     ylim=c(log(0.001),0))
abline(h=log(0.05))
