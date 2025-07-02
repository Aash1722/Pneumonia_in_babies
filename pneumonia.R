rm(list=ls())

options(width=60,length=200,digits=3)   # Set number of digits to be printed                                                             

library(survival)
library(KMsurv)  # Read data from the R package KMsurv;
library(MASS)

data(pneumon)

pdata <- pneumon[,c("hospital","mthage","alcohol","smoke","bweight","race","agepn")] #Only including factors relevent to our research

head(pdata)

pdata$alcohol[pdata$alcohol<=1] <- 0 # If mother has less than 1 drink per month, assign value 0 (no alcohol use)
pdata$alcohol[pdata$alcohol>1] <- 1 # If mother has 1 drink or more per month, assign value 1 (alcohol use)

pdata$smoke[pdata$smoke>0] <- 1 # If mother has more than 0 packs of cigarettes per day, assign value 1 (cigarette use)

attach(pdata)

## KAPLAN-MEIER ESTIMATES ##
# Non-parametric survival function estimate including confidence bands
fit1 <- survfit(Surv(agepn, hospital)~1, conf.type = "log-log", data = pdata)
# Non-parametric survival curve and confidence bands
plot(fit1,  lwd=2, mark.time=TRUE, cex=2, cex.lab=1.4,xlab="age child in hospital",ylab="probability",
     main="Survival Function Estimates (Kaplan-Meier Plot)", conf.int=T)


## NELSON-AALEN ESTIMATES ##
#Including confidence bands
fitNA1 <- survfit(Surv(agepn, hospital) ~ 1,ctype = 1, conf.type = "plain", data = pdata)
# Plotting
plot(fitNA1, lwd=2, lty=2, conf.int=FALSE, mark.time=TRUE, cex=2, cex.lab=1.4,
     cumhaz = TRUE, col="black",
     xlab="age child in hospital", ylab="Cumulative Hazard",main="The Cumulative Hazard estimates")

# CHECKING FOR NO TIES
table(pdata$agepn[pdata$alcohol==1])
# This survival data has ties since tied events exist,
# i.e. there are instances where two or more events occur at the same time.

## FITTING SEMI-PARAMETRIC COXPH MODEL ##
# Null model
coxph.fit.null <- coxph(Surv(agepn, hospital)~1, data = pdata)
sum.coxph.fit.null <- summary(coxph.fit.null)
print(sum.coxph.fit.null)
# Model with race
coxph.fit.r <- coxph(Surv(agepn, hospital)~factor(race), data = pdata)
sum.coxph.fit.r <- summary(coxph.fit.r)
print(sum.coxph.fit.r)
# Using Breslow's method instead of Efron's which is the default
coxph.fit.r2 <- coxph(Surv(agepn, hospital)~factor(race), method = "breslow", data = pdata)
sum.coxph.fit.r2 <- summary(coxph.fit.r2)
print(sum.coxph.fit.r2)
## Differences between values among methods are due to the fact that there exist tied events
# Model with all predictors
coxph.fit.all <- coxph(Surv(agepn, hospital)~mthage+factor(alcohol)+factor(smoke)+factor(bweight)+factor(race), data = pdata)
sum.coxph.fit.all <- summary(coxph.fit.all)
print(sum.coxph.fit.all)

## VARIABLE SELECTION ##
# Based on AIC through backward selection - recall smaller AIC is preferred

# Main effects only
coxph.fitbackward <- stepAIC(coxph.fit.all, direction="backward",trace=-1)
coxph.fitbackward
print(coxph.fitbackward$anova)
# From this variable selection method, the only fixed-time factors to be included
# in the final model are cigarette use during pregnancy and whether the baby had
# a normal birthweight.

# Main effects and two-way interactions
coxph.fit.allsq <- coxph(Surv(agepn, hospital)~(mthage+factor(alcohol)+factor(smoke)+factor(bweight)+factor(race))**2, data = pdata)
coxph.fitbackward2 <- stepAIC(coxph.fit.allsq, direction="backward",trace=-1)
coxph.fitbackward2
print(coxph.fitbackward2$anova)
# From this,the final model should include the age of the mother at birth,
# cigarette use during pregnancy, whether the baby had a normal birthweight, and
# the interaction between the mother's age at birth and whether the baby had
# a normal birthweight.

cox.final <- coxph(Surv(agepn, hospital)~mthage+factor(smoke)+factor(bweight)+(mthage*factor(bweight)), data = pdata)
sum.cox.final <- summary(cox.final)
print(sum.cox.final)

## CHECKING MODEL ASSUMPTIONS ##

# 1. Non-informative (or random) censoring
# There are opportunities for informative censoring, i.e. having participants
# drop out of the study for reasons related to the study, which could lead
# to biased results.

# 2. Survival times are independent
# We found the presence of ties earlier in the code. So this assumption is violated.

# 3. Hazards are proportional, i.e. hazard ratios are constant over time
# Note: There are no time-varying covariates

#To include time-dependent covariates in the Cox PH model, we use the following function
#expand.breakpoints () to expand the data set first; 

expand.breakpoints <-function(dataset, index = "patnum", status = "status", tevent = "time", 
                              Zvar = F, breakpoints = NULL)
{
  # Expand <dataset> so that each individual has a row for each 
  # unique failure time that is <= his/her/its own failure time.
  #
  # ERROR checking
  onceonly <- paste("must be a character string matching exactly", 
                    "one name in the first parameter, dataset") 
  # Require dataset to be of type data.frame
  if((missing(dataset) || !is.data.frame(dataset)))
    stop("\n Parameter, dataset, must be a data frame")
  varset <- names(dataset)
  covset <- varset
  lvset <- 1:length(varset) # Require dataset to have unique names
  if(any(duplicated(varset))) {
    stop(paste("\n Parameter, dataset, must have uniquely defined", 
               "column names"))
  }
  # Require index to match exactly 1 dataset name
  if(length((indexloc <- lvset[varset == index])) != 1)
    stop(paste("\n Parameter, index,", onceonly))
  covset <- covset[covset != index] 
  # Require status to match exactly 1 dataset name
  if(length((statusloc <- lvset[varset == status])) != 1)
    stop(paste("\n Parameter, status,", onceonly))
  covset <- covset[covset != status] 
  # Require tevent to match exactly 1 dataset name
  if(length((teventloc <- lvset[varset == tevent])) != 1)
    stop(paste("\n Parameter, tevent,", onceonly))
  covset <- covset[covset != tevent] # -----------------------------
  # Form vector of breakpoints, if necessary
  if(is.null(breakpoints)) {
    times <- dataset[, tevent]
    breakpoints <- sort(unique(times[dataset[, status] == 1]))
  }
  # Require breakpoints to be a vector of length >= 1
  if((is.null(breakpoints)) || (!(is.numeric(
    breakpoints))) || (!(is.vector(
      breakpoints))) || (length(breakpoints) < 1)) stop(paste(
        "\n Parameter, breakpoints, must be a numeric vector", 
        "with at least 1 element")) #*****************
  #Begin
  #*****************
  n <- nrow(dataset)
  temp.stop <- c(breakpoints, 0) # The 0 is a place-filler
  if(breakpoints[1] > 0)
    temp.start <- c(0, breakpoints)
  else temp.start <- c(-1, breakpoints)
  n.break <- length(breakpoints) + 1
  t.event <- dataset[, tevent] # ---------------------------
  ## Begin calculate n.epochs
  n.epochs <- sapply(1:n, function(m, breaks, times)
    sum(breaks < times[m]) + 1, breaks = breakpoints, times = t.event) 
  # End n.epochs
  id <- rep(1:n, n.epochs) # Index into supplied dataset
  last.epoch <- cumsum(n.epochs)
  epoch <- unlist(sapply(n.epochs, function(x)
    1:x)) # Index into vectors of interval start & stop points
  if(Zvar) {
    Zmat <- diag(rep(1, n.break))[epoch, ]
    Z.name <- paste("Z", seq(1, n.break), sep = "")
  }
  Tstop <- temp.stop[epoch]
  Tstop[last.epoch] <- dataset[, tevent]
  status2 <- rep(0, length(epoch))
  status2[last.epoch] <- dataset[, status]
  new.dataset <- data.frame(dataset[id, index, drop = F], temp.start[epoch], 
                            Tstop, status2, epoch, dataset[id, covset, drop = F])
  if(Zvar)
    new.dataset <- data.frame(new.dataset, Zmat)
  nam <- c(index, "Tstart", "Tstop", status, "epoch", covset)
  if(Zvar)
    nam <- c(nam, Z.name)
  dimnames(new.dataset) <- list(1:length(epoch), nam)
  return(new.dataset)
}

nobs<-length(agepn)
patnum<-seq(1:nobs) 
new.data<-data.frame(agepn, hospital, patnum, mthage, smoke, bweight)
  # Only including covariates in the final model
data.expand<-expand.breakpoints(new.data, index="patnum", status="hospital", tevent="agepn", Zvar=F)

#Create time-dependent covariate Z2t=Z1*ln(t); Z1=bweight here;

Z2t<-data.expand$bweight*log(data.expand$Tstop)


cox.tdc<-coxph(Surv(Tstart,Tstop, hospital)~mthage+factor(smoke)+factor(bweight)+(mthage*factor(bweight))+Z2t, 
               data=data.expand)
print(cox.tdc)
# Seeing that p-value=0.602 from the Wald test, we have strong evidence that the
# proportional hazards assumption is not met.

# 4. log(Hazard) is a linear function of of the X's/covariates (that are numeric and not categorical/factors)
# Using Martingale residuals

cox.num <- coxph(Surv(agepn, hospital)~mthage, data = pdata)
sum.cox.num <- summary(cox.num)
print(sum.cox.num)

plot(predict(cox.num), residuals(cox.num, type = "martingale"),
     xlab = "fitted values", ylab = "Martingale residuals",
     main = "Residual Plot", las = 1)
# add a line ax y=residual=0
abline(h=0)
# fit a smoother through the points
lines(smooth.spline(predict(cox.num), 
                    residuals(cox.num, type = "martingale")), col="red")

# Looking the residual plot (specifically the red line), the linearity assumption
# appears to be met.

# 5. Values of X (covariates) don't change over time
# Note: There are no time-varying covariates.

# 6. Baseline hazard (h_0(t)) is unspecified
# h_0(t) was not given for this dataset.

## MODEL DIAGNOSTICS ##

# Cox Snell residual plots to access different Cox PH model fitting for covariate bweight

bw0 <- (bweight==0)
bw1 <- (bweight==1)
fit2 <- coxph(Surv(agepn, hospital) ~ factor(bweight), data = pdata)
mres2 = resid(cox.final, type="martingale")
csres2 = hospital-mres2

# Stratified by bweight
r.surv10 <- survfit(Surv(csres2[bw0],hospital[bw0])~1,type="fleming-harrington")
r.surv11 <- survfit(Surv(csres2[bw1],hospital[bw1])~1,type="fleming-harrington")
# Plotting
par(las = 1, mfrow = c(1, 1), mai = c(0.5, 1.0, 1, 0.1), omi = c(1.0, 0, 0.5, 0))
plot(0, 0, lty = 1, type = 'n', xlim = c(0, 0.15), ylim = c(0, 0.12), xlab = "Residual",
     ylab = "Estimated Cum Hazards")
box()
lines(r.surv10$time, -log(r.surv10$surv), type = 's', lty = 1)
lines(r.surv11$time, -log(r.surv11$surv), type = 's', lty = 2)
lines(c(0, 3), c(0, 3), lty = 1)
# Dashed line for babies that had a normal birthweight; 
# solid line for babies not having a normal birthweight
# The estimates being close to the 45-degree line suggests the stratified model 
# fits better than the unstratified model, as a potential way to remedy the violation
# of the proportional hazards assumption.

