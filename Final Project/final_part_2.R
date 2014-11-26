# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
# Final Project Part #2
# 
# author: Brett Berry
# 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library("BRugs")

setwd("~/Bayes/Final Project")
data = read.table("seismic.txt", header=TRUE)

#data["d2"] = c(data["d"]^2)
#data["m2"] = c(data["m"]^2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Part a: 
#
# Fit a linear model of the response y in terms of first and second order predictors
# d, m, s, r, d2, m2.  Then use variable elimination to get down to three of these.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Linear model with all 5 parameters, manual elimination

model1 = lm(formula = y ~ m + d + s + r + I(d^2) + I(m^2), data = data)
summary(model1)

# Eliminate s
model2 = lm(formula = y ~ m + d + r + I(d^2) + I(m^2), data = data)
summary(model2)

# Eliminate d
model3 = lm(formula = y ~ m + r + I(d^2) + I(m^2), data = data)
summary(model3)

# Eliminate m
model4 = lm(formula = y ~ r + I(d^2) + I(m^2), data = data)
summary(model4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Part b: 
#
# Use the step command for step-wise variable selection to confirm the model from
# variable elimination. (Set scale to log(n) for BIC or 2 for AIC, and if there are
# differences among the three methods go with majority rule in the end.)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Variable elimination to 3 predictors using forward step function
fwdModel = step(lm(y ~ 1, data = data),
                direction="forward", 
                scope =(~m + d + s + r + I(d^2) + I(m^2)), steps=3)

#BIC Method for variable elimination
modelBIC = step(model1, k = log(100))

#AIC Method for variable elimination
modelAIC = step(model1, k = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Part c: 
#
# Look at the residuals and choose one of log, √ or no transformation on y to
# make model assumptions reasonable
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Original residuals without transformations

boxplot(modelAIC$residuals)
qqnorm(modelAIC$residuals)

# Square root transformation and plots

sqrtmodel = lm(formula = sqrt(y) ~ r + I(m^2) + I(d^2), data = data)

boxplot(sqrtmodel$residuals)
qqnorm(sqrtmodel$residuals)

# Log transformation and plots

logmodel =  lm(formula = log(y) ~ r + I(m^2) + I(d^2), data = data)

boxplot(logmodel$residuals)
qqnorm(logmodel$residuals)

# The square root transformation produces the best plot results.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Part d: 
#
# Use the linear model to predict y at Portland if a magnitude 8 earthquake hits
# the Cascadia fault along the coast (say 100 km away). At Portland, the rock is
# a combination of r = 0 (newer rock formation) and r = 1 (tertiary rock). Do the
# prediction at both values of r and include a prediction interval. Also the soil s is
# alluvial so s = 1 if that is in your model. (Note: if you transformed y say with
# √y for better diagnostics, be sure to untransform the prediction interval.)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Rock type 0, newer rock formation
Portland = predict(modelAIC,data.frame(m = 8.0, d = 100, r = 0))

# Rock type 1, tertiary rock
Portland1 = predict(modelAIC,data.frame(m = 8.0, d = 100, r = 1))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Part e: 
#
# Use BUGS on three predictors to get parameter estimates on β using a g-prior
# with g = n = 100, σ^2 = s^2, and centered at the least-squares estimates, and use
# a Gamma(a=.25, b=.25) prior on 1/σ2. (Note: BUGS uses the precision matrix (XtX)/(s
# 2n) in dmnorm.)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write the BUGS models into the file "model.txt"
cat("model
{
    beta0 ~ dnorm(0,100)
    beta1 ~ dnorm(0,100)
    beta2 ~ dnorm(0,100)
    
    for(i in 1:N)
{
    Derogatory.reports[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1*Income.per.dependent[i] + beta2*Monthly.credit.card.exp[i]
} 
    
}", file="model.txt")

BRugsFit(
  modelFile = "model.txt",
  data = list(
        "Income.per.dependent"=data$Income.per.dependent,
        "Monthly.credit.card.exp"=data$Monthly.credit.card.exp,
        "Derogatory.reports"=data$Derogatory.reports,
        "N"=100
        ), 
  inits =  list(beta0=rnorm(1,0,1),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1)),
  numChains = 3,
  parametersToSave = c("beta0","beta1","beta2"),
  nBurnin = 1000, 
  nIter = 10000,
  nThin = 10
)