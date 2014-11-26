library("MASS")
library("BRugs")
library("sandwich")

setwd("~/Bayes/Final Project")
data = read.csv("creditscore.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part a: Frequentist Poisson Regression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#General linearized model of derogatory reports in terms of age, income, monthly credit card expenditure, home ownersip, and self-employment status

fullReg = glm(Derogatory.reports ~ Age + Income.per.dependent + Monthly.credit.card.exp + Own.home + Self.employed, family = "poisson", data = data)

summary(fullReg)

# Linear model reduced to two predictors 

reducedReg = glm(Derogatory.reports ~ 1, family = "poisson", data = data)

fwdModel = step(reducedReg, 
                direction="forward", 
                scope =(~Age + Income.per.dependent + Monthly.credit.card.exp + Own.home + Self.employed), steps=2)

# Alternative method using BIC 

step(fullReg , k = log(100))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part b: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Wald Frequentist confidence intervals for the two predictor parameter estimates, 
# using "confint" function from the MASS package for simplification.

confint(fwdModel, level = 0.95)

# Alternative method: Using Wald Frequentist (glm) formula: theta_hat +/- 1.96*se(theta_hat)

partialReg = glm(Derogatory.reports ~ Income.per.dependent + Monthly.credit.card.exp, family = "poisson", data = data)

cov.model <- vcovHC(partialReg, type="HC0") # numerically approximate covariance matrix
std.err <- sqrt(diag(cov.model)) # compute standard error
confidenceInt <- cbind( LL = coef(partialReg) - 1.96 * std.err, UL = coef(partialReg) + 1.96 * std.err) #compute confidence intervals

confidenceInt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part c: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# new data parameters
newPerson = data.frame(Income.per.dependent = 5.0, Monthly.credit.card.exp = 1000)

#prediction of new data
predNew = predict(fwdModel, newPerson, type = "response")

sim = replicate(1000, rpois(rep(1, length(predNew)), predNew))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part d: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
  nIter = 100000,
  nThin = 10
)


# rjags way of doing it
# 
# library("rjags")
# 
# inits = function() {
#   list(beta0=rnorm(1,0,1),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1))
# }
# 
# fit = jags.model(
#   file="model.txt",
#   data=list(
#     "Income.per.dependent"=data$Income.per.dependent,
#     "Monthly.credit.card.exp"=data$Monthly.credit.card.exp,
#     "Derogatory.reports"=data$Derogatory.reports,
#     "N"=100
#     ), 
#   inits=inits,
#   n.chains=3,
#   n.adapt=1000)
# 
# update(fit, 1000)
# 
# samples = coda.samples(fit,5000,variable.names=c("beta0","beta1","beta2"),10)

 


