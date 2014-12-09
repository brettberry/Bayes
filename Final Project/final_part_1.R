library("MASS")
library("BRugs")
library("sandwich")

setwd("~/Bayes/Final Project")
data = read.csv("creditscore.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# Part A  Do frequentist Poisson regression on the five predictors with glm, then use step
# to get the best model with two predictors
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#
# Part B: For a new person with Income.per.dependent at 5.0 and Monthly.credit.card.exp
# at 1000, estimate the probability that y > 0 by plugging parameter estimates into
# the Poisson probability model (here y denotes the Poisson number of derogatory
# reports).
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Wald Frequentist confidence intervals for the two predictor parameter estimates, 
# using "confint" function from the MASS package for simplification.

confint(fwdModel, level = 0.95)

# Alternative method: Using Wald Frequentist (glm) formula: theta_hat +/- 1.96*se(theta_hat)

partialReg = glm(Derogatory.reports ~ Income.per.dependent + Monthly.credit.card.exp, family = "poisson", data = data)

cov.model <- vcovHC(partialReg, type="HC0") # numerically approximate covariance matrix
std.err <- sqrt(diag(cov.model)) # compute standard error
confidenceInt <- cbind( LL = coef(partialReg) - 1.96 * std.err, UL = coef(partialReg) + 1.96 * std.err) #compute confidence intervals

confidenceInt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# Part c: For a new person with Income.per.dependent at 5.0 and Monthly.credit.card.exp
# at 1000, estimate the probability that y > 0 by plugging parameter estimates into
# the Poisson probability model (here y denotes the Poisson number of derogatory
# reports).
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# new data parameters
newPerson = data.frame(Income.per.dependent = 5.0, Monthly.credit.card.exp = 1000)

#prediction of new data
predNew = predict(fwdModel, newPerson, type = "response")

sim = replicate(1000, rpois(rep(1, length(predNew)), predNew))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#
# Part d:  Do Poisson regression in BUGS with the two predictors, using independent normal
# priors centered at 0 and with variances equal to 100. Report confidence or
# credibility intervals using the posterior distributions and draw the DAG.
# 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# write the BUGS models into the file "model.txt"
cat("model
    {
      beta0 ~ dnorm(0,0.01)
      beta1 ~ dnorm(0,0.01)
      beta2 ~ dnorm(0,0.01)

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
  inits =  list(beta0=0,beta1=0,beta2=0),
  numChains = 3,
  parametersToSave = c("beta0","beta1","beta2"),
  nBurnin = 1000, 
  nIter = 100000,
  nThin = 10
)


# rjags way of doing it
# 
library("rjags")

inits = function() {
  list(beta0=0,beta1=0,beta2=0)
}

fit = jags.model(
  file="model.txt",
  data=list(
    "Income.per.dependent"=data$Income.per.dependent,
    "Monthly.credit.card.exp"=data$Monthly.credit.card.exp,
    "Derogatory.reports"=data$Derogatory.reports,
    "N"=100
    ), 
  inits=inits,
  n.chains=3,
  n.adapt=1000)

# burn in
update(fit, 1000)

samples = coda.samples(fit,100000,variable.names=c("beta0","beta1","beta2"),10)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# Part E: For the new person, get the posterior predictive probability of one or more derogatory
# reports in the Bayesian model P(Y > 0 | y1, . . . , y100) and compare with
# the frequentist number above.
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# write the BUGS models into the file "model.txt"
cat("
model {

	for(i in 1:n) {
		log(m[i]) <- inprod(X[i,],beta[])
		y[i] ~ dpois(m[i])
		ypp[i] ~ dpois(m[i])
	}

	log(meannew) <- beta[1] + 5*beta[2] + 1000*beta[3]
	ynew ~ dpois(meannew)

	beta[1] ~ dnorm(0,0.01)
	beta[1] ~ dnorm(0,0.01)
	beta[1] ~ dnorm(0,0.01)

	pynew <- step(ynew - 1)

	for ( i in 1:100) {
		N[i] <- 1 - step(ypp[i] - 1)
	}

	Nsum <- sum(N[])
	pN <- step(Nsum - 82)

}", file="model.txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Part F:
# 
# This data set shows a common issue with the Poisson model, namely too many
# 0s. Sometimes the data is modeled as a ZIP (zero-inflated Poisson) model to
# deal with this, where the population is viewed as a mixture of people with no
# chance of bad reports and a population at risk. Here we will focus on seeing how
# bad the problem is.
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# i. Do a scatter plot of the fitted values (found using predict.glm on the fitted
# model) against the actual y values and interpret the plot.
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


plot(predict(partialReg, type="response"), data$Derogatory.reports)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ii. The number of 0s in the data is 82. To see how typical this is, simulate 2000
# new data sets with the same predictors using the Poisson distribution at the
# fitted parameter values from glm. Compute the proportion of the data sets
# where the number of 0s is at least 82. (This is a parametric boostrap.)
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


yNew = matrix(0,2000,100)
nZeros = rep(0,2000)
meanY = mean(data$Derogatory.report) 

for (i in 1:2000) {
	yNew[i,] = rpois(100,meanY)
	nZeros[i] = sum(yNew[i,] == 0)
}

nAtleast82 = sum(nZeros >= 82)
proportionAtleast82 = nAtleast82/2000


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
# iii. For the Bayesian version of the parametric bootstrap, use the posterior predictive
# distribution to estimate the probability P(N ≥ 82 | y1, . . . , y100),
# where N is the number of 0s in a new data set Y˜ = (Y1, . . . , Y100). (Note:
# this problem is similar in spirit to Exercise 4.3 in Hoff and posterior predictive
# distributions show up on p. 60.)
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write the BUGS models into the file "model.txt"
cat("
model {

	for(i in 1:n) {
		log(m[i]) <- inprod(X[i,],beta[])
		y[i] ~ dpois(m[i])
		ypp[i] ~ dpois(m[i])
	}

	log(meannew) <- beta[1] + 5*beta[2] + 1000*beta[3]
	ynew ~ dpois(meannew)

	beta[1] ~ dnorm(0,0.01)
	beta[2] ~ dnorm(0,0.01)
	beta[3] ~ dnorm(0,0.01)

	pynew <- step(ynew - 1)

	for ( i in 1:100) {
		N[i] <- 1 - step(ypp[i] - 1)
	}

	Nsum <- sum(N[])
	pN <- step(Nsum - 82)

}", file="model.txt")


library("rjags")

X = as.matrix(data[,2:3])
X = cbind(rep(1,100),X)

fit = jags.model(
  file="model.txt",
  data=list(y=data$Derogatory.reports,X=X,n=100),
  inits=list(beta=rep(0,3)),
  n.chains=1,
  n.adapt=1000)

# burn in
update(fit, 1000)

samples = coda.samples(fit,100000,variable.names=c("Nsum","pN"),10)

# 
a = as.array(samples[[1]])
sum(a[,2])
