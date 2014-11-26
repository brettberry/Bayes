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

data["d2"] = c(data["d"]^2)
data["m2"] = c(data["m"]^2)

X = as.matrix(data[,4:6])
XtX = t(X) %*% X

yhat.ols = (X %*% sqrtmodel$coeff[2:4])


# write the BUGS models into the file "model.txt"
cat("
model
{
	for(i in 1:n){
		m[i] <- inprod(X[i,],beta[]) 
		y[i] ~ dnorm(m[i], precError)
	}

	beta[1:3] ~ dmnorm(mu0[],precMV[,]) #prior on coefficients, mu0=0 
	
	for(k in 1:3){
		for (j in 1:3){
			precMV[k,j]<-XtX[j,k]*precError/n 


			# note: this with precError below are the g-prior with g=n
			# and in terms of precision, giving Sigma0 (Hoff's notation) (prec^-1)* XtX^-1 * g
		}
	}

	precError ~ dgamma(0.25,0.25) #prior on 1/sigma^2 for errors 
	sigma <- pow(precError,-.5) 
}

", file="model.txt" )

library(BRugs)

brugs.fit = BRugsFit(
	modelFile="model.txt",
	data=list(y=data$y,n=100,XtX=XtX,mu0=rep(0,3),X=X),
	inits=list(beta=c(0,0,0),precError=1),
	parametersToSave=c("beta","precError","sigma"),
	numChains=1,
	# nBurnin = 1000,
	# nThin = 10,
	nIter=1000)

brugs.fit$Stats

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Part f: 
#
# Compute the mean-square predictive error of both the Bayesian model and the
# regression model on the test set seismic.test.txt and compare. Compute the
# error on y itself, not the transformed value.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data.test = read.table("seismic.test.txt", header=TRUE)

data["d2"] = c(data["d"]^2)
data["m2"] = c(data["m"]^2)

X = as.matrix(data[,4:6])
XtX = t(X) %*% X

yhat.ols = (X %*% sqrtmodel$coeff[2:4])

yhat.bayes=X%*%(brugs.fit$Stats[1:3,1])
mse.bayes=mean( (data$y-yhat.bayes)^2 ) 
mse.ols=mean((data$y-yhat.ols)^2 ) 
mse.ols
mse.bayes
