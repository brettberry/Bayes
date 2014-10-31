# Midterm Project: October 29, 2014
# Problem 5.1abc
# Brett Berry

# Set given parameters
mu_0 = 5
s2_0 = 4
k_0  = 1
nu_0 = 2

# Part a: Calculate posterior means and 95% CI for each school

# Set directory to locate data files
setwd("~/Downloads")

# Read in data files
school1 = read.delim("school1.dat", sep="\n", header=FALSE)[[1]]
school2 = read.delim("school2.dat", sep="\n", header=FALSE)[[1]]
school3 = read.delim("school3.dat", sep="\n", header=FALSE)[[1]]

# Sample means for each school
ybar1 = mean(school1)
ybar2 = mean(school2)
ybar3 = mean(school3)

# Variance of each sample set
var1 = var(school1)
var2 = var(school2)
var3 = var(school3)

# Sample size for each school
n1 = length(school1)
n2 = length(school2)
n3 = length(school3)

# Create function for posterior mean
post_mean = function(k_0, n, mu_0, ybar){
  
  return(
    (k_0*mu_0 + n*ybar)/(k_0 + n)
  )
}

# Posterior mean for each sample mean
postMeanM1 = post_mean(k_0, n1, mu_0, ybar1)
postMeanM2 = post_mean(k_0, n2, mu_0, ybar2) 
postMeanM3 = post_mean(k_0, n3, mu_0, ybar3)

# Posterior mean for each sample standard deviation
postMeanSD1 = post_mean(k_0, n1, mu_0, var1^.5)
postMeanSD2 = post_mean(k_0, n2, mu_0, var2^.5)
postMeanSD3 = post_mean(k_0, n3, mu_0, var3^.5)

# Create a function to calculate posterior variance
post_var = function(n, nu_0, s2_0, v, k_0, ybar, mu_0){
  return(
    ((nu_0*s2_0 + (n-1)*v + k_0*n/(k_0 + n)*(ybar - mu_0)^2)/(nu_0 + n))
  )
}

# Posterior variance for each sample
postVar1 = post_var(n1, nu_0, s2_0, var1, k_0, ybar1, mu_0)
postVar2 = post_var(n2, nu_0, s2_0, var2, k_0, ybar2, mu_0)
postVar3 = post_var(n3, nu_0, s2_0, var3, k_0, ybar3, mu_0)

# Monte Carlo procedure for post sample variance
postSampleVar1 = 1/rgamma(10000, (nu_0 + n1)/2, postVar1*(nu_0 + n1)/2)
postSampleVar2 = 1/rgamma(10000, (nu_0 + n2)/2, postVar2*(nu_0 + n2)/2)
postSampleVar3 = 1/rgamma(10000, (nu_0 + n3)/2, postVar3*(nu_0 + n3)/2)

# Create a matrix from posterior samples variances
postSampleVar = matrix(0,10000,3)
postSampleVar[,1] = postSampleVar1
postSampleVar[,2] = postSampleVar2
postSampleVar[,3] = postSampleVar3

# Monte Carlo procedure for post sample mean
postSampleTheta1 = rnorm(10000, postMeanM1, sqrt(postSampleVar1/(k_0 + n1)))
postSampleTheta2 = rnorm(10000, postMeanM2, sqrt(postSampleVar2/(k_0 + n2)))
postSampleTheta3 = rnorm(10000, postMeanM3, sqrt(postSampleVar3/(k_0 + n3)))

# Calculate confidence interval for each school's mean
cIntervalMean1 = quantile(postSampleTheta1, c(0.025,0.975))
cIntervalMean2 = quantile(postSampleTheta2, c(0.025,0.975))
cIntervalMean3 = quantile(postSampleTheta3, c(0.025,0.975))

# Calculate confidence interval for each school's standard deviation
cIntervalVar1 = quantile(postSampleVar1^.5, c(0.025,0.975))
cIntervalVar2 = quantile(postSampleVar2^.5, c(0.025,0.975))
cIntervalVar3 = quantile(postSampleVar3^.5, c(0.025,0.975))




# Part b: Compute the posterior probability that theta_i < theta_j < theta_k 
# for all six permutations {i,j,k} of {1,2,3}


# Create a matrix from posterior samples
postsample = matrix(0,10000,3)
postsample[,1] = postSampleTheta1
postsample[,2] = postSampleTheta2
postsample[,3] = postSampleTheta3

# Compute probability of each permutation
p123 = sum(postsample[,1] < postsample[,2] & postsample[,2] < postsample[,3])/10000
p132 = sum(postsample[,1] < postsample[,3] & postsample[,3] < postsample[,2])/10000
p213 = sum(postsample[,2] < postsample[,1] & postsample[,1] < postsample[,3])/10000
p231 = sum(postsample[,2] < postsample[,3] & postsample[,3] < postsample[,1])/10000
p312 = sum(postsample[,3] < postsample[,1] & postsample[,1] < postsample[,2])/10000
p321 = sum(postsample[,3] < postsample[,2] & postsample[,2] < postsample[,1])/10000




# Part C: Compute the posterior probability that ~Y_i < ~Y_j < ~Y_k 
# for all six permutations {i,j,k} of {1,2,3}

# Define posterior variables
kn = c(k_0 + n1, k_0 + n2, k_0 + n3)
mu = c(
  (k_0 * mu_0 + n1 * ybar1)/kn[1],
  (k_0 * mu_0 + n2 * ybar2)/kn[2],
  (k_0 * mu_0 + n3 * ybar3)/kn[3])

# Posterior Prediction Distribution 
postPredictDist = function(i){
  return(
    rnorm(10000, mu[i], sqrt(postSampleVar[i]*(1+1/kn[i])))
  )
}

predictPostSample = matrix(0,10000,3)


for(i in 1:3){
  predictPostSample[,i] = postPredictDist(i) #say that ten times fast!
}

# Compute probabilities for each permutation
pp123 = sum(predictPostSample[,1] < predictPostSample[,2] & predictPostSample[,2] < predictPostSample[,3])/10000
pp132 = sum(predictPostSample[,1] < predictPostSample[,3] & predictPostSample[,3] < predictPostSample[,2])/10000
pp213 = sum(predictPostSample[,2] < predictPostSample[,1] & predictPostSample[,1] < predictPostSample[,3])/10000
pp231 = sum(predictPostSample[,2] < predictPostSample[,3] & predictPostSample[,3] < predictPostSample[,1])/10000
pp312 = sum(predictPostSample[,3] < predictPostSample[,1] & predictPostSample[,1] < predictPostSample[,2])/10000
pp321 = sum(predictPostSample[,3] < predictPostSample[,2] & predictPostSample[,2] < predictPostSample[,1])/10000
