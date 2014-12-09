
cat("

model
{
for(i in 1:n){
m[i]<-inprod(X[i,],beta[]) 
y[i] ~dnorm(m[i], precError)
} 
beta[1:4] ~ dmnorm(mu0[],precMV[,]) #prior on coefficients, mu0=0 
for(k in 1:4){
 for (j in 1:4){
precMV[k,j]<-XtX[j,k]*precError/n #note: this with precError below are the g-prior with g=n
# and in terms of precision, giving Sigma0 (Hoff's notation) (prec^-1)* XtX^-1 * g
}
} 

a<-nu0/2   # general parameters in gamma Hoff notation
b<-nu0*var0/2
precError ~ dgamma(a,b) #prior on 1/sigma^2 for errors 
sigma<-pow(precError,-.5) 
}
", file="model.txt")


y=c(-.87,-10.74,-3.27,-1.97,7.50,-7.25,17.05,4.96,10.40,11.05,.26,2.51)
x1=rep(1,12)
x2=c(rep(0,6),rep(1,6)) 
x3=c(23,22,22,25,27,20,31,23,27,28,22,24)
x4=x2*x3
oxy=data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4) 

X=as.matrix(oxy[,2:5])
XtX=t(X)%*%X 


# lm(y~ x1 + x2*x3,data=oxy)
# lm(y~ -1 + x1 + x2*x3,data=oxy)
lm.fit=lm(y~ x2+x3+x4,data=oxy)   #x1 is the intercept, so don't need it twice, x2*x3 gives all interactions, same as x2+x3+x4

yhat.ols=X%*%lm.fit$coefficients

library(BRugs)
brugs.fit=BRugsFit(
	modelFile="model.txt",
	data=list(y=oxy$y,n=12,XtX=XtX,nu0=1,var0=8.54,mu0=rep(0,4),X=X),
	inits=list(beta=c(0,0,0,0),precError=1),
	parametersToSave=c("beta","precError","sigma"),
	numChains=1,
	nIter=10000)

yhat.bayes=X%*%(brugs.fit$Stats[1:4,1])
mse.bayes=mean( (y-yhat.bayes)^2 ) 
mse.ols=mean((y-yhat.ols)^2 ) 
mse.ols
mse.bayes
# q()