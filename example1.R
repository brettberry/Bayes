rm(list=ls())
set.seed(2391)

posterior<-function(theta, N=50, k=20){
  
  like<-(theta^(k)*(1-theta)^(N-k)) #likelihood function: bernoulli 
  pr<-ifelse(abs(theta)<0.2, 1, 0) #prior
  out<-pr*like #according to bayes rule
  ifelse(theta>1, 0, out) #ensure not out of range
  ifelse(theta<0, 0, out)
  return(out)
}

theta.seq<-seq(from=0, to=0.4, by=0.001) #generates a sequence

post.den<-posterior(theta.seq) #passes the generated sequence through the posteroir fn

post.den.plot<-(post.den/sum(post.den))*length(theta.seq) #adjust so that it's denom=1

plot(post.den.plot~theta.seq, type="l") #plot with a smooth line

#This highlights Crowwell's Rule: ie if the prior assigns no probabity
#to an event there is no amount of evidence that can change this.
#Hence the importance of an appropriate prior

#Next run again, but with a flat (improper, but forced) prior (=5)

rm(list=ls())
set.seed(2391)

posterior<-function(theta, N=50, k=20){
  
  like<-(theta^(k))*(1-theta)^(N-k)
  pr<-5 #this time we change the posterior to an improper flat function
  out<-pr*like
  ifelse(theta>1, 0, out)
  ifelse(theta<0, 0, out)
  return(out)
}

theta.seq<-seq(from=0, to=0.8, by=0.001)

post.den<-posterior(theta.seq)

post.den.plot<-post.den*(length(theta.seq)/sum(post.den))

plot(post.den.plot~theta.seq, type="l")

#recompute with an informative prior

rm(list=ls())
set.seed(2391)

posterior<-function(theta, N=50, k=20){
  
  like<-(theta^(k))*(1-theta)^(N-k)
  pr<-dnorm(theta, mean=0.3, sd=0.05)
  out<-pr*like
  ifelse(theta>1, 0, out)
  ifelse(theta<0, 0, out)
  return(out)
}

theta.seq<-seq(from=0, to=0.8, by=0.001)

post.den<-posterior(theta.seq)

post.den.plot<-post.den*(length(theta.seq)/sum(post.den))

plot(post.den.plot~theta.seq, type="l")


#Beta prior density

theta.seq=seq(from=0, to=1, by=0.001)
plot(dbeta(theta.seq, shape1=30, shape2=30)~theta.seq, type="l")

rm(list=ls())
set.seed(2391)

posterior<-function(theta, N=50, k=20, alpha=10, beta=30{
  
  like<-(theta^(k)*(1-theta)^(N-k))
  pr<-(theta^alpha)*(1-theta)^beta
  out<-pr*like
  ifelse(theta>1, 0, out)
  ifelse(theta<0, 0, out)
  return(out)
})
  

theta.seq<-seq(from=0, to=0.8, by=0.001)

post.den<-posterior(theta.seq)

post.den.plot<-post.den*(length(theta.seq)/sum(post.den))

plot(post.den.plot~theta.seq, type="l")
lines(dbeta(theta.seq, shape1=10, shape2=30)~theta.seq, lty=2)
  
