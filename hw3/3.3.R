
library(rjags)

# samples

yA=c(12,9,12,14,13,13,15,8,15,6)
yB=c(11,11,10,9,9,8,7,10,6,8,8,9,7)


# see http://www.inside-r.org/packages/cran/BRugs/docs/BRugsFit
# fit.1=BRugsFit(modelFile="model.txt",
# 	data=list(yA=yA,yB=yB),
# 	inits=list(list(thetaA=10,thetaB=10)),
# 	parametersToSave=c("thetaA","thetaB","D"),
# 	numChains=1,
# 	nIter=10000)

# fit.1=BRugsFit(modelFile="model.txt",inits=list(list(thetaA=10,thetaB=10)),data=list(yA=yA,yB=yB),numChains=1,nIter=10000,parametersToSave=c("thetaA","thetaB","D")) 


jags = jags.model('model.txt',
				data=list('yA'=yA,'yB'=yB),
				n.chains = 1,
				n.adapt = 100)
 
update(jags, 1000)
 
jags.samples(jags,c('thetaA', 'thetaB','D'),1000)
