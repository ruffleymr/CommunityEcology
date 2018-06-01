######################################################################################################
## Perform posterior predictive simulations
######################################################################################################
######################################################################################################
## 1. Determine 250 posterior samples, some are from neutral OU and some are from filtering OU
library(abc)

## load logged Summary Stats of observed data
setwd("/Users/Megan/Documents/ComAssembly_MassSims")
load(file="K.SumStatsLogged.Rdata")
load(file="Kipuka_ComTraits.Rdata")	#com.traits
load(file="Kipuka_RegTraits.Rdata")
com.traits <- log(com.traits/100)
reg.traits <- log(reg.traits/100)

##load 10000 simulations from empircal based sims that have same sample size
setwd("/Users/Megan/Documents/EcolLetters_Data/EmpiricalBased_Sims/")
load(file="ComParams.Rdata")		#Com.Params
load(file="RegTraits.Rdata")
load(file="SumStats.Rdata")			#SumStats 
load(file="RegPhyParams.Rdata") #regphy.params

setwd("/Users/Megan/Documents/EcolLetters_Data/EmpSims_52918/")
load(file="ComParams.Rdata")		#Com.Params
load(file="RegTraits.Rdata")
load(file="SumStats.Rdata")			#SumStats 
load(file="RegParams.Rdata")


## ABC
obs.ss <- rbind(K.SumStats, K.SumStats)
mod.index=Com.Params[,1]
SS <- SumStats[,c(3,4,10,12,13,16,17,18,19,20)]
obs.ss <- obs.ss[,c(3,4,10,12,13,16,17,18,19,20)]
K.modsel <- postpr(target=obs.ss[1,], index=mod.index, sumstat=SS, tol=0.025, method="rejection" )
summary(K.modsel)

##plot posterior of parameters for model 4

params.post <- Com.Params[as.numeric(rownames(K.modsel$ss)),]
params.post.reg <- reg.params[as.numeric(rownames(K.modsel$ss)),]

######May 31, Thursday, get post of Tau and P(a)
head(params.post)

#Tau Post and prior
tau.post <- as.numeric(params.post[params.post[,1]=="mod5", 5])
unif <- runif(1000, 7, 133)
plot(density(tau.post, bw=6), type="l", xlab="Tau", main="", lwd=2, bty="n")
polygon(density(tau.post, bw=6), col="skyblue")
points(density(unif, bw=8), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.5)


cp.M5 <- Com.Params[Com.Params[,1]=="mod5",]
nrow(cp.M5)
hist(as.numeric(cp.M5[,5]))

cp.M5 <- rlnorm(1660, meanlog=3.5, sd=.5)
hist(cp.M5)

cp.M5 <- rnorm(1660, mean=75, sd=20)
hist(cp.M5)


rt.M5 <- Reg.Traits[Com.Params[,1]=="mod5"]
post.rt <- Reg.Traits[as.numeric(rownames(K.modsel$ss[params.post[,1]=="mod5",]))]

##determine prior of PA
#estimate the posterior for the average probability of acceptance
mn.dif <- rep(NA, length(rt.M5))
for (i in 1:length(mn.dif)){
	ob.tr <- rt.M5[[i]]
	dif <- rep(NA, length(ob.tr))
	for (p in 1:length(dif)){
		dif[p] <- abs(ob.tr[p] - mean(ob.tr))
	}
	mn.dif[i] <- mean(dif)
}

Pa <- rep(NA, length(cp.M5))
for (i in 1:length(Pa)){
	wght <- as.numeric(cp.M5[i])
	Pa[i] <- exp(-((mn.dif[i])^2)/wght)
}


obs.mean.dif <- rep(NA, length(post.rt))
for (i in 1:length(obs.mean.dif)){
	ob.tr <- post.rt[[i]]
	dif <- rep(NA, length(ob.tr))
	for (p in 1:length(dif)){
		dif[p] <- abs(ob.tr[p] - mean(ob.tr))
	}
	obs.mean.dif[i] <- mean(dif)
}

obs.Pa <- rep(NA, length(tau.post))
for (i in 1:length(obs.Pa)){
	wght <- as.numeric(tau.post[i])
	obs.Pa[i] <- exp(-((obs.mean.dif[i])^2)/wght)
}


plot(density(obs.Pa,bw =.013), xlab="P(acceptance)", main="", lwd=2, ylim=c(0,10), bty="n")
polygon(density(obs.Pa,bw =.013), col="pink2")
abline(v=median(obs.Pa), lty=2, lwd=2)




plot(density(Pa, bw=.02), xlab="P(acceptance)", main="", lwd=2, ylim=c(0,7))
points(density(obs.Pa, bw=.02), type="l", lwd=2, lty=3)
legend("topright",  legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)













##look at sigma^2
plot(density(as.numeric(params.post[,3]), bw=0.2), xlab="sigma^2", main="", lwd=2)
legend("topright", legend=c("KS.test p-value < 0.0001"), bty="n", cex=1.3)
points(density(sample(as.numeric(Com.Params[,3]), 250)), type="l")
#abline(v=median(as.numeric(params.post[,3])))  ##so simulate with rate param between 0 and 2/ or not?
length(as.numeric(params.post[,3]))
min(as.numeric(params.post[,3]))
max(as.numeric(params.post[,3]))
sig.ks <- ks.test(as.numeric(params.post[,3]), runif(250, 0.05, 5))
sig.ks.prior <- ks.test(as.numeric(params.post[,3]), sample(as.numeric(Com.Params[,3]), 250))

#prior v. post plot
plot(density(as.numeric(Com.Params[,2]), bw=0.19), xlab="sigma^2", main="", lwd=2, ylim=c(0,.7))
points(density(as.numeric(params.post[,2]), bw=0.15), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)


##log at alpha
plot(density(as.numeric(params.post[,7]), bw=0.01), xlab="alpha", main="", lwd=2)
legend("topright", legend=c("KS.test p-value 0.16"), bty="n", cex=1.3)
#abline(v=median(as.numeric(params.post[,7])))  ##so simulate with rate param between 0 and 2/ or not?
length(as.numeric(params.post[,7]))
min(as.numeric(params.post[,7]))
max(as.numeric(params.post[,7]))
sig.ks <- ks.test(as.numeric(params.post[,7]), runif(250, 0.01, .2))
sig.ks <- ks.test(as.numeric(params.post[,7]), sample(as.numeric(Com.Params[,7]), 500))
sig.ks
plot(density(sample(as.numeric(Com.Params[,7]), 250)))

#prior v. post plot
plot(density(as.numeric(Com.Params[,7]), bw=0.009), xlab="alpha", main="", lwd=2, ylim=c(0,8))
points(density(as.numeric(params.post[,7]), bw=0.007), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)


##look at speciation
plot(density(as.numeric(params.post.reg[,2]), bw=0.01), xlab="speciation", main="", lwd=2)
legend("topright", legend=c("KS.test p-value < 0.0001"), bty="n", cex=1.3)
#abline(v=median(as.numeric(params.post[,7])))  ##so simulate with rate param between 0 and 2/ or not?
length(as.numeric(params.post.reg[,2]))
min(as.numeric(params.post.reg[,2]))
max(as.numeric(params.post.reg[,2]))
sig.ks <- ks.test(as.numeric(params.post.reg[,2]), runif(250, 0.05, .2))
sig.ks <- ks.test(as.numeric(params.post.reg[,2]), sample(as.numeric(regphy.params[,2]), 250))
plot(density(sample(as.numeric(regphy.params[,2]), 250)))
plot(density(runif(250, 0.05, .2)))

#prior v. post plot
plot(density(as.numeric(regphy.params[,2]), bw=0.008), xlab="speciation", main="", lwd=2, ylim=c(0,17))
points(density(as.numeric(params.post.reg[,2]), bw=0.008), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)


##look at extinction
plot(density(as.numeric(params.post.reg[,4]), bw=0.001), xlab="speciation", main="", lwd=2)
legend("topright", legend=c("KS.test p-value < 0.0001"), bty="n", cex=1.3)
#abline(v=median(as.numeric(params.post[,7])))  ##so simulate with rate param between 0 and 2/ or not?
length(as.numeric(params.post.reg[,2]))
min(as.numeric(params.post.reg[,2]))
max(as.numeric(params.post.reg[,2]))
sig.ks <- ks.test(as.numeric(params.post.reg[,2]), runif(250, 0.05, .2))
sig.ks <- ks.test(as.numeric(params.post.reg[,4]), sample(as.numeric(regphy.params[,4]), 250), alternative="greater")
sig.ks

#speciation fraction
plot(density(as.numeric(regphy.params[,3]), bw=0.02), xlab="extinction fraction", main="", lwd=2, ylim=c(0,4))
points(density(as.numeric(params.post.reg[,3]), bw=0.03), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)


#############Make Conensus plot
par(mfrow=c(2,2))


params.M5 <- params.post[params.post[,1]=="mod5",]
cp.M5 <- Com.Params[Com.Params[,1]=="mod5",]
rt.M5 <- Reg.Traits[Com.Params[,1]=="mod5"]
post.rt <- Reg.Traits[as.numeric(rownames(K.modsel$ss))]

plot(density(as.numeric(cp.M5[,6]), bw=0.011), xlab="Sf", main="", lwd=2, ylim=c(0,10))
points(density(as.numeric(params.M5[,6]), bw=0.01), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)

plot(density(as.numeric(cp.M5[,9]), bw=10), xlab="Tau.Constant (Sigma^2 * TreeDepth)", main="", lwd=2, ylim=c(0,0.006))
points(density(as.numeric(params.M5[,9]), bw=20), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)


plot(density(as.numeric(cp.M5[,5]), bw=5), xlab="Tau", main="", lwd=2)
points(density(as.numeric(params.M5[,5]), bw=5), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)

##determine prior of PA
#estimate the posterior for the average probability of acceptance
mn.dif <- rep(NA, length(rt.M5))
for (i in 1:length(mn.dif)){
	ob.tr <- rt.M5[[i]]
	dif <- rep(NA, length(ob.tr))
	for (p in 1:length(dif)){
		dif[p] <- abs(ob.tr[p] - mean(ob.tr))
	}
	mn.dif[i] <- mean(dif)
}

Pa <- rep(NA, nrow(cp.M5))
for (i in 1:length(Pa)){
	wght <- as.numeric(cp.M5[i,5])
	Pa[i] <- exp(-((mn.dif[i])^2)/wght)
}


obs.mean.dif <- rep(NA, length(post.rt))
for (i in 1:length(obs.mean.dif)){
	ob.tr <- post.rt[[i]]
	dif <- rep(NA, length(ob.tr))
	for (p in 1:length(dif)){
		dif[p] <- abs(ob.tr[p] - mean(ob.tr))
	}
	obs.mean.dif[i] <- mean(dif)
}

obs.Pa <- rep(NA, nrow(params.post))
for (i in 1:length(obs.Pa)){
	wght <- as.numeric(params.post[i,5])
	obs.Pa[i] <- exp(-((obs.mean.dif[i])^2)/wght)
}
obs.Pa <- obs.Pa[params.post[,1]=="mod5"]


plot(density(Pa, bw=.02), xlab="P(acceptance)", main="", lwd=2, ylim=c(0,7))
points(density(obs.Pa, bw=.02), type="l", lwd=2, lty=3)
legend("topright",  legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)


plot(density(mn.dif, bw=.3), xlab="mean pairwise difference between regional traits", main="", lwd=2, ylim=c(0,.5))
points(density(obs.mean.dif, bw=.3), type="l", lwd=2, lty=3)
legend("topright",  legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)


plot(params.M5[,3], obs.Pa)

par(mfrow=c(2,2))

plot(density(as.numeric(Com.Params[,3]), bw=0.19), xlab="sigma^2", main="", lwd=2, ylim=c(0,.7))
points(density(as.numeric(params.post[,3]), bw=0.15), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)

plot(density(as.numeric(Com.Params[,7]), bw=0.009), xlab="alpha", main="", lwd=2, ylim=c(0,8))
points(density(as.numeric(params.post[,7]), bw=0.007), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)

plot(density(as.numeric(regphy.params[,2]), bw=0.008), xlab="speciation", main="", lwd=2, ylim=c(0,17))
points(density(as.numeric(params.post.reg[,2]), bw=0.008), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)

plot(density(as.numeric(regphy.params[,3]), bw=0.02), xlab="extinction fraction", main="", lwd=2, ylim=c(0,4))
points(density(as.numeric(params.post.reg[,3]), bw=0.03), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)






params.M4 <- params.post[params.post[,1]=="mod4",]

##look at sigma^2
plot(density(as.numeric(params.M4[,3]), bw=0.2), breaks=50, xlab="sigma^2", main="", lwd=2)
legend("topright", legend=c("KS.test p-value < 0.0001"), bty="n", cex=1.3)
abline(v=median(as.numeric(params.M4[,3])))  ##so simulate with rate param between 0 and 2/ or not?
length(as.numeric(params.M4[,3]))
min(as.numeric(params.M4[,3]))
max(as.numeric(params.M4[,3]))
sig.ks <- ks.test(as.numeric(params.M4[,3]), runif(114, 0.05, 5)


#look at alpha
hist(as.numeric(params.M4[,7]), breaks=50)
abline(v=median(as.numeric(params.M4[,7])))  ##I guess keep between 0.01 and 0.2

##plot posterior of parameters for model 5
params.M5 <- params.post[params.post[,1]=="mod5",]

##look at sigma^2
hist(as.numeric(params.M5[,3]), breaks=50)
abline(v=median(as.numeric(params.M5[,3])))  ##also unclear

#look at alpha
hist(as.numeric(params.M5[,7]), breaks=50)
abline(v=median(as.numeric(params.M5[,7])))  ##I guess keep between 0.01 and 0.2

#look at tau
hist(as.numeric(params.M5[,10]), breaks=10, main="", xlab=" Tau ")
abline(v=median((as.numeric(params.M5[,10]))), col="red", lwd=1.5, lty=2)

#look at tau logged
hist(log(as.numeric(params.M5[,10])), breaks=10, main="", xlab="log( Tau )")
abline(v=log(median((as.numeric(params.M5[,10])))), col="red", lwd=1.5, lty=2)

plot(density(as.numeric(params.M5[,10])), lwd=3, main="", bty="l", xlab="Tau", ylim=c(0,0.03))
polygon(density(as.numeric(params.M5[,10])), col="skyblue")
abline(v=median(as.numeric(params.M5[,10])), lty=2, lwd=2)
points(density(as.numeric(cp.M5[,10]), bw=5), lwd=2, type="l", lty=3)
legend("topright", legend=c("simulated", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.5)


#estimate the posterior for the average probability of acceptance
dif<- rep(NA, length(reg.traits))
for (i in 1:length(dif)){
	dif[i] <- reg.traits[i] - mean(reg.traits)
}

Pa <- rep(NA, nrow(params.M5))
for (i in 1:length(Pa)){
	wght <- as.numeric(params.M5[i,10])
	Pa[i] <- (exp(-((mean(abs(dif)))^2)/wght))
}


plot(density(as.numeric(params.M5[,10])), lwd=3, main="", bty="l", xlab="Tau", ylim=c(0,0.03))
polygon(density(as.numeric(params.M5[,10])), col="skyblue")
abline(v=median(as.numeric(params.M5[,10])), lty=2, lwd=2)
points(density(as.numeric(cp.M5[,10]), bw=5), lwd=2, type="l", lty=3)
legend("topright", legend=c("simulated", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.5)

plot(as.numeric(params.M5[,6]), obs.Pa )
plot(as.numeric(params.post[params.post[,1]=="mod5",3]), obs.Pa )
plot(as.numeric(Com.Params[Com.Params[,1]=="mod5",3]), Pa)
plot(as.numeric(cp.M5[,6]), Pa)
plot(lm(obs.Pa ~ as.numeric(params.post[params.post[,1]=="mod5",3]) + as.numeric(params.M5[,6])))


plot(density(obs.Pa, bw = .03), lwd=3, main="", xlab="Mean P(acceptance)",  ylim=c(0,8.5), xlim=c(-.1,1.05), bty="n")
polygon(density(obs.Pa, bw = .03), col="pink2")
points(density(Pa, bw=.03), xlab="P(acceptance)", lwd=2, lty=3, type="l")
legend("topright", legend=c("simulated", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.3)
segments(x0=median(obs.Pa), y0=0, x1=median(obs.Pa), y1=6.7, lwd=2, lty=2)


obs.mean.dif <- obs.mean.dif[params.post[,1]=="mod5"]
plot(density(obs.mean.dif, bw=.15), type="l", lwd=3, main="", xlab="Regional Traits Mean Pairwise Differences",
			ylim=c(0,0.7), xlim=c(0,12))
polygon(density(obs.mean.dif, bw = .15), col="darkseagreen3")
points(density(mn.dif, bw=.15), lwd=2, type="l", lty=3)
segments(x0=median(obs.mean.dif), y0=0, x1=median(obs.mean.dif), y1=0.6, lwd=2, lty=2)
legend("topright", legend=c("simulated", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.5)


plot(density(as.numeric(params.post[params.post[,1]=="mod5",3]), bw=0.2), type="l", lwd=2,  xlab="sigma^2", main="", ylim=c(0,.7), bty="n", xlim=c(0,6))
polygon(density(as.numeric(params.post[params.post[,1]=="mod5",3]), bw=0.2), col="darkseagreen3")
points(density(as.numeric(Com.Params[Com.Params[,1]=="mod5",3]), bw=0.19), type="l", lty=3, lwd=2)
segments(x0=median(as.numeric(params.post[params.post[,1]=="mod5",3])), y0=0, x1=median(as.numeric(params.post[params.post[,1]=="mod5",3])), y1=0.7, lwd=2, lty=2)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.5)


plot(density(as.numeric(params.post[,3]), bw=0.15), type="l", lwd=2,  xlab="sigma^2", main="", ylim=c(0,.7), bty="n")
polygon(density(as.numeric(params.post[,3]), bw=0.15), col="darkseagreen3")
points(density(as.numeric(Com.Params[,3]), bw=0.19), type="l", lty=3, lwd=2)
segments(x0=median(as.numeric(params.post[,3])), y0=0, x1=median(as.numeric(params.post[,3])), y1=0.7, lwd=2, lty=2)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.5)


plot(density(as.numeric(params.M5[,6]), bw=0.013), type="l", xlab="Sf", main="", lwd=2, ylim=c(0,10), bty="n")
polygon(density(as.numeric(params.M5[,6]), bw=0.013), col="skyblue")
points(density(as.numeric(cp.M5[,6]), bw=0.013), type="l", lwd=2, lty=3)
legend("topright", legend=c("prior", "posterior"), lwd=c(2,2), lty=c(3,1), bty="n", cex=1.5)


plot(log(as.numeric(params.M5[,10])), Pa, xlab="log( Tau )", ylab="P(Acceptance into Local)")
pa.med <- (exp(-((mean(abs(dif)))^2)/median(as.numeric(params.M5[,10]))))
abline(h=pa.med, col="red", lty=2, lwd=1.5)
points(log(median(as.numeric(params.M5[,10]))), pa.med, bg="red", pch=21)  ##limit to make sure prob is between .6 and 1

range(as.numeric(params.M5[,7]))
######################################################################################################
##Simulate data for post predicitve sims
setwd("/Users/Megan/Documents/EcolLetters_Data/EmpiricalBased_Sims/")
load(file="ComParams5.Rdata")		#Com.Params
load(file="SumStats5.Rdata")

##ABC param estimate
obs.ss <- rbind(K.SumStats, K.SumStats)
mod.index=Com.Params[,1]
SS <- SumStats[,c(3,4,10,12,13,16,17,18,19,20)]
obs.ss <- obs.ss[,c(3,4,10,12,13,16,17,18,19,20)]
Kpp.modsel <- abc(target=obs.ss[1,], param=Com.Params, sumstat=SS, tol=0.025, method="rejection" )
attributes(Kpp.modsel)

Kpp.modsel$ss

mod.fit <- gfit(target=obs.ss[1,], sumstat=Kpp.modsel$ss, nb.replicate=100, tol=.01, statistic=median, subset=NULL, trace=FALSE)

summary(mod.fit)


Kpp.modsel$unadj.values
Kpp.modsel$ss

hist(Kpp.modsel$ss[,1])
abline(v=obs.ss[1,1], lwd=3)

hist(Kpp.modsel$ss[,2])
abline(v=obs.ss[1,2], lwd=3)

hist(Kpp.modsel$ss[,3])
abline(v=obs.ss[1,3], lwd=3)

hist(Kpp.modsel$ss[,4])
abline(v=obs.ss[1,4], lwd=3)

hist(Kpp.modsel$ss[,5])
abline(v=obs.ss[1,5], lwd=3)

hist(Kpp.modsel$ss[,6])
abline(v=obs.ss[1,6], lwd=3)

hist(Kpp.modsel$ss[,7])
abline(v=obs.ss[1,7], lwd=3)

hist(Kpp.modsel$ss[,8])
abline(v=obs.ss[1,8], lwd=3)

hist(Kpp.modsel$ss[,9])
abline(v=obs.ss[1,9], lwd=3)

hist(Kpp.modsel$ss[,10])
abline(v=obs.ss[1,10], lwd=3)

fit.data <- rbind(obs.ss[1,], Kpp.modsel$ss)
pca <- prcomp(fit.data)
plot(pca$x[,1], pca$x[,2])
points(pca$x[1,1], pca$x[1,2], pch=15, col="red")


hist(log(as.numeric(Kpp.modsel$unadj.values[,10])), breaks=20, main="", xlab="Log(Tau)")
hist(log(as.numeric(params.M5[,10])), breaks=10, main="", xlab="log( Tau )", add=T, col="red")
legend("topright", c("posterior", "post. preditive"), fill=c("red", "white"))

hist(as.numeric(Kpp.modsel$unadj.values[,7]))


Pa <- rep(NA, nrow(Kpp.modsel$unadj.values))
for (i in 1:length(Pa)){
	wght <- as.numeric(Kpp.modsel$unadj.values[i,10])
	Pa[i] <- (exp(-((mean(abs(dif)))^2)/wght))
}

plot(log(as.numeric(Kpp.modsel$unadj.values[,10])), Pa, xlab="log( Tau )", ylab="P(Acceptance into Local)", cex=1.5)
pa.med <- (exp(-((mean(abs(dif)))^2)/median(as.numeric(Kpp.modsel$unadj.values[,10]))))
abline(h=pa.med, col="black", lty=2, lwd=1.5)

points(log(as.numeric(params.M5[,10])), Pa.obs,  col="red", cex=.8)
pa.med.obs <- (exp(-((mean(abs(dif)))^2)/median(as.numeric(params.M5[,10]))))
abline(h=pa.med.obs, col="red", lty=2, lwd=1.5)
legend("bottomright", c("Posterior", "Post.Predictive"), col=c("red","black"), pch=c(1,1), pt.cex=c(.8, 1.5))



######################################################################################################
##Simulate data for post predicitve sims, check and see how sum stats are coming out. 
setwd("/Users/Megan/Documents/EcolLetters_Data/EmpiricalBased_Sims/")
load(file="ComParams5.Rdata")		#Com.Params
load(file="SumStats5.Rdata")

K.modsel
ss <- K.modsel$ss[K.modsel$values=="mod5",]
head(ss)
sample.index <- sample(seq(1,10000), 72)
PP.ss <- SumStats[1:1000,c(3,4,10,12,13,16,17,18,19,20)]

head(PP.ss)

#mean com trait
plot(density(PP.ss[,1]), ylim=c(0,2))
points(density(ss[,1]), type="l")

#var com trait
plot(density(PP.ss[,2]), ylim=c(0,2))
points(density(ss[,2]), type="l")

#Msig
plot(density(PP.ss[,3]), ylim=c(0,3))
points(density(ss[,3]), type="l")

#Svar
plot(density(PP.ss[,4]), ylim=c(0,60), xlim=c(-1.5, .5))
points(density(ss[,4]), type="l")

#Shgt
plot(density(PP.ss[,5]), ylim=c(0,50), xlim=c(-.05, 2))
points(density(ss[,5]), type="l")

#krts
plot(density(PP.ss[,6]), ylim=c(0,1))
points(density(ss[,6]), type="l")

#mean reg traits
plot(density(PP.ss[,7]), ylim=c(0,1))
points(density(ss[,7]), type="l")

#var reg traits
plot(density(PP.ss[,8]),  ylim=c(0,4), xlim=c(0,50))
points(density(ss[,8]), type="l")

#mean trait dif
plot(density(PP.ss[,9]), ylim=c(0,1))
points(density(ss[,9]), type="l")

#var trait dif
plot(density(PP.ss[,10]), ylim=c(0,1))
points(density(ss[,10]), type="l")

