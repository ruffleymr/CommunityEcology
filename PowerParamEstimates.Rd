######################################################################################################
## This script investigates the power to estimate the tau parameter
######################################################################################################
######################################################################################################
## 1. Load needed packages and data
library(abc)

setwd("/Users/Megan/Documents/EcolLetters_Data")
load(file="RegPars40.Rdata")
load(file="ComParams40.Rdata")
load(file="ComTraits40.Rdata")
load(file="RegTraits40.Rdata")
load(file="SumStats40.Rdata")

n <-   		   c(100, 	200,	300,		400,	500,		600)

#Local Communities
MSE.bm.filt <- c(0.026,	0.027,	0.021,	0.029,	0.023,	0.019)
MSE.bm.comp <- c(0.008, 0.0085,	0.022,	0.019,	0.014,	0.016)
MSE.ou.filt <- c(0.017,	0.017,	0.013,	0.012,	0.023,	0.014)
MSE.ou.comp <- c(0.017,	0.016, 	0.020,	0.027,	0.025,	0.018)



plot(n, MSE.bm.filt, ylim=c(0,.1), type="l", ylab="MSE of P(Acceptance)", xlab="Local Community Size", lwd=2)
points(n, MSE.bm.comp, type="l", lty=2, lwd=2)
points(n, MSE.ou.filt, type="l", lty=3, lwd=2)
points(n, MSE.ou.comp, type="l", lty=4, lwd=2)
legend("topright", c("bm.filt","bm.comp", "ou.filt", "ou.comp"), lty = c(1,2,3,4), lwd=c(3,3,3,3), cex=1.5)

head(Com.Params)

##Investigate Model 2: BM habitat filtering.
#Select 100 simulations, and perform param est. and calculate the MSE of the estimate

#get sims that have sample size less than 100 and are model 2

index <- Com.Params[as.numeric(Com.Params[,5]) > 500 & as.numeric(Com.Params[,5]) <= 600 & Com.Params[,1]=="mod6", 1]
ss <- SumStats[as.numeric(Com.Params[,5]) > 500 & as.numeric(Com.Params[,5]) <= 600 & Com.Params[,1]=="mod6", c(3,4,10,12,13,16,17,18,19,20)]
coms <- Com.Params[as.numeric(Com.Params[,5]) > 500 & as.numeric(Com.Params[,5]) <= 600 & Com.Params[,1]=="mod6", ]

psuedo.obs.index <- sample(nrow(ss), 100)
pseudo.obs.ss <- ss[psuedo.obs.index, ]

sq.error <- rep(NA, 100)
pq.mat <- matrix(NA, 100, 2)
Pa.error <- matrix(NA, 100, 3)
i=50
for (i in 1:100) {
	psuedo.abc <- abc(target=pseudo.obs.ss[i,], param=coms, sumstat=ss[-psuedo.obs.index[i],], tol=0.05, method="rejection" )
	p <- median(as.numeric(psuedo.abc$unadj.values[,10]))/100
	q <- as.numeric(coms[psuedo.obs.index[i],10])/100
	pq.mat[i,] <- c(p, q)
	sq.error[i] <- (p-q)^2


	dif <- rep(NA, length(Com.Traits[psuedo.obs.index[i]]))
	obs.comtr <- Com.Traits[[psuedo.obs.index[i]]]
	for (k in 1:length(dif)){
		dif[k] <- obs.comtr[k] - mean(obs.comtr)
	}

	Pa <- rep(NA, nrow(psuedo.abc$unadj.values))
	for (p in 1:length(Pa)){
		wght <- as.numeric(psuedo.abc$unadj.values[p,10])
		Pa[p] <- 1 - (exp(-((mean(abs(dif)))^2)/wght))
	}
	
	
	plot(as.numeric(psuedo.abc$unadj.values[,10]), Pa)
	abline(h=median(Pa), lty=2, lwd=1.5, col="red")
	Pa.obs <- 1 - (exp(-((mean(abs(dif)))^2)/as.numeric(coms[psuedo.obs.index[i],10])))
	abline(h=Pa.obs, lwd=2)
	coms[psuedo.obs.index[i],]
	
	Pa.error[i,] <- c(median(Pa), Pa.obs, (median(Pa)-Pa.obs)^2)
	
}

hist(Pa.error[,3])
mean(Pa.error[,3])
plot(Pa.error[,1], Pa.error[,2])

