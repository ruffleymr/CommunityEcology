######################################################################################################
## This script determines the power of RF, ABC, and MPD/MNTD at infering community assembly models
######################################################################################################
######################################################################################################
## 1. Load needed packages and data
library(picante)
library(randomForest)
library(abc)
library(geiger)

setwd("/Users/Megan/Documents/EcolLetters_Data")
load(file="RegPars40.Rdata")
load(file="RegTrees40.Rdata")
load(file="ComParams40.Rdata")
load(file="SumStats40.Rdata")
load(file="PhyDisp40.Rdata")
load(file="TraitMod40.Rdata")


setwd("/Users/Megan/Documents/EcolLetters_Data/Sims_52318/")
load(file="ComTrees.Rdata")
load(file="ComParams.Rdata")
load(file="RegTraits.Rdata")
load(file="ComTraits.Rdata")
load(file="SumStats.Rdata")
load(file="rf.Rdata")

index <- Com.Params[,1]
ref <- data.frame(SumStats, index)
rf <- randomForest(index ~., data=ref, ntree=1000, importance=T)
rf
######################################################################################################
##2. Run Random Forest  
n <- seq(200, 600, by=50)
sims <- c()
rf.err <- c()
range(as.numeric(Com.Params[,4])
for (i in n){
	if (i == 600){
		index <- Com.Params[as.numeric(Com.Params[,4]) > i-50, 1]
		ss <- SumStats[as.numeric(Com.Params[,4]) > i-50, ]
	}else{
		index <- Com.Params[as.numeric(Com.Params[,4]) > i-50 & as.numeric(Com.Params[,4]) <= i, 1]
		ss <- SumStats[as.numeric(Com.Params[,4]) > i-50 & as.numeric(Com.Params[,4]) <= i, ]
	}
	ref <- data.frame(ss, index)
	rf <- randomForest(index ~., data=ref, ntree=500, importance=T)
	rf.err <- c(rf.err, rf$err.rate[500,1])
	sims <- c(sims, length(index))
	print(i)
	print(length(index))
	print(rf$err.rate[500,1])
}

##Error Rate Plot
plot(n, rf.err, ylab="Error Rate %", xlab="Taxa in Local Community",
		pch=16, col="black", type="b", lwd=1.3, cex=1.2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,.5))
axis(side=1, at=seq(0,600, by=50))
axis(side=2, at=seq(0.0,0.5, by=0.1))

##Power Plot
plot(n, 1-rf.err, ylab="Power", xlab="Taxa in Local Community",
		pch=16, col="black", type="b", lwd=1.3, cex=1.2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=50))
axis(side=2, at=seq(0.0,1, by=0.2))

######################################################################################################
##3. Run ABC  

abc.power <- c()
abc.mpp <- c()
abc.mpp.whenT <- c()
n <- seq(200, 600, by=50)

for (i in n){
	if (i == 600){
		index <- Com.Params[as.numeric(Com.Params[,4]) > i-50, 1]
		ss <- SumStats[as.numeric(Com.Params[,4]) > i-50, c(3,4,10,12,13,16,17,18,19,20)]
	}else{
		index <- Com.Params[as.numeric(Com.Params[,4]) > i-50 & as.numeric(Com.Params[,4]) <= i, 1]
		ss <- SumStats[as.numeric(Com.Params[,4]) > i-50 & as.numeric(Com.Params[,4]) <= i,
																	c(3,4,10,12,13,16,17,18,19,20)]
	}
	nval <- length(index)/6 - (length(index)/6*.2)
	print(paste("nval = ", nval, sep=""))
	cv <- cv4postpr(index, ss, nval=nval, tols=.01, method="rejection") 
	mods <- c()
	mpp <- c()
	mpp.t <- c()
	for (k in 1:6){
		mat <- cv$model.probs$tol0.01[rownames(cv$model.probs$tol0.01)==paste("mod", k, sep=""),k]
		est <- cv$estim$tol0.01[names(cv$estim$tol0.01)==paste("mod", k, sep="")]
		
		mods <- c(mods, sum(est==paste("mod", k, sep=""))/nval)
		mpp <- c(mpp, mean(mat))
		mpp.t <- c(mpp.t, mean(mat[est==paste("mod", k, sep="")]))	
	}
	abc.power <- c(abc.power, mean(mods))
	abc.mpp <- c(abc.mpp, mean(mpp))
	abc.mpp.whenT <- c(abc.mpp.whenT, mean(mpp.t))
	
	print(i)
	print(paste("sims =  ",length(index), sep=""))
	print(paste("Power = ",mean(mods), sep=""))
	print(paste("MPP = ",mean(mpp), sep=""))
	print(paste("MPP.true = ", mean(mpp.t), sep=""))
}

points(n, abc.power, ylab="Power", xlab="Taxa in Local Community",
		pch=19, col="black", type="b", lwd=1.3, cex=1.2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=50))
axis(side=2, at=seq(0.0,1, by=0.2))

plot(n, 1-abc.power, ylab="Error Rate %", xlab="Taxa in Local Community",
		pch=16, col="black", type="b", lwd=1.3, cex=1.2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,.5))
axis(side=1, at=seq(0,600, by=50))
axis(side=2, at=seq(0.0,0.5, by=0.1))

plot(n, abc.mpp.whenT, ylab="Power", xlab="Taxa in Local Community",
		pch=16, col="black", type="b", lwd=1.3, cex=1.2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=50))
axis(side=2, at=seq(0.0,1, by=0.2))

abc.mpp.whenT

######################################################################################################
## 4. MPD and MNTD Power

phy.power <- c()
tr.mod <- c("BM", "BM", "BM", "OU", "OU", "OU")

for (i in n){
	if (i == 600){
		p.mat <- PhyDisp[as.numeric(Com.Params[,5]) > i-50, ]
		true <- Com.Params[as.numeric(Com.Params[,5]) > i-50, 1]
		tr <- trait.mod[as.numeric(Com.Params[,5]) > i-50, 1]
	}else{
		p.mat <- PhyDisp[as.numeric(Com.Params[,5]) > i-50 & as.numeric(Com.Params[,5]) <= i, ]
		true <- Com.Params[as.numeric(Com.Params[,5]) > i-50 & as.numeric(Com.Params[,5]) <= i, 1]
		tr <- trait.mod[as.numeric(Com.Params[,5]) > i-50 & as.numeric(Com.Params[,5]) <= i, 1]
	}
	mods <-c()
	for (k in 1:6){
		 if (k == 1 | k == 4){
		 	est <- sum(true==paste("mod",k,sep="") & tr == tr.mod[k] & p.mat[,5]>0.05 & p.mat[,5]<0.95 & p.mat[,6]>0.05 & p.mat[,6]<0.95)
		 }
		 if (k == 2 | k == 5){
		 	est <- sum(true==paste("mod",k,sep="") & tr == tr.mod[k] & (p.mat[,5]<=0.05 | p.mat[,6]<=0.05))
		 }
		 if (k == 3 | k == 6){
		 	est <- sum(true==paste("mod",k,sep="") & tr == tr.mod[k] & (p.mat[,5]>=0.95 | p.mat[,6]>=0.95))
		 }
		correct.prop <- est/sum(true==paste("mod",k,sep="")) 
		mods <- c(mods, correct.prop)
	}
	phy.power <- c(phy.power, mean(mods))
}

plot(n, phy.power, ylab="Power", xlab="Taxa in Local Community",
		pch=16, col="black", type="l", lwd=1.7, cex=2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=100))
axis(side=2, at=seq(0.0,1, by=0.2))



######################################################################################################
## 5. Make Figures

## RF, ABC, and MPD/MNTD Power WITH POINTS
plot(n, 1-rf.err, ylab="Power", xlab="Taxa in Local Community",
		pch=16, col="black", type="b", lwd=1.7, cex=1, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=100))
axis(side=2, at=seq(0.0,1, by=0.2))
points(n, abc.power, ylab="Power", xlab="Taxa in Local Community",
		pch=0, col="black", type="b", lwd=1.7, cex=1, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, ylim=c(0,1))
points(n, phy.power, ylab="", xlab="",
		pch=4, col="black", type="b", lwd=1.7, cex=1, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, ylim=c(0,1))
		 

## RF, ABC, and MPD/MNTD Power WITH LINES
plot(n, 1-rf.err, ylab="Power", xlab="Taxa in Local Community",
		 col="black", type="l", lwd=3.5, cex=2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1.2, axes=FALSE, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=100))
axis(side=2, at=seq(0.0,1, by=0.2))
points(n, abc.power, ylab="Power", xlab="Taxa in Local Community",
		col="black", type="l", lwd=3.5, cex=2, lty=3, bty="o", 
		 cex.axis=1, cex.lab=1, ylim=c(0,1))
points(n, phy.power, ylab="", xlab="",
		col="black", type="l", lwd=3.5, cex=2, lty=5, bty="o", 
		 cex.axis=1, cex.lab=1, ylim=c(0,1))
legend(x=430, y=.2, c("RF", "ABC", "Mpd/Mntd"), lty=c(1,3,5), lwd=2.3, cex=1)

		 
## RF, ABC, and MPD/MNTD Error Rates WITH LINES
plot(n, rf.err, ylab="Error Rate", xlab="Taxa in Local Community",
		 col="black", type="l", lwd=3.5, cex=2, lty=1, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=100))
axis(side=2, at=seq(0.0,1, by=0.2))
points(n, 1-abc.power, ylab="Power", xlab="Taxa in Local Community",
		col="black", type="l", lwd=3.5, cex=2, lty=3, bty="o", 
		 cex.axis=1, cex.lab=1, ylim=c(0,1))
points(n, 1-phy.power, ylab="", xlab="",
		col="black", type="l", lwd=3.5, cex=2, lty=5, bty="o", 
		 cex.axis=1, cex.lab=1, ylim=c(0,1))
legend("topright", c("RF", "ABC", "Mpd/Mntd"), lty=c(1,3,5), lwd=2, cex=1.2)

		
## Plot MPPwhen true
plot(n, abc.mpp.whenT, ylab="Mean Posterior Probability", xlab="Taxa in Local Community",
		col="black", type="l", lwd=3.4, cex=1.2, lty=3, bty="o", 
		 cex.axis=1, cex.lab=1, axes=FALSE, ylim=c(0,1))
axis(side=1, at=seq(0,600, by=100))
axis(side=2, at=seq(0.0,1, by=0.2))
