######################################################################################################
## This script determines the power of RF, ABC, and MPD/MNTD at infering community assembly models
######################################################################################################
######################################################################################################
## 1. Load needed packages and data
library(diversitree)
library(randomForest)
library(abc)
library(picante)
library(phytools)

setwd("/Users/Megan/Documents/ComAssembly_MassSims")
load(file="Kipuka_ComTree.Rdata")  #com.tre
load(file="Kipuka_ComTraits.Rdata")	#com.traits
load(file="Kipuka_RegTree.Rdata") #reg.tre
load(file="Kipuka_RegTraits.Rdata") #reg.traits

load(file="K.SumStatsLogged.Rdata")
load(file="K.SumStats.Rdata")
load(file="K.SumStats.noGym.Rdata")

setwd("/Users/Megan/Documents/EcolLetters_Data/EmpiricalBased_Sims/")
load(file="ComParams.Rdata")		#Com.Params
load(file="SumStats.Rdata")			#SumStats
load(file="rf.Rdata")

com.traits <- log(com.traits/100)
reg.traits <- log(reg.traits/100)


com.traits <- com.traits/100
reg.traits <- reg.traits/100

######################################################################################################
## 2. Run RF, ABC, MPD, MNTD

## 2.1 RF
obs.ss <- rbind(K.SumStats, K.SumStats)
rownames(obs.ss) <- c("ss1", "ss2")
predict(rf, obs.ss, type="prob") 
varImpPlot(rf)

## 2.2 ABC
obs.ss <- rbind(K.SumStats, K.SumStats)
mod.index=Com.Params[,1]
SS <- SumStats[,c(3,4,10,12,13,16,17,18,19,20)]
obs.ss <- obs.ss[,c(3,4,10,12,13,16,17,18,19,20)]

SS <- SumStats[,c(4,12,13,18,19,20)]
obs.ss <- obs.ss[,c(4,12,13,18,19,20)]

K.modsel <- postpr(target=obs.ss[1,], index=mod.index, sumstat=SS, tol=0.02, method="rejection" )
summary(K.modsel)

K.modsel <- postpr(target=obs.ss[1,], index=mod.index, sumstat=SS, tol=0.05, method="rejection" )
summary(K.modsel)


SS <- SumStats[,c(3,4,10,12,13,16,17,18,19,20)]
obs.ss <- obs.ss[,c(3,4,10,12,13,16,17,18,19,20)]
K.modsel <- abc(target=obs.ss[1,], param=Com.Params, sumstat=SS, tol=0.025, method="rejection" )
attributes(K.modsel)

hist(as.numeric(K.modsel$unadj.values[,10]))


reg.tre <- drop.tip(reg.tre, tip=c("Pseudotsuga_menziesii","Juniperus_scopulorum","Pinus_flexilis"))

## 2.3 & 2.4
mat <- matrix(NA, 1, length(reg.tre$tip.label))
colnames(mat) <- reg.tre$tip.label
for (j in 1:length(reg.tre$tip.label)) {
	if (is.element(reg.tre$tip.label[j], com.tre$tip.label)){
		mat[j] <- 1
	} else {
		mat[j] <- 0 
	}
}
mat <- rbind(mat, mat)
ses.mpd <- ses.mpd(mat, cophenetic(reg.tre), runs=1000, iterations=1, null.model="taxa.labels")[1,]
ses.mntd <- ses.mntd(mat, cophenetic(reg.tre), runs=1000, iterations=1, null.model="taxa.labels")[1,]
ses.mpd	
ses.mntd

## 2.5 Evaluate models of trait evoltion

bm.fit.com <- fitContinuous(com.tre, dat=com.traits, model="BM", control=list(niter=100))
ou.fit.com <- fitContinuous(com.tre, dat=com.traits, model="OU", control=list(niter=100), bounds=list(alpha = c(min = exp(-500), max = .3)))

bm.fit.reg <- fitContinuous(reg.tre, dat=reg.traits, model="BM", control=list(niter=100))
ou.fit.reg <- fitContinuous(reg.tre, dat=reg.traits, model="OU", control=list(niter=100), bounds=list(alpha = c(min = exp(-500), max = .3)))

var(com.traits)
var(reg.traits)
var(reg.traits[names(reg.traits)!="Pseudotsuga_menziesii"])
hist(reg.traits)

######################################################################################################
## 3. asses model fit and plot

SS.mod2 <- SS[mod.index=="mod5",]

SS.mod2 <- K.modsel$ss[K.modsel$values=="mod5",]
nrow(SS.mod2)
mod2.fit <- gfit(target=obs.ss[1,], sumstat=SS.mod2, nb.replicate=100, tol=.01, statistic=median, subset=NULL, trace=FALSE)

hist(SS.mod2[,2])
abline(v=obs.ss[1,2])


summary(mod2.fit) ##ang+gym, only k.modsel, logged: mod4:0, mod5:0
 ##ang+gym, logged: p-val = 0.66, 0.75, 0.445, 0.40, 0.52, 0.38
 ## Ang+GYm, and not logged:p-value = 0.045, 0.82, 0.11, 0.0, 0.005, 0.023
attributes(mod2.fit)


plot(density(mod2.fit$dist.sim), main="Angiosperms+Gymnosperms", xlab="", bty="n",
		 col="black", lwd=3.5, cex=2, lty=1, bty="o", 
		  cex.lab=1.2, ylim=c(0,.6), xlim=c(2,15), axes=FALSE)
axis(side=1, at=seq(2,16, by=4))
axis(side=2, at=seq(0.0 ,0.6, by=0.2))
abline(v=mod2.fit$dist.obs, col="red", lwd=4, lty=3)
legend(x=7, y=.4, legend=c("Model Fit", "P-value = 0.82"), bty="n", cex=1.6)

attributes(K.modsel)
fit.data <- rbind(obs.ss[1,],K.modsel$ss[K.modsel$values=="mod4",])

, K.modsel$ss[K.modsel$values=="mod5",])

fit.data <- rbind(obs.ss[1,],K.modsel$ss[K.modsel$values=="mod5",])

fit.data <- rbind(obs.ss[1,], K.modsel$ss)


hist(fit.data[,1])
abline(v=obs.ss[1,1], col="red")
hist(fit.data[,2])
abline(v=obs.ss[1,2], col="red")
hist(fit.data[,3])
abline(v=obs.ss[1,3], col="red")
hist(fit.data[,4])
abline(v=obs.ss[1,4], col="red")
hist(fit.data[,5])
abline(v=obs.ss[1,5], col="red")
hist(fit.data[,6])
abline(v=obs.ss[1,6], col="red")
x11()
pca <- prcomp(fit.data)
attributes(pca)
plot(pca$x[,1], pca$x[,2], ylab="PC2", xlab="PC1", pch=16)
points(pca$x[1,1], pca$x[1,2], bg="red", pch=23, cex=1.6)

plot(pca$x[,3], pca$x[,2], ylab="PC2", xlab="PC1", pch=16)
points(pca$x[1,3], pca$x[1,2], bg="red", pch=23, cex=1.6)
######################################################################################################
## 4. Plot Phylogeny with traits

reg.traits[reg.traits>6] <-6

color <- c()
for (j in length(reg.tre$tip.label):1) {
	if (is.element(reg.tre$tip.label[j], com.tre$tip.label)){
		color <- c("grey10", color)
	} else {
		color <- c("white", color)
	}
}
sum(color=="black")

plotTree.barplot(tree=reg.tre, x=reg.traits, args.plotTree=list(edge.width=2, fsize=.3, type="fan", lwd=1.2, ftype="bi"), args.barplot=list(col=color, xlab="height", border="black"))

plot.phylo(reg.tre, type="fan", fsize=.2)
plotTree.wBars(reg.tre, reg.traits , fsize=.4, type="fan", tip.labels=F, scale=30, lwd=.6, col=color, 
				args.plotTree=list(adj=.1))




######################################################################################################
##figure out priors for predictive simulations

params.post <- Com.Params[as.numeric(rownames(K.modsel$ss)),]
params.M4 <- params.post[params.post[,1]=="mod4",]
params.M5 <- params.post[params.post[,1]=="mod5",]


K.modsel$unadj.values
hist(as.numeric(K.modsel$unadj.values[,3]))


hist(log(as.numeric(K.modsel$unadj.values[,10])))
abline(v=log(median((as.numeric(K.modsel$unadj.values[,10])))), col="red", lwd=1.5, lty=2)

dif <- rep(NA, length(com.traits))
for (i in 1:length(dif)){
	dif[i] <- com.traits[i] - mean(com.traits)
}


Pa <- rep(NA, nrow(K.modsel$unadj.values))
for (i in 1:length(Pa)){
	wght <- as.numeric(K.modsel$unadj.values[i,10])
	Pa[i] <- (exp(-((mean(abs(dif)))^2)/wght))
}

plot(log(as.numeric(K.modsel$unadj.values[,10])), Pa, xlab="log( Tau )", ylab="P(Acceptance into Local)")
pa.med <- (exp(-((mean(abs(dif)))^2)/median(as.numeric(params.M5[,10]))))
abline(h=pa.med, col="red", lty=2, lwd=1.5)


range(as.numeric(params.M4[,3]))
range(as.numeric(params.M4[,6]))
range(as.numeric(params.M4[,7]))
range(as.numeric(params.M4[,8]))
range(as.numeric(params.M4[,10]))

range(as.numeric(params.M5[,3]))
range(as.numeric(params.M5[,6]))
range(as.numeric(params.M5[,7]))
range(as.numeric(params.M5[,8]))
range(as.numeric(params.M5[,10]))

hist(as.numeric(params.M4[,3]), breaks=50)
abline(v=median(as.numeric(params.M4[,3])))
hist(as.numeric(params.M4[,7]), breaks=50)
abline(v=median(as.numeric(params.M4[,7])))


hist(as.numeric(params.M5[,2]), breaks=50)
abline(v=median(as.numeric(params.M5[,2])))


hist(as.numeric(params.M5[,5]), breaks=50)
abline(v=median(as.numeric(params.M5[,5])))

plot(density(as.numeric(params.M5[,5])), ylim=c(0,0.04))
points(density(tau.prior.mod5), type="l", lty=2)
tau.prior.mod5 <- sample(as.numeric(Com.Params[Com.Params[,1]=="mod2",5]), 400)

hist(as.numeric(params.M5[,3]), breaks=50)
abline(v=median(as.numeric(params.M5[,3])))
hist(as.numeric(params.M5[,6]), breaks=50) 
abline(v=median(as.numeric(params.M5[,6])), lwd=3)
hist(as.numeric(params.M5[,7]), breaks=50)
abline(v=median(as.numeric(params.M5[,7])))
hist(as.numeric(params.M5[,8]), breaks=50)

hist(as.numeric(params.M5[,10]))

hist(log(as.numeric(params.M5[,10])), breaks=10, main="", xlab="log( Tau )", add=T, col="blue")
abline(v=log(median((as.numeric(params.M5[,10])))), col="red", lwd=1.5, lty=2)

dif <- rep(NA, length(com.traits))
for (i in 1:length(dif)){
	dif[i] <- com.traits[i] - mean(com.traits)
}


Pa <- rep(NA, nrow(params.M5))
for (i in 1:length(Pa)){
	wght <- as.numeric(params.M5[i,5])
	Pa[i] <- (exp(-((mean(abs(dif)))^2)/wght))
}

plot(log(as.numeric(params.M5[,5])), Pa, xlab="log( Tau )", ylab="P(Acceptance into Local)")
pa.med <- (exp(-((mean(abs(dif)))^2)/median(as.numeric(params.M5[,10]))))
abline(h=pa.med, col="red", lty=2, lwd=1.5)
points(log(median(as.numeric(params.M5[,10]))), pa.med, bg="red", pch=21)
median(Pa)
as.numeric(rownames(K.modsel$ss))


plot(density(Pa))