library(ape)			#ibest cluster
library(phytools)
library(geiger)	
library(psych)
library(RPANDA)
library(phylobase)
library(phylosignal)
library(picante)
library(PhyloMeasures)
library(coda)
library(apTreeshape)
library(TotalCopheneticIndex)
library(nLTT)
library(modes)
library(randomForest)
library(abc)

load(file="RegTrees3.Rdata")		#reg.trees.1
load(file="ComTrees.Rdata")
Com.tre1 <- Com.Trees
load(file="ComTrees2.Rdata")
Com.Trees <- rbind(Com.tre1, Com.Trees)



load(file="RegTrees.Rdata")		#reg.trees.1
load(file="RegPars2.Rdata")			#reg.pars
load(file="RegTraits.Rdata")		#Reg.Traits
load(file="ComTrees.Rdata")			#Com.Trees
load(file="ComParams.Rdata")		#Com.Params
load(file="ComTraits.Rdata")		#Com.Traits

load(file="ComTrees.Rdata")			#Com.Trees
ct <- Com.Trees
load(file="ComTrees2.Rdata")	
Com.Trees <- c(ct, Com.Trees)

sims=12000

SumStats <- matrix(NA, sims, 31)

for (i in 1:sims){ 
	SumStats[i,1] <- mean(Com.Trees[[i]]$edge.length)
	SumStats[i,2] <- var(Com.Trees[[i]]$edge.length)
	SumStats[i,3] <- mean(Com.Traits[[i]])
	SumStats[i,4] <- var(Com.Traits[[i]])
	print(i)
}

for (i in 1:sims){ 
	w <- 1/cophenetic(Com.Trees[[i]])
	diag(w) <- 0
	SumStats[i,5] <-  Moran.I(Com.Traits[[i]], w)$observed
	print(i)
}

for (i in 1:sims){
	SumStats[i,6] <- max(node.depth.edgelength(Com.Trees[[i]]))								
	SumStats[i,7] <- colless(as.treeshape(Com.Trees[[i]]))
	SumStats[i,8] <- sackin(as.treeshape(Com.Trees[[i]]))
	#SumStats[i,9] <- tci(Com.Trees[[i]])
	SumStats[i,10] <- nLTTstat_exact(reg.trees[[i]], Com.Trees[[i]])
	print(i)
}

##running on TCI screen
load(file="ComTrees.Rdata")
tci.list <- c()
for (i in 1:sims){
	tci.list <- c(tci.list, tci(Com.Trees[[i]]))
	print(i)
}
load(file="TCI.Rdata")

##running on PenStats screen

Pen.Stats <- matrix(NA, sims, 5)

for (i in 1:sims){
	PIC <- pic(Com.Traits[[i]], Com.Trees[[i]], var.contrasts=T)
	Pen.Stats[i,1] <- mean(PIC[,1]^2)
	Pen.Stats[i,2] <- sd(PIC[,1])/mean(PIC[,1])
	Svar.lm <- lm(abs(PIC[,1]) ~ PIC[,2])
	Pen.Stats[i,3] <- Svar.lm$coefficients[2]

	tr.df <- data.frame(Com.Traits[[i]])
	tr.df <- tr.df[Com.Trees[[i]]$tip.label, ]
	tr.df <- data.frame(tr.df)
	rownames(tr.df) <- Com.Trees[[i]]$tip.label
	
	PIC <- pic(Com.Traits[[i]], Com.Trees[[i]], var.contrasts=T)
	#PICreconstruction <- ace(tr.df[,1], Com.Trees[[i]], type="continuous", method="ML", CI=F)
	#Sasr.lm <- lm(abs(PIC[,1]) ~ PICreconstruction$ace)
	
	Shgt.lm <- nh.test(Com.Trees[[i]], tr.df[,1], regression.type="lm", show.plot=FALSE)
	
	expCont.BM <- rnorm(n=length(Com.Traits[[i]]), mean=0, sd=sqrt(mean(PIC[,1]^2)))
	Dcdf <- ks.test(PIC[,1], expCont.BM)

	#Pen.Stats[i,4] <- Sasr.lm$coefficients[2]
	Pen.Stats[i,4] <- Shgt.lm$coefficients[2]
	Pen.Stats[i,5] <- Dcdf$statistic
	print(i)
}

save(Pen.Stats, file="PenStats.Rdata")
load(file="PenStats.Rdata")


for (i in 1:sims){
	SumStats[i,15] <- skewness(Com.Traits[[i]])
	SumStats[i,16] <- kurtosis(Com.Traits[[i]], finite=T)
	SumStats[i,17] <- mean(Reg.Traits[[i]])
	SumStats[i,18] <- var(Reg.Traits[[i]])
	SumStats[i,19] <- mean(Reg.Traits[[i]]) - mean(Com.Traits[[i]])
	SumStats[i,20] <- var(Reg.Traits[[i]]) - var(Com.Traits[[i]])
	print(i)
	
}	

for (i in 1:sims){
	SumStats[i,21] <- mean(reg.trees[[i]]$edge.length)
	SumStats[i,22] <- var(reg.trees[[i]]$edge.length)
	SumStats[i,23] <- mean(reg.trees[[i]]$edge.length) - mean(Com.Trees[[i]]$edge.length)
	SumStats[i,24] <- var(reg.trees[[i]]$edge.length) - var(Com.Trees[[i]]$edge.length)
	print(i)
}

for (i in 1:sims){						
	SumStats[i,25] <- amps(Com.Traits[[i]])$Peaks[1]
	SumStats[i,26] <- amps(Com.Traits[[i]])$Peaks[2]
	SumStats[i,27] <- bimodality_coefficient(Com.Traits[[i]])
	SumStats[i,28] <- bimodality_ratio(Com.Traits[[i]])
	if (is.na(SumStats[i,28])){
		SumStats[i,28] <- 0
	}
	SumStats[i,29] <- modes(Com.Traits[[i]])[1]
	SumStats[i,30] <- modes(Com.Traits[[i]])[2]
	print(i)
}


load(file="SumStats.Rdata")
head(SumStats)
ncol(SumStats)

############################################
spectR.stats <- matrix(NA, sims, 5)

for (i in 1:sims){
	p <- spectR(Com.Trees[[i]])
	spectR.stats[i,] <- c(p$principal_eigenvalue,
						p$asymmetry,
						p$peakedness1,
						p$peakedness2,
						p$eigengap)
	print(i)
}
save(spectR.stats, file="spectRstats.Rdata")

load(file="spectRstats.Rdata")

############################################
PhySignal <- matrix(NA, sims, 4)

for (i in 1:sims){
	fourD <- phylo4d(Com.Trees[[i]], tip.data=Com.Traits[[i]])
	physig <- phyloSignal(fourD, reps=1, methods="all")$stat
	PhySignal[i,] <- c(physig[,1],
						physig[,2],
						physig[,3],
						physig[,4])
	print(i)
}
save(PhySignal, file="PhySignal.Rdata")

load(file="PhySignal.Rdata")

############################################
PhyDisp <- matrix(NA, sims, 6)

for (i in 2001:sims){
	mat <- matrix(NA, 1, length(reg.trees[[i]]$tip.label))
	colnames(mat) <- reg.trees[[i]]$tip.label
	for (j in 1:length(reg.trees[[i]]$tip.label)) {
		if (is.element(reg.trees[[i]]$tip.label[j], Com.Trees[[i]]$tip.label)){
			mat[j] <- 1
		} else {
			mat[j] <- 0 
		}
	}
	mat <- rbind(mat, mat)
	ses.mpd <- ses.mpd(mat, cophenetic(reg.trees[[i]]), runs=50, iterations=1, null.model="taxa.labels")[1,]
	ses.mntd <- ses.mntd(mat, cophenetic(reg.trees[[i]]), runs=50, iterations=1, null.model="taxa.labels")[1,]
	
	PhyDisp[i,] <- c(ses.mpd$mpd.obs, ses.mntd$mntd.obs, ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z, 
					ses.mpd$mpd.obs.p,ses.mntd$mntd.obs.p)
	print(i)
}	
save(PhyDisp, file="PhyDisp.Rdata")

load(file="PhyDisp.Rdata")
				

nrow(SumStats)

length(tci.list)
nrow(spectR.stats)
nrow(PhySignal)
nrow(PhyDisp)

SumStats <- cbind(SumStats, tci.list, spectR.stats, PhySignal, PhyDisp)
					
############################################
Sasr <- c()						
for (i in 1:sims){
	PIC <- pic(Com.Traits[[i]], Com.Trees[[i]], var.contrasts=T)
	tr.df <- data.frame(Com.Traits[[i]])
	tr.df <- tr.df[Com.Trees[[i]]$tip.label, ]
	tr.df <- data.frame(tr.df)
	rownames(tr.df) <- Com.Trees[[i]]$tip.label
	PICreconstruction <- ace(tr.df[,1], Com.Trees[[i]], type="continuous", method="ML")
	Sasr.lm <- lm(abs(PIC[,1]) ~ PICreconstruction$ace)
	Sasr <- c(Sasr, Sasr.lm$coefficients[2])
	print(i)
}
save(Sasr, file="Sasr.Rdata")	
 