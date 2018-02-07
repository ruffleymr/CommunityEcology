
###################################################################
#This script takes the second 10,000 Regional Species Phylogenies (simulated using the 'Simulate_RegionalSpeciesTrees.R' script) and simulates a trait under brownian motion onto the tips of the tree, and then assembles the community under habitat filtering processes 
###################################################################
#################################
#1.) Load needed packages
require(geiger)											#local machine
require(RPANDA)
require(phylobase)
require(phylosignal)
require(picante)
require(PhyloMeasures)
require(coda)
require(apTreeshape)
require(TotalCopheneticIndex)
require(nLTT)
require(phytools)
require(ape)


library(ape)			#ibest cluster
library(phylobase)
library(phytools)
library(geiger)	
library(picante)
library(PhyloMeasures)
library(phylosignal)
library(apTreeshape)
library(TotalCopheneticIndex)
library(nLTT)
library(RPANDA)


#will need this keep.tip function
keep.tip <- function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

#################################
#2.) Load Regional Trees, and pull out the first 10,000; trees[[1:10000]]

#setwd("/Users/Megan/Documents/ComPhy_SpecificScripts")
load(file="RegionalTrees.Rdata")
reg.trees <- trees[10001:20000]

#################################
#3.) Simulate Brownian Motion and Neutral Community; will need to draw rate parameter from a prior.
sims <- 10000			#Will need 10,000 rate parameters drawn from prior
mod <- rep("mod2", sims)		#Label each row as from model 1
sig2 <- runif(sims, 0.1, 5.0)	#rate parameters
tau.e <- runif(sims, 0, 1)
BMhf.params <- cbind(mod, sig2, tau.e, NA) ##Have two columns of 'NA' because other models with have parameters here, and all matrices need to have to 
colnames(BMhf.params) <- c("mod", "sig2", "tau.e", "alpha")
save(BMhf.params, file="BMhf_parameters.Rdata")

BMhf.RegTraits <- list()
BMhf.ComTraits <- list()
BMhf.ComTrees <- list()
Rej <- list()

for (i in 1:sims){
	#simulate traits on the regional phylogeny
	traits <- sim.char(reg.trees[[i]], par=sig2[i], nsim=1, model="BM", root=0)[,,1]
	BMhf.RegTraits[[i]] <- traits
	com.traits <- c()
	tr.ord <- sort(traits)
	wght <- tau.e[i] * (sqrt(sig2[i])*max(node.depth.edgelength(reg.trees[[i]])))
	
	if (runif(1,0,1) <= 0.5){
		Xi <- sample(tr.ord[1:as.integer(length(traits)/3)],1)
	}else{
		Xi <- sample(tr.ord[(length(traits)-as.integer(length(traits)/3)):length(traits)],1)
	}
	com.traits <- c(com.traits, Xi)
	traits <- traits[!traits==Xi]
	rej <- 0
	
	while (length(com.traits) < as.integer(length(traits)/2)) {
		Xj <- sample(traits,1)
		Probabilities <- rep(NA, length(com.traits))
		for (z in 1:length(com.traits)){
			Probabilities[z] <- (exp(-(((com.traits[z]-Xj)^2))/wght))
		}
		Pj <- mean(Probabilities)
		if (Pj > runif(1,0,1)) {
			com.traits <- c(com.traits, Xj)
			traits <- traits[!traits==Xj]
		}else{
			rej <- rej +1
		}
	}
	
	hf.tree <- keep.tip(reg.trees[[i]], names(com.traits))
	
	#store trait and phylogenetic information for the community
	BMhf.ComTraits[[i]] <- com.traits
	BMhf.ComTrees[[i]] <- hf.tree
	Rej[[i]] <- rej
}

save(BMhf.RegTraits, file="BMhf_RegTraits.Rdata")
save(BMhf.ComTraits, file="BMhf_traitData.Rdata")
save(BMhf.ComTrees, file="BMhf_ComTrees.Rdata")
save(Rej, file="BMhf_Rej.Rdata")
#################################
#4.) Calculate all summary statistics

SumStats <- matrix(NA, sims, 31)

for (k in 1:sims){
	
	spectR.stats <- spectR(BMhf.ComTrees[[k]])
	fourD <- phylo4d(BMhf.ComTrees[[k]], tip.data=BMhf.ComTraits[[k]])
	physig <- phyloSignal(fourD, reps=1, methods="all")$stat
	w <- 1/cophenetic(BMhf.ComTrees[[k]])
	diag(w) <- 0
	MI <-  Moran.I(BMhf.ComTraits[[k]], w)$observed
	
	mat <- matrix(NA, 1, length(reg.trees[[k]]$tip.label))
	colnames(mat) <- reg.trees[[k]]$tip.label
	for (j in 1:length(reg.trees[[k]]$tip.label)) {
		if (is.element(reg.trees[[k]]$tip.label[j], BMhf.ComTrees[[k]]$tip.label)){
			mat[j] <- 1
		} else {
			mat[j] <- 0 
		}
	}
	mat <- rbind(mat, mat)
	ses.mpd <- ses.mpd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	ses.mntd <- ses.mntd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	
	PIC <- pic(BMhf.ComTraits[[k]], BMhf.ComTrees[[k]], var.contrasts=T)
	Svar.lm <- lm(abs(PIC[,1]) ~ PIC[,2])
	
	tr.df <- data.frame(BMhf.ComTraits[[k]])
	tr.df <- tr.df[BMhf.ComTrees[[k]]$tip.label, ]
	tr.df <- data.frame(tr.df)
	rownames(tr.df) <- BMhf.ComTrees[[k]]$tip.label
	PICreconstruction <- ace(tr.df[,1], BMhf.ComTrees[[k]], type="continuous", method="ML")

	Sasr.lm <- lm(abs(PIC[,1]) ~ PICreconstruction$ace)
	Shgt.lm <- nh.test(BMhf.ComTrees[[k]], tr.df[,1], regression.type="lm", show.plot=FALSE)
	
	expCont.BM <- rnorm(n=length(BMhf.ComTraits[[k]]), mean=0, sd=sqrt(mean(PIC[,1]^2)))
	Dcdf <- ks.test(PIC[,1], expCont.BM)
	
	SumStats[k,] <-  c(spectR.stats$principal_eigenvalue,
						spectR.stats$asymmetry,
						spectR.stats$peakedness1,
						spectR.stats$peakedness2,
						spectR.stats$eigengap,
						mean(BMhf.ComTrees[[k]]$edge.length),
						var(BMhf.ComTrees[[k]]$edge.length),
						mean(BMhf.ComTraits[[k]]),
						var(BMhf.ComTraits[[k]]),
						physig[,1],
						physig[,2],
						physig[,3],
						physig[,4],
						MI,
						max(node.depth.edgelength(BMhf.ComTrees[[k]])),
						colless(as.treeshape(BMhf.ComTrees[[k]])),
						sackin(as.treeshape(BMhf.ComTrees[[k]])),
						tci(BMhf.ComTrees[[k]]),
						nLTTstat_exact(reg.trees[[k]],BMhf.ComTrees[[k]]),
						ses.mpd$mpd.obs,
						ses.mntd$mntd.obs,
						ses.mpd$mpd.obs.z, 
						ses.mntd$mntd.obs.z,
						ses.mpd$mpd.obs.p,
						ses.mntd$mntd.obs.p,
						mean(PIC[,1]^2),
						sd(PIC[,1])/mean(PIC[,1]),
						Svar.lm$coefficients[2],
						Sasr.lm$coefficients[2],
						Shgt.lm$coefficients[2],
						Dcdf$statistic
						)

}


colnames(SumStats) <-   c("Pr.Eigen",
							"Assymetry",
							"Peak1",
							"Peak2",
							"Eigen.Gap",
							"Mean.BL", 
							"Var.BL",
							"Mean.Tr",
							"Var.Tr",
							"Cmean",
							"I",
							"Blom.K",
							"K*",
							"Moran.I",
							"Age",
							"Colless",
							"Sackin",
							"tci",
							"nLTT",
							"MPD",
							"MNTD",
							"MPD.z",
							"MNTD.z",
							"MPD.p",
							"MNTD.p",
							"Msig",
							"Cvar",
							"Svar",
							"Sasr",
							"Shgt",
							"Dcdf"
							)

save(SumStats, file="BMhf_SumStats.Rdata")
