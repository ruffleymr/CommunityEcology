
###################################################################
#This script takes the last 10,000 Regional Species Phylogenies (simulated using the 'Simulate_RegionalSpeciesTrees.R' script) and simulates a trait under brownian motion onto the tips of the tree, and then assembles the community under competetive exclusion processes 
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

setwd("/Users/Megan/Documents/ComPhy_SpecificScripts")
load(file="RegionalTrees.Rdata")
reg.trees <- trees[50001:60000]

#################################
#3.) Simulate Brownian Motion and Neutral Community; will need to draw rate parameter from a prior.
sims <- 10000			#Will need 10,000 rate parameters drawn from prior
mod <- rep("mod6", sims)		#Label each row as from model 1
sig2 <- runif(sims, 0.1, 5.0)	#rate parameters
tau.c <- runif(sims, 0, 1)
alpha <- runif(sims, 0.01, 1)
OUce.params <- cbind(mod, sig2, tau.c, alpha) ##Have two columns of 'NA' because other models with have parameters here, and all matrices need to have to 
colnames(OUce.params) <- c("mod", "sig2", "tau.e", "alpha")
save(OUce.params, file="OUce_parameters.Rdata")

OUce.RegTraits <- list()
OUce.ComTraits <- list()
OUce.ComTrees <- list()
Rej <- list()

for (i in 1:sims){
	#simulate traits on the regional phylogeny
	traits <- fastBM(reg.trees[[i]], a=0, theta=0, alpha=alpha[i], sig2=sig2[i], internal=FALSE)
	OUce.RegTraits[[i]] <- traits
	com.traits <- c()
	wght <- tau.c[i] * ((sqrt(sig2[i])*max(node.depth.edgelength(reg.trees[[i]]))))
	Xi <- sample(traits, 1)
	com.traits <- c(com.traits, Xi)
	traits <- traits[!traits==Xi]
	rej <- 0
	
	while (length(com.traits) < as.integer(length(traits)/2)) {
		Xj <- sample(traits,1)
		Probabilities <- rep(NA, length(com.traits))
		for (z in 1:length(com.traits)){
			Probabilities[z] <- (1 - exp(-(((com.traits[z]-Xj)^2)/wght)))
		}
		Pj <- mean(Probabilities)
		if (Pj > runif(1,0,1)) {
			com.traits <- c(com.traits, Xj)
			traits <- traits[!traits==Xj]
		}else{
			rej <- rej +1
		}
	}
	
	ce.tree <- keep.tip(reg.trees[[i]], names(com.traits))
	
	#store trait and phylogenetic information for the community
	OUce.ComTraits[[i]] <- com.traits
	OUce.ComTrees[[i]] <- ce.tree
	Rej[[i]] <- rej
}

save(OUce.RegTraits, file="OUce_RegTraits.Rdata")
save(OUce.ComTraits, file="OUce_traitData.Rdata")
save(OUce.ComTrees, file="OUce_ComTrees.Rdata")
save(Rej, file="OUce_Rej.Rdata")
#################################
#4.) Calculate all summary statistics

SumStats <- matrix(NA, sims, 31)

for (k in 1:sims){
	
	spectR.stats <- spectR(OUce.ComTrees[[k]])
	fourD <- phylo4d(OUce.ComTrees[[k]], tip.data=OUce.ComTraits[[k]])
	physig <- phyloSignal(fourD, reps=1, methods="all")$stat
	w <- 1/cophenetic(OUce.ComTrees[[k]])
	diag(w) <- 0
	MI <-  Moran.I(OUce.ComTraits[[k]], w)$observed
	
	mat <- matrix(NA, 1, length(reg.trees[[k]]$tip.label))
	colnames(mat) <- reg.trees[[k]]$tip.label
	for (j in 1:length(reg.trees[[k]]$tip.label)) {
		if (is.element(reg.trees[[k]]$tip.label[j], OUce.ComTrees[[k]]$tip.label)){
			mat[j] <- 1
		} else {
			mat[j] <- 0 
		}
	}
	mat <- rbind(mat, mat)
	ses.mpd <- ses.mpd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	ses.mntd <- ses.mntd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	
	PIC <- pic(OUce.ComTraits[[k]], OUce.ComTrees[[k]], var.contrasts=T)
	Svar.lm <- lm(abs(PIC[,1]) ~ PIC[,2])
	
	tr.df <- data.frame(OUce.ComTraits[[k]])
	tr.df <- tr.df[OUce.ComTrees[[k]]$tip.label, ]
	tr.df <- data.frame(tr.df)
	rownames(tr.df) <- OUce.ComTrees[[k]]$tip.label
	PICreconstruction <- ace(tr.df[,1], OUce.ComTrees[[k]], type="continuous", method="ML")

	Sasr.lm <- lm(abs(PIC[,1]) ~ PICreconstruction$ace)
	Shgt.lm <- nh.test(OUce.ComTrees[[k]], tr.df[,1], regression.type="lm", show.plot=FALSE)
	
	expCont.OU <- rnorm(n=length(OUce.ComTraits[[k]]), mean=0, sd=sqrt(mean(PIC[,1]^2)))
	Dcdf <- ks.test(PIC[,1], expCont.OU)
	
	SumStats[k,] <-  c(spectR.stats$principal_eigenvalue,
						spectR.stats$asymmetry,
						spectR.stats$peakedness1,
						spectR.stats$peakedness2,
						spectR.stats$eigengap,
						mean(OUce.ComTrees[[k]]$edge.length),
						var(OUce.ComTrees[[k]]$edge.length),
						mean(OUce.ComTraits[[k]]),
						var(OUce.ComTraits[[k]]),
						physig[,1],
						physig[,2],
						physig[,3],
						physig[,4],
						MI,
						max(node.depth.edgelength(OUce.ComTrees[[k]])),
						colless(as.treeshape(OUce.ComTrees[[k]])),
						sackin(as.treeshape(OUce.ComTrees[[k]])),
						tci(OUce.ComTrees[[k]]),
						nLTTstat_exact(reg.trees[[k]],OUce.ComTrees[[k]]),
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

save(SumStats, file="OUce_SumStats.Rdata")
