###################################################################
#This script takes the fourth 10,000 Regional Species Phylogenies (simulated using the 'Simulate_RegionalSpeciesTrees.R' script) and simulates a trait under brownian motion onto the tips of the tree, and then assembles the community under a neutral process 
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
reg.trees <- trees[30001:40000]

#################################
#3.) Simulate Brownian Motion and Neutral Community; will need to draw rate parameter from a prior.
sims <- 10000				#Will need 10,000 rate parameters drawn from prior
mod <- rep("mod4", sims)		#Label each row as from model 1
sig2 <- runif(sims, 0.01, 4.0)	#rate parameters
alpha <- runif(sims, 0.01, 0.1)
OUneut.params <- cbind(mod, sig2, NA, alpha) ##Have two columns of 'NA' because other models with have parameters here, and all matrices need to have to 
colnames(OUneut.params) <- c("mod", "sig2", "tau", "alpha")
save(OUneut.params, file="OUneut_parameters.Rdata")

OUneut.traits <- list()
OUneut.ComTrees <- list()
for (i in 1:sims){
	#simulate traits on the regional phylogeny
	traits <- fastBM(reg.trees[[i]], a=0, theta=0, alpha=alpha[i], sig2=sig2[i], internal=FALSE)
	
	#Because this is the neutral model, sample the traits randomly to build a community that is roughly half the size of the regional
	n <- length(reg.trees[[i]]$tip.label)
	com.traits <- sample(traits, as.integer(n/2))
	neut.tree <- keep.tip(reg.trees[[i]], names(com.traits))
	
	#store trait and phylogenetic information for the community
	OUneut.traits[[i]] <- com.traits
	OUneut.ComTrees[[i]] <- neut.tree
}

save(OUneut.traits, file="OUneut_traitData.Rdata")
save(OUneut.ComTrees, file="OUneut_ComTrees.Rdata")

#load(file="OUneut_parameters.Rdata")
#load(file="OUneut_traitData.Rdata")
#load(file="OUneut_ComTrees.Rdata")
#load(file="OUneut_SumStats.Rdata")

#################################
#4.) Calculate all summary statistics

SumStats <- matrix(NA, sims, 31)

for (k in 1:sims){
	
	spectR.stats <- spectR(OUneut.ComTrees[[k]])
	fourD <- phylo4d(OUneut.ComTrees[[k]], tip.data=OUneut.traits[[k]])
	physig <- phyloSignal(fourD, reps=1, methods="all")$stat
	w <- 1/cophenetic(OUneut.ComTrees[[k]])
	diag(w) <- 0
	MI <-  Moran.I(OUneut.traits[[k]], w)$observed
	
	mat <- matrix(NA, 1, length(reg.trees[[k]]$tip.label))
	colnames(mat) <- reg.trees[[k]]$tip.label
	for (j in 1:length(reg.trees[[k]]$tip.label)) {
		if (is.element(reg.trees[[k]]$tip.label[j], OUneut.ComTrees[[k]]$tip.label)){
			mat[j] <- 1
		} else {
			mat[j] <- 0 
		}
	}
	mat <- rbind(mat, mat)
	ses.mpd <- ses.mpd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	ses.mntd <- ses.mntd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	
	PIC <- pic(OUneut.traits[[k]], OUneut.ComTrees[[k]], var.contrasts=T)
	Svar.lm <- lm(abs(PIC[,1]) ~ PIC[,2])
	
	tr.df <- data.frame(OUneut.traits[[k]])
	tr.df <- tr.df[OUneut.ComTrees[[k]]$tip.label, ]
	tr.df <- data.frame(tr.df)
	rownames(tr.df) <- OUneut.ComTrees[[k]]$tip.label
	PICreconstruction <- ace(tr.df[,1], OUneut.ComTrees[[k]], type="continuous", method="ML")

	Sasr.lm <- lm(abs(PIC[,1]) ~ PICreconstruction$ace)
	Shgt.lm <- nh.test(OUneut.ComTrees[[k]], tr.df[,1], regression.type="lm", show.plot=FALSE)
	
	expCont.OU <- rnorm(n=length(OUneut.traits[[k]]), mean=0, sd=sqrt(mean(PIC[,1]^2)))
	Dcdf <- ks.test(PIC[,1], expCont.OU)
	
	SumStats[k,] <-  c(spectR.stats$principal_eigenvalue,
						spectR.stats$asymmetry,
						spectR.stats$peakedness1,
						spectR.stats$peakedness2,
						spectR.stats$eigengap,
						mean(OUneut.ComTrees[[k]]$edge.length),
						var(OUneut.ComTrees[[k]]$edge.length),
						mean(OUneut.traits[[k]]),
						var(OUneut.traits[[k]]),
						physig[,1],
						physig[,2],
						physig[,3],
						physig[,4],
						MI,
						max(node.depth.edgelength(OUneut.ComTrees[[k]])),
						colless(as.treeshape(OUneut.ComTrees[[k]])),
						sackin(as.treeshape(OUneut.ComTrees[[k]])),
						tci(OUneut.ComTrees[[k]]),
						nLTTstat_exact(reg.trees[[k]],OUneut.ComTrees[[k]]),
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
							
							
save(SumStats, file="OUneut_SumStats.Rdata")

