###################################################################
#This script takes the first 10,000 Regional Species Phylogenies (simulated using the 'Simulate_RegionalSpeciesTrees.R' script) and simulates a trait under brownian motion onto the tips of the tree, and then assembles the community under a neutral process 
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
loading(phytools)
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
load(file="Large_RegionalTrees.Rdata")
reg.trees <- trees[1:10000]

#################################
#3.) Simulate Brownian Motion and Neutral Community; will need to draw rate parameter from a prior.
sims <- 10000				#Will need 10,000 rate parameters drawn from prior
mod <- rep("mod1", sims)		#Label each row as from model 1
sig2 <- runif(sims, 0.01, 4.0)	#rate parameters
BMneut.params <- cbind(mod, sig2, NA, NA) ##Have two columns of 'NA' because other models with have parameters here, and all matrices need to have to 
colnames(BMneut.params) <- c("mod", "sig2", "tau", "alpha")
#save(BMneut.params, file="BMneut_parameters.Rdata")

BMneut.traits <- list()
BMneut.ComTrees <- list()
for (i in 1:sims){
	#simulate traits on the regional phylogeny
	traits <- sim.char(reg.trees[[i]], par=sig2[i], nsim=1, model="BM", root=0)[,,1]
	
	#Because this is the neutral model, sample the traits randomly to build a community that is roughly half the size of the regional
	n <- length(reg.trees[[i]]$tip.label)
	com.traits <- sample(traits, as.integer(n/2))
	neut.tree <- keep.tip(reg.trees[[i]], names(com.traits))
	
	#store trait and phylogenetic information for the community
	BMneut.traits[[i]] <- com.traits
	BMneut.ComTrees[[i]] <- neut.tree
}

#save(BMneut.traits, file="BMneut_traitData.Rdata")
#save(BMneut.ComTrees, file="BMneut_ComTrees.Rdata")

load(file="BMneut_parameters.Rdata")
load(file="BMneut_traitData.Rdata")
load(file="BMneut_ComTrees.Rdata")
load(file="BMneut_SumStats.Rdata")

#################################
#4.) Calculate all summary statistics

#SumStats <- matrix(NA, sims, 31)

SumStats[3730,]
start <- 3730
for (k in start:sims){
	
	spectR.stats <- spectR(BMneut.ComTrees[[k]])
	fourD <- phylo4d(BMneut.ComTrees[[k]], tip.data=BMneut.traits[[k]])
	physig <- phyloSignal(fourD, reps=1, methods="all")$stat
	w <- 1/cophenetic(BMneut.ComTrees[[k]])
	diag(w) <- 0
	MI <-  Moran.I(BMneut.traits[[k]], w)$observed
	
	mat <- matrix(NA, 1, length(reg.trees[[k]]$tip.label))
	colnames(mat) <- reg.trees[[k]]$tip.label
	for (j in 1:length(reg.trees[[k]]$tip.label)) {
		if (is.element(reg.trees[[k]]$tip.label[j], BMneut.ComTrees[[k]]$tip.label)){
			mat[j] <- 1
		} else {
			mat[j] <- 0 
		}
	}
	mat <- rbind(mat, mat)
	ses.mpd <- ses.mpd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	ses.mntd <- ses.mntd(mat, cophenetic(reg.trees[[k]]), null.model="taxa.labels")[1,]
	
	PIC <- pic(BMneut.traits[[k]], BMneut.ComTrees[[k]], var.contrasts=T)
	Svar.lm <- lm(abs(PIC[,1]) ~ PIC[,2])
	
	tr.df <- data.frame(BMneut.traits[[k]])
	tr.df <- tr.df[BMneut.ComTrees[[k]]$tip.label, ]
	tr.df <- data.frame(tr.df)
	rownames(tr.df) <- BMneut.ComTrees[[k]]$tip.label
	PICreconstruction <- ace(tr.df[,1], BMneut.ComTrees[[k]], type="continuous", method="ML")

	Sasr.lm <- lm(abs(PIC[,1]) ~ PICreconstruction$ace)
	Shgt.lm <- nh.test(BMneut.ComTrees[[k]], tr.df[,1], regression.type="lm", show.plot=FALSE)
	
	expCont.BM <- rnorm(n=length(BMneut.traits[[k]]), mean=0, sd=sqrt(mean(PIC[,1]^2)))
	Dcdf <- ks.test(PIC[,1], expCont.BM)
	
	SumStats[k,] <-  c(spectR.stats$principal_eigenvalue,
						spectR.stats$asymmetry,
						spectR.stats$peakedness1,
						spectR.stats$peakedness2,
						spectR.stats$eigengap,
						mean(BMneut.ComTrees[[k]]$edge.length),
						var(BMneut.ComTrees[[k]]$edge.length),
						mean(BMneut.traits[[k]]),
						var(BMneut.traits[[k]]),
						physig[,1],
						physig[,2],
						physig[,3],
						physig[,4],
						MI,
						max(node.depth.edgelength(BMneut.ComTrees[[k]])),
						colless(as.treeshape(BMneut.ComTrees[[k]])),
						sackin(as.treeshape(BMneut.ComTrees[[k]])),
						tci(BMneut.ComTrees[[k]]),
						nLTTstat_exact(reg.trees[[k]],BMneut.ComTrees[[k]]),
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

}

save(SumStats, file="BMneut_SumStats.Rdata")




