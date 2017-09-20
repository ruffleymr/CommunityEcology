#################################################################################################################
##Generate Reference Table and Tables of Summary Statistics for cross validation of ABC
##Currently, simulating under 2 community assembly models; neutral and habitat similarity with 5 summary stats
#################################################################################################################
##1.) Load needed packages
library(abc)
library(phytools)
library(RPANDA)
library(TreeSim)
library(phytools)
library(geiger)
library(OUwie)
library(phylobase)
library(phylosignal)
library(picante)
library(PhyloMeasures)

##loading packages on cluster at UI
library(abc)
library(phytools, lib="/mnt/lfs2/ruff6699/Rpackages")
library(RPANDA, lib="/mnt/lfs2/ruff6699/Rpackages")
library(TreeSim, lib="/mnt/lfs2/ruff6699/Rpackages")
library(geiger, lib="/mnt/lfs2/ruff6699/Rpackages")
#library(OUwie, lib="/mnt/lfs2/ruff6699/Rpackages")
library(phylobase, lib="/mnt/lfs2/ruff6699/Rpackages")
library(phylosignal, lib="/mnt/lfs2/ruff6699/Rpackages")
library(picante, lib="/tmp/RtmpsHUpxd/downloaded_packages")
library(PhyloMeasures)

#################################################################################################################
##Simulate Birth Death Trees; draw parameters for speciation rate, extinction rate, mean of random deviates, and variance of random deviates, and total time.
##on local machine
setwd("/Users/Megan/Documents/CommunityPhylogenetics/")
##on cluster
setwd("/mnt/lfs2/ruff6699/ComPhy")
##Variables that will be held constant throughout the simulations
#Create Matrices to hold all parameters simulated under, and the summary stats from those simulations

keep.tip <- function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))
sims <- 1000

#sink(file=paste(paste(getwd(),"/GenRefTable_OU_BM_ABC_4262017.Rout",sep="")))
par.sim <- matrix(NA, sims, 9)
stat.sim <- matrix(NA, sims, 15)
trait.mod <- matrix(NA, sims, 5)
com.mod	<- matrix(NA, sims, 7)


for (i in 1:nrow(par.sim)){
	u <- floor(runif(1,1, 7))
	n <- as.integer(runif(1,60,100))
	numbsim<-1
	lambda <- runif(1,0.5, 8.0)
	eps <- runif(1, 0.1, 0.9) ##eps is death over birth
	mu <- eps*lambda
	frac <- runif(1,0.3,0.8)
	age <- as.integer(runif(1,80,120))
	par <- runif(1,0,2.0)
	sig2 <- runif(1,0,2.0)
	alpha <- runif(1,0.1,4.0)
	tree <- sim.bd.taxa.age(n, numbsim, lambda, mu, frac=frac, age, mrca = FALSE)[[1]]
	phydist <-  cophenetic(tree)
	vec <- tree$tip.label
	mat <- matrix(NA, 2, length(tree$tip.label))
	colnames(mat) <- tree$tip.label
	rownames(mat) <- c("com1", "com2")
	
	######################################################################################################################
	if (u == 1){ ##neutral model BM
		#determine what all simulate paramters will be
		dtips <- sample(tree$tip.label, as.integer(n/2))
		neut.tree <- drop.tip(tree, dtips)
		
		##Simulate Traits
		x <- sim.char(neut.tree, par, nsim=1,model="BM", root=0)[,,1]
		sts <- spectR(neut.tree)	
		fourD <- phylo4d(neut.tree, tip.data=x)
		phsig <- phyloSignal(fourD, reps=1, methods="all")$stat
		w <- 1/cophenetic(neut.tree)
		diag(w) <- 0
		MI <-  Moran.I(x, w)$observed
		
		##Record all sum stats and parameters
		stat.sim[i,] <- c(sts$principal_eigenvalue, sts$asymmetry, sts$peakedness1, sts$peakedness2, sts$eigengap, mean(neut.tree$edge.length), var(neut.tree$edge.length), mean(x), var(x), phsig[,1], phsig[,2], phsig[,3], phsig[,4], phsig[,5], MI)
		par.sim[i,] <- c(u, n, lambda, eps, mu, frac, age, par, "NA")
		
		###########################################################
		##Determine best fitting model of trait evoltion and record
		bm.fit <- fitContinuous(neut.tree, dat=x, model="BM")
		ou.fit <- fitContinuous(neut.tree, dat=x, model="OU")
		mod = ""
		dif <- abs(ou.fit$opt$aicc - bm.fit$opt$aicc)
		if (bm.fit$opt$aicc < ou.fit$opt$aicc){
			mod = "BM"
		}else{
			mod= "OU"
		}
		trait.mod[i,] <- c("BM", mod, bm.fit$opt$aicc, ou.fit$opt$aicc, dif)
		
		##Determine MPD and MNTD and record p-values compared to null
		com.vec <- neut.tree$tip.label
		for (b in seq(1,length(tree$tip.label))) {
			if (is.element(vec[b], com.vec)) {
				mat[1,b] <- 1
				mat[2,b] <- 1
			}else{
				mat[1,b] <- 0
				mat[2,b] <- 0
			}
		}
		ses.mpd <- ses.mpd(mat, phydist, null.model="taxa.labels")[1,]
		ses.mntd <- ses.mntd(mat, phydist, null.model="taxa.labels")[1,]
		com.mod[i,] <- c("neut", mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z,  ses.mpd$mpd.obs.p, ses.mntd$mntd.obs.p)
		
	}
	######################################################################################################################
	if (u == 2){ ##habitat filtering BM
		x <- sim.char(tree, par, nsim=1,model="BM", root=0)[,,1] 
		x.order <- sort(x)
		#p <- length(tree$tip.label)
		dtips <- names(x.order[as.integer((n/2)+1):n])
		xtips <- x.order[1:as.integer(n/2)]
		habitat.tree <- drop.tip(tree, dtips)
		sts <- spectR(habitat.tree)	
		fourD <- phylo4d(habitat.tree, tip.data=xtips)
		phsig <- phyloSignal(fourD, reps=1, methods="all")$stat
		w <- 1/cophenetic(habitat.tree)
		diag(w) <- 0
		MI <-  Moran.I(xtips, w)$observed
		stat.sim[i,1:15] <- c(sts$principal_eigenvalue, sts$asymmetry, sts$peakedness1, sts$peakedness2, sts$eigengap, mean(habitat.tree$edge.length), var(habitat.tree$edge.length), mean(xtips), var(xtips), phsig[,1], phsig[,2], phsig[,3], phsig[,4], phsig[,5], MI)
		par.sim[i,] <- c(u, n, lambda, eps, mu, frac, age, par, "NA")
		
		###########################################################
		##Determine best fitting model of trait evoltion and record
		bm.fit <- fitContinuous(habitat.tree, dat=xtips, model="BM")
		ou.fit <- fitContinuous(habitat.tree, dat=xtips, model="OU")
		mod = ""
		dif <- abs(ou.fit$opt$aicc - bm.fit$opt$aicc)
		if (bm.fit$opt$aicc < ou.fit$opt$aicc){
			mod = "BM"
		}else{
			mod= "OU"
		}
		trait.mod[i,] <- c("BM", mod, bm.fit$opt$aicc, ou.fit$opt$aicc, dif)
		
		##Determine MPD and MNTD and record p-values compared to null
		com.vec <- habitat.tree$tip.label
		for (b in seq(1,length(tree$tip.label))) {
			if (is.element(vec[b], com.vec)) {
				mat[1,b] <- 1
				mat[2,b] <- 1
			}else{
				mat[1,b] <- 0
				mat[2,b] <- 0
			}
		}
		ses.mpd <- ses.mpd(mat, phydist, null.model="taxa.labels")[1,]
		ses.mntd <- ses.mntd(mat, phydist, null.model="taxa.labels")[1,]
		com.mod[i,] <- c("habitat", mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z,  ses.mpd$mpd.obs.p, ses.mntd$mntd.obs.p)
		
		
	}
	######################################################################################################################
	if (u == 3){ ##limited similarity
		x <- sim.char(tree, par, nsim=1,model="BM", root=0)[,,1] 
		limx <- c()
		prop <- sample(x,1)
		limx <- c(prop)
		while (length(limx) < length(unique(round(x, digit=1))) && length(limx) < n/2) {
			prop <- sample(x,1)
			if (round(prop, digit=1) %in% round(limx, digit=1)){
				prop <- sample(x,1)
				}
			else{
				limx <- c(limx, prop)
				}
		}	
		limsim.tree <- keep.tip(tree, names(limx))
		xtips <- limx
		sts <- spectR(limsim.tree)	
		fourD <- phylo4d(limsim.tree, tip.data=xtips)
		phsig <- phyloSignal(fourD, reps=1, methods="all")$stat
		w <- 1/cophenetic(limsim.tree)
		diag(w) <- 0
		MI <-  Moran.I(xtips, w)$observed
		stat.sim[i,1:15] <- c(sts$principal_eigenvalue, sts$asymmetry, sts$peakedness1, sts$peakedness2, sts$eigengap, mean(limsim.tree$edge.length), var(limsim.tree$edge.length), mean(xtips), var(xtips), phsig[,1], phsig[,2], phsig[,3], phsig[,4], phsig[,5], MI)	
		par.sim[i,] <- c(u, n, lambda, eps, mu, frac, age, par, "NA")
		
		###########################################################
		##Determine best fitting model of trait evoltion and record
		bm.fit <- fitContinuous(limsim.tree, dat=xtips, model="BM")
		ou.fit <- fitContinuous(limsim.tree, dat=xtips, model="OU")
		mod = ""
		dif <- abs(ou.fit$opt$aicc - bm.fit$opt$aicc)
		if (bm.fit$opt$aicc < ou.fit$opt$aicc){
			mod = "BM"
		}else{
			mod= "OU"
		}
		trait.mod[i,] <- c("BM", mod, bm.fit$opt$aicc, ou.fit$opt$aicc, dif)
		
		##Determine MPD and MNTD and record p-values compared to null
		com.vec <- limsim.tree$tip.label
		for (b in seq(1,length(tree$tip.label))) {
			if (is.element(vec[b], com.vec)) {
				mat[1,b] <- 1
				mat[2,b] <- 1
			}else{
				mat[1,b] <- 0
				mat[2,b] <- 0
			}
		}
		ses.mpd <- ses.mpd(mat, phydist, null.model="taxa.labels")[1,]
		ses.mntd <- ses.mntd(mat, phydist, null.model="taxa.labels")[1,]
		com.mod[i,] <- c("limsim", mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z,  ses.mpd$mpd.obs.p, ses.mntd$mntd.obs.p)

		
	}
	######################################################################################################################
	if (u == 4){ ##neutral model OU
		#determine what all simulate paramters will be
		dtips <- sample(tree$tip.label, as.integer(n/2))
		neut.tree <- drop.tip(tree, dtips)
		x <- fastBM(neut.tree, a=0,theta=0, alpha=alpha, sig2=sig2, internal=FALSE)
		sts <- spectR(neut.tree)	
		fourD <- phylo4d(neut.tree, tip.data=x)
		phsig <- phyloSignal(fourD, reps=1, methods="all")$stat
		w <- 1/cophenetic(neut.tree)
		diag(w) <- 0
		MI <-  Moran.I(x, w)$observed
		stat.sim[i,1:15] <- c(sts$principal_eigenvalue, sts$asymmetry, sts$peakedness1, sts$peakedness2, sts$eigengap, mean(neut.tree$edge.length), var(neut.tree$edge.length), mean(x), var(x), phsig[,1], phsig[,2], phsig[,3], phsig[,4], phsig[,5], MI)
		par.sim[i,] <- c(u, n, lambda, eps, mu, frac, age, sig2, alpha)
		
		###########################################################
		##Determine best fitting model of trait evoltion and record
		bm.fit <- fitContinuous(neut.tree, dat=x, model="BM")
		ou.fit <- fitContinuous(neut.tree, dat=x, model="OU")
		mod = ""
		dif <- abs(ou.fit$opt$aicc - bm.fit$opt$aicc)
		if (bm.fit$opt$aicc < ou.fit$opt$aicc){
			mod = "BM"
		}else{
			mod= "OU"
		}
		trait.mod[i,] <- c("OU", mod, bm.fit$opt$aicc, ou.fit$opt$aicc, dif)
		
		##Determine MPD and MNTD and record p-values compared to null
		com.vec <- neut.tree$tip.label
		for (b in seq(1,length(tree$tip.label))) {
			if (is.element(vec[b], com.vec)) {
				mat[1,b] <- 1
				mat[2,b] <- 1
			}else{
				mat[1,b] <- 0
				mat[2,b] <- 0
			}
		}
		ses.mpd <- ses.mpd(mat, phydist, null.model="taxa.labels")[1,]
		ses.mntd <- ses.mntd(mat, phydist, null.model="taxa.labels")[1,]
		com.mod[i,] <- c("neut", mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z,  ses.mpd$mpd.obs.p, ses.mntd$mntd.obs.p)
		
		
	}
	######################################################################################################################	
	if (u == 5){ ##habitat filtering OU
		x <- fastBM(tree, a=0,theta=0, alpha=alpha, sig2=sig2, internal=FALSE)
		x.order <- sort(x)
		#p <- length(tree$tip.label)
		dtips <- names(x.order[as.integer((n/2)+1):n])
		xtips <- x.order[1:as.integer(n/2)]
		habitat.tree <- drop.tip(tree, dtips)
		sts <- spectR(habitat.tree)	
		fourD <- phylo4d(habitat.tree, tip.data=xtips)
		phsig <- phyloSignal(fourD, reps=1, methods="all")$stat
		w <- 1/cophenetic(habitat.tree)
		diag(w) <- 0
		MI <-  Moran.I(xtips, w)$observed
		stat.sim[i,1:15] <- c(sts$principal_eigenvalue, sts$asymmetry, sts$peakedness1, sts$peakedness2, sts$eigengap, mean(habitat.tree$edge.length), var(habitat.tree$edge.length), mean(xtips), var(xtips), phsig[,1], phsig[,2], phsig[,3], phsig[,4], phsig[,5], MI)
		par.sim[i,] <- c(u, n, lambda, eps, mu, frac, age, sig2, alpha)
		
		###########################################################
		##Determine best fitting model of trait evoltion and record
		bm.fit <- fitContinuous(habitat.tree, dat=xtips, model="BM")
		ou.fit <- fitContinuous(habitat.tree, dat=xtips, model="OU")
		mod = ""
		dif <- abs(ou.fit$opt$aicc - bm.fit$opt$aicc)
		if (bm.fit$opt$aicc < ou.fit$opt$aicc){
			mod = "BM"
		}else{
			mod= "OU"
		}
		trait.mod[i,] <- c("OU", mod, bm.fit$opt$aicc, ou.fit$opt$aicc, dif)
		
		##Determine MPD and MNTD and record p-values compared to null
		com.vec <- habitat.tree$tip.label
		for (b in seq(1,length(tree$tip.label))) {
			if (is.element(vec[b], com.vec)) {
				mat[1,b] <- 1
				mat[2,b] <- 1
			}else{
				mat[1,b] <- 0
				mat[2,b] <- 0
			}
		}
		ses.mpd <- ses.mpd(mat, phydist, null.model="taxa.labels")[1,]
		ses.mntd <- ses.mntd(mat, phydist, null.model="taxa.labels")[1,]
		com.mod[i,] <- c("habitat", mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z,  ses.mpd$mpd.obs.p, ses.mntd$mntd.obs.p)
		
		
	}
	######################################################################################################################
	if (u == 6){ ##limited similarity OU
		x <- fastBM(tree, a=0,theta=0, alpha=alpha, sig2=sig2, internal=FALSE)
		limx <- c()
		prop <- sample(x,1)
		limx <- c(prop)
		while (length(limx) < length(unique(round(x, digit=1))) && length(limx) < n/2) {
			prop <- sample(x,1)
			if (is.element(round(prop, digit=1), round(limx, digit=1))){
				prop <- sample(x,1)
				}
			else{
				limx <- c(limx, prop)
				}
		}	
		limsim.tree <- keep.tip(tree, names(limx))
		xtips <- limx
		sts <- spectR(limsim.tree)	
		fourD <- phylo4d(limsim.tree, tip.data=xtips)
		phsig <- phyloSignal(fourD, reps=1, methods="all")$stat
		w <- 1/cophenetic(limsim.tree)
		diag(w) <- 0
		MI <-  Moran.I(xtips, w)$observed
		stat.sim[i,1:15] <- c(sts$principal_eigenvalue, sts$asymmetry, sts$peakedness1, sts$peakedness2, sts$eigengap, mean(limsim.tree$edge.length), var(limsim.tree$edge.length), mean(xtips), var(xtips), phsig[,1], phsig[,2], phsig[,3], phsig[,4], phsig[,5], MI)	
		par.sim[i,] <- c(u, n, lambda, eps, mu, frac, age, sig2, alpha)
	
		###########################################################
		##Determine best fitting model of trait evoltion and record
		bm.fit <- fitContinuous(limsim.tree, dat=xtips, model="BM")
		ou.fit <- fitContinuous(limsim.tree, dat=xtips, model="OU")
		mod = ""
		dif <- abs(ou.fit$opt$aicc - bm.fit$opt$aicc)
		if (bm.fit$opt$aicc < ou.fit$opt$aicc){
			mod = "BM"
		}else{
			mod= "OU"
		}
		trait.mod[i,] <- c("OU", mod, bm.fit$opt$aicc, ou.fit$opt$aicc, dif)
		
		##Determine MPD and MNTD and record p-values compared to null
		com.vec <- limsim.tree$tip.label
		for (b in seq(1,length(tree$tip.label))) {
			if (is.element(vec[b], com.vec)) {
				mat[1,b] <- 1
				mat[2,b] <- 1
			}else{
				mat[1,b] <- 0
				mat[2,b] <- 0
			}
		}
		ses.mpd <- ses.mpd(mat, phydist, null.model="taxa.labels")[1,]
		ses.mntd <- ses.mntd(mat, phydist, null.model="taxa.labels")[1,]
		com.mod[i,] <- c("limsim", mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z,  ses.mpd$mpd.obs.p, ses.mntd$mntd.obs.p)

	}
	
}
