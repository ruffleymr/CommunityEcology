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
library(coda)
library(apTreeshape)
library(TotalCopheneticIndex)
library(nLTT)

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
library(apTreeshape, lib="/mnt/lfs2/ruff6699/Rpackages")
library(TotalCopheneticIndex, lib="/mnt/lfs2/ruff6699/Rpackages")
library(nLTT, lib="/mnt/lfs2/ruff6699/Rpackages")


#################################################################################################################
##Simulate Birth Death Trees; draw parameters for speciation rate, extinction rate, mean of random deviates, and variance of random deviates, and total time.
##on local machine
setwd("/Users/Megan/Documents/CommunityPhylogenetics/")
##on cluster
setwd("/mnt/lfs2/ruff6699/ComPhy_Sep21")
##Variables that will be held constant throughout the simulations
#Create Matrices to hold all parameters simulated under, and the summary stats from those simulations

keep.tip <- function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))
sims <- 10

#sink(file=paste(paste(getwd(),"/GenRefTable_OU_BM_ABC_4262017.Rout",sep="")))
par.sim <- matrix(NA, sims, 9)
stat.sim <- matrix(NA, sims, 24)
trait.mod <- matrix(NA, sims, 5)
com.mod	<- matrix(NA, sims, 7)


for (i in 1:nrow(par.sim)){
	u <- floor(runif(1,1, 7))
	n <- as.integer(runif(1,200,500))
	numbsim<-1
	lambda <- runif(1,0.5, 8.0)
	eps <- runif(1, 0.1, 0.9) ##eps is death over birth
	mu <- eps*lambda
	frac <- runif(1,0.3,0.8)
	age <- as.integer(runif(1,80,120))
	par <- runif(1,0,5.0)
	sig2 <- runif(1,0,5.0)
	tau <- sig2
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
		print(1)
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
		stat.sim[i,1:15] <- c(sts$principal_eigenvalue, sts$asymmetry, sts$peakedness1, sts$peakedness2, sts$eigengap, mean(neut.tree$edge.length), var(neut.tree$edge.length), mean(x), var(x), phsig[,1], phsig[,2], phsig[,3], phsig[,4], phsig[,5], MI)
		par.sim[i,] <- c(u, n, lambda, eps, mu, frac, age, par, "NA")
		
		###########################################################
		##Determine best fitting model of trait evoltion and record
		bm.fit <- fitContinuous(neut.tree, dat=x, model="BM", control=list(niter=10))
		ou.fit <- fitContinuous(neut.tree, dat=x, model="OU", control=list(niter=10), bounds=list(alpha = c(min = exp(-500), max = 5.0)))
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
		stat.sim[i,16:24] <- c(mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z, colless(as.treeshape(neut.tree)), sackin(as.treeshape(neut.tree)), tci(neut.tree), nLTTstat_exact(tree, neut.tree), nLTTstat_exact(tree, neut.tree, ignore_stem=F))

	}
	######################################################################################################################
	if (u == 2){ ##habitat filtering BM
		print(2)
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
		bm.fit <- fitContinuous(habitat.tree, dat=xtips, model="BM", control=list(niter=10))
		ou.fit <- fitContinuous(habitat.tree, dat=xtips, model="OU", control=list(niter=10), bounds=list(alpha = c(min = exp(-500), max = 5.0)))
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
		stat.sim[i,16:24] <- c(mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z, colless(as.treeshape(habitat.tree)), sackin(as.treeshape(habitat.tree)), tci(habitat.tree), nLTTstat_exact(tree, habitat.tree), nLTTstat_exact(tree, habitat.tree, ignore_stem=F))
		
	}
	######################################################################################################################
	if (u == 3){ ##limited similarity
		print(3)
		#par <- runif(1,0,5.0)
		x <- sim.char(tree, par, nsim=1,model="BM", root=0)[,,1] 
		#hist(x, col="blue", breaks=10)
		alph<- abs(range(x)[1]-range(x)[2])/2
		limx <- c()
		Xi <- sample(x,1)
		limx <- c(Xi)
		x <- x[!x==Xi]
		#acc<-0
		#rej<-0
		while (length(limx) < n/2) {
			Xj <- sample(x,1)
			Probabilities <- rep(NA, length(limx))
			for (z in 1:length(limx)){
				Probabilities[z] <- (1 - exp(-(((limx[z]-Xj)^2)/(alph*tau))))
			}
			Pj <- mean(Probabilities)
			if (Pj > runif(1,0,1)) {
				limx <- c(limx, Xj)
				x <- x[!x==Xj]
				#acc <- acc +1
			}else{
				#x <- x[!x==Xj]
				#rej <- rej +1
			}	
		}
		#par
		#alph
		#acc/sum(acc,rej)
		#rej/sum(acc,rej)
		#par(new=T)
		#hist(limx, col="gray", xlim=c(-8,4), ylim=c(0,60))
		
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
		bm.fit <- fitContinuous(limsim.tree, dat=xtips, model="BM", control=list(niter=10))
		ou.fit <- fitContinuous(limsim.tree, dat=xtips, model="OU", control=list(niter=10), bounds=list(alpha = c(min = exp(-500), max = 5.0)))
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
		stat.sim[i,16:24] <- c(mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z, colless(as.treeshape(limsim.tree)), sackin(as.treeshape(limsim.tree)), tci(limsim.tree), nLTTstat_exact(tree, limsim.tree), nLTTstat_exact(tree, limsim.tree, ignore_stem=F))
		
	}
	######################################################################################################################
	if (u == 4){ ##neutral model OU
		#determine what all simulate paramters will be
		print(4)
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
		bm.fit <- fitContinuous(neut.tree, dat=x, model="BM", control=list(niter=10))
		ou.fit <- fitContinuous(neut.tree, dat=x, model="OU", control=list(niter=10), bounds=list(alpha = c(min = exp(-500), max = 5.0)))
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
		stat.sim[i,16:24] <- c(mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z, colless(as.treeshape(neut.tree)), sackin(as.treeshape(neut.tree)), tci(neut.tree), nLTTstat_exact(tree, neut.tree), nLTTstat_exact(tree, neut.tree, ignore_stem=F))
		
	}
	######################################################################################################################	
	if (u == 5){ ##habitat filtering OU
		print(5)
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
		bm.fit <- fitContinuous(habitat.tree, dat=xtips, model="BM", control=list(niter=10))
		ou.fit <- fitContinuous(habitat.tree, dat=xtips, model="OU", control=list(niter=10), bounds=list(alpha = c(min = exp(-500), max = 5.0)))
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
		
		stat.sim[i,16:24] <- c(mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z, colless(as.treeshape(habitat.tree)), sackin(as.treeshape(habitat.tree)), tci(habitat.tree), nLTTstat_exact(tree, habitat.tree), nLTTstat_exact(tree, habitat.tree, ignore_stem=F))
	}
	######################################################################################################################
	if (u == 6){ ##limited similarity OU
		print(6)
		#par <- runif(1,0,5.0)
		x <- fastBM(tree, a=0,theta=0, alpha=alpha, sig2=par, internal=FALSE)
		#hist(x, col="blue", breaks=10)
		alph<- abs(range(x)[1]-range(x)[2])/2
		limx <- c()
		Xi <- sample(x,1)
		limx <- c(Xi)
		x <- x[!x==Xi]
		#acc<-0
		#rej<-0
		while (length(limx) < n/2) {
			Xj <- sample(x,1)
			Probabilities <- rep(NA, length(limx))
			for (z in 1:length(limx)){
				Probabilities[z] <- (1 - exp(-(((limx[z]-Xj)^2)/(alph*tau))))
			}
			Pj <- mean(Probabilities)
			if (Pj > runif(1,0,1)) {
				limx <- c(limx, Xj)
				x <- x[!x==Xj]
				#acc <- acc +1
			}else{
				#x <- x[!x==Xj]
				#rej <- rej +1
			}	
		}
		#par
		#alph
		#acc/sum(acc,rej)
		#rej/sum(acc,rej)
		#par(new=T)
		#hist(limx, col="gray", xlim=c(-1,.5), ylim=c(0,50))
		
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
		bm.fit <- fitContinuous(limsim.tree, dat=xtips, model="BM", control=list(niter=10))
		ou.fit <- fitContinuous(limsim.tree, dat=xtips, model="OU", control=list(niter=10), bounds=list(alpha = c(min = exp(-500), max = 5.0)))
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
		stat.sim[i,16:24] <- c(mpd(mat,phydist)[1], mntd(mat, phydist)[1], ses.mpd$mpd.obs.z, ses.mntd$mntd.obs.z, colless(as.treeshape(limsim.tree)), sackin(as.treeshape(limsim.tree)), tci(limsim.tree), nLTTstat_exact(tree, limsim.tree), nLTTstat_exact(tree, limsim.tree, ignore_stem=F))
	}	
}


write.csv(par.sim, file="params_Sept28")
write.csv(stat.sim, file="stats_Sept28")
write.csv(com.mod, file="commod_Sept28")
write.csv(trait.mod, file="traitmod_Sept28")

