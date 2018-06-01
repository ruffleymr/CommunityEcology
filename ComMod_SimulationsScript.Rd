#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## 1. load packages and data
library(ape)			
library(phytools)
library(geiger)	
library(psych)

setwd("/Users/Megan/Documents/EcolLetters_Data")
load(file="RegTrees40.Rdata") ##reg.pars or regphy.params
load(file="RegPars40.Rdata") ##reg.trees

setwd("/Users/Megan/Documents/EcolLetters_Data/EmpSims_52918")
load(file="RegParams.Rdata")
load(file="RegTrees.Rdata")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## 2. set prior distribution on parameters and create storage lists

#number of total simulations
sims <- 10000

#random generator of model number
mod <- round(runif(sims, .5, 6.5))

#prior on rate of trait evolution for BM and OU
sig2 <- runif(sims, .05, 5)

#prior on sixe of the community
#com.size <- runif(sims, 0.3, 0.7) 
com.size <- rep(0.53, sims)

#prior on alpha for OU simulations
alpha <- runif(sims, 0.01, 0.2)

#prior on strength of filtering/competition
tau <- rlnorm(sims, meanlog=3.5, sd=.5)

##Global Lists that keep track of parameters, variable, trees, traits
Com.Params <- matrix(NA, sims, 7) 
colnames(Com.Params) <- c("Model", "sig^2", "com.size", "n", "tau", "alpha", "rej")

Reg.Traits <- list()  #to be filled by 'reg.traits'
Com.Trees <- list()	  #to be filled by 'com.tree'
Com.Traits <- list()  #to be filled by 'com.traits'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## 3. Simulate data

#will need this keep.tip function
keep.tip <- function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

for (i in 5179:sims) {
	n <- round(length(reg.trees[[i]]$tip.label) * com.size[i])
	rej <- 0
	
	##Mod 1 Brownian Motion, neutral
	if (mod[i] == 1){
		traits <- sim.char(reg.trees[[i]], par=sig2[i], nsim=1, model="BM", root=0)[,,1]
		Reg.Traits[[i]] <- traits
		com.traits <- sample(traits, n)
	}
	
	##Mod 2 Brownian Motion, habitat filtering
	if (mod[i] == 2){
		traits <- sim.char(reg.trees[[i]], par=sig2[i], nsim=1, model="BM", root=0)[,,1]
		Reg.Traits[[i]] <- traits

		tr.ord <- sort(traits)
		if (runif(1,0,1) <= 0.5){
			Xi <- sample(tr.ord[1:as.integer(length(traits)/3)],1)
		}else{
			Xi <- sample(tr.ord[(length(traits)-as.integer(length(traits)/3)):length(traits)],1)
		}
		com.traits <- Xi
		traits <- traits[!traits==Xi]
		while (length(com.traits) < n) {
			Xj <- sample(traits,1)
			Probabilities <- rep(NA, length(com.traits))
			for (z in 1:length(com.traits)){
				Probabilities[z] <- (exp(-((com.traits[z]-Xj)^2)/tau[i]))
			}
			Pj <- mean(Probabilities)
			if (Pj > runif(1,0,1)) {
				com.traits <- c(com.traits, Xj)
				traits <- traits[!traits==Xj]
			}else{
				rej=rej+1
			}
		}
	}
	
	##Mod 3 Brownian motion, competitive exclusion
	if (mod[i] == 3){
		traits <- sim.char(reg.trees[[i]], par=sig2[i], nsim=1, model="BM", root=0)[,,1]
		Reg.Traits[[i]] <- traits

		Xi <- sample(traits, 1)
		com.traits <- Xi
		traits <- traits[!traits==Xi]
		while (length(com.traits) < n) {
			Xj <- sample(traits,1)
			Probabilities <- rep(NA, length(com.traits))
			for (z in 1:length(com.traits)){
				Probabilities[z] <- (1 - exp(-(((com.traits[z]-Xj)^2)/tau[i])))
			}
			Pj <- mean(Probabilities)
			if (Pj > runif(1,0,1)) {
				com.traits <- c(com.traits, Xj)
				traits <- traits[!traits==Xj]
			}else{
				rej=rej+1
			}
		}		
	}
	
	##Mod 4 Ornstein-Uhlenbeck, neutral
	if (mod[i] == 4){
		traits <- fastBM(reg.trees[[i]], a=0, theta=0, alpha=alpha[i], sig2=sig2[i], internal=FALSE)
		Reg.Traits[[i]] <- traits
		com.traits <- sample(traits, n)
	}
	
	##Mod 5 Ornstein-Uhlenbeck, habitat filtering
	if (mod[i] == 5){
		traits <- fastBM(reg.trees[[i]], a=0, theta=0, alpha=alpha[i], sig2=sig2[i], internal=FALSE)
		Reg.Traits[[i]] <- traits

		tr.ord <- sort(traits)
		if (runif(1,0,1) <= 0.5){
			Xi <- sample(tr.ord[1:as.integer(length(traits)/3)],1)
		}else{
			Xi <- sample(tr.ord[(length(traits)-as.integer(length(traits)/3)):length(traits)],1)
		}
		com.traits <- Xi
		traits <- traits[!traits==Xi]
		while (length(com.traits) < n) {
			Xj <- sample(traits,1)
			Probabilities <- rep(NA, length(com.traits))
			for (z in 1:length(com.traits)){
				Probabilities[z] <- (exp(-((com.traits[z]-Xj)^2)/tau[i]))
			}
			Pj <- mean(Probabilities)
			if (Pj > runif(1,0,1)) {
				com.traits <- c(com.traits, Xj)
				traits <- traits[!traits==Xj]
			}else{
				rej=rej+1
			}
		}
	}
	
	if (mod[i] == 6){
		traits <- fastBM(reg.trees[[i]], a=0, theta=0, alpha=alpha[i], sig2=sig2[i], internal=FALSE)
		Reg.Traits[[i]] <- traits

		Xi <- sample(traits, 1)
		com.traits <- Xi
		traits <- traits[!traits==Xi]
		while (length(com.traits) < n) {
			Xj <- sample(traits,1)
			Probabilities <- rep(NA, length(com.traits))
			for (z in 1:length(com.traits)){
				Probabilities[z] <- (1 - exp(-(((com.traits[z]-Xj)^2)/tau[i])))
			}
			Pj <- mean(Probabilities)
			if (Pj > runif(1,0,1)) {
				com.traits <- c(com.traits, Xj)
				traits <- traits[!traits==Xj]
			}else{
				rej=rej+1
			}
		}
	}
	
	print(i)
	Com.Trees[[i]] <- keep.tip(reg.trees[[i]], names(com.traits))
	Com.Params[i,] <- c(paste("mod",mod[i], sep=""), sig2[i], com.size[i], n, tau[i], alpha[i], rej)
	Com.Traits[[i]] <- com.traits
}

#setwd("/Users/Megan/Documents/EcolLetters_Data/Sims_52318/")

save(Com.Trees, file="ComTrees.Rdata")
save(Com.Params, file="ComParams.Rdata")
save(Reg.Traits, file="RegTraits.Rdata")
save(Com.Traits, file="ComTraits.Rdata")

load(file="SumStats.Rdata")