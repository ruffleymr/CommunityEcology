###################################################################
#This script simulates all of the regional species pool phylogenies. 
###################################################################
#################################
#1.) Load needed packages
#require(TreeSim)										#local machine
require(TreeSim, lib="/mnt/lfs2/ruff6699/Rpackages")	#ibest luster

#################################
#2.) Simulate Parameters by drawing from prior distributions
sims <- 60000	#Set the total number of trees to be simulated/number
n <- as.integer(runif(sims, 500, 1000))		#Vector of number of tips for each tree
lambda <- runif(sims, 0.05, 0.2)			#Vector of speciation rate, per million years
eps <- runif(sims, 0.2, 0.8) 				#Vector of eps; where eps is death over birth	
mu <- eps*lambda							#Vector of extinction rate, per million years
frac <- runif(sims, 0.5, 0.8)				#Vector of incomplete sampling fraction
parameters <- cbind(n, lambda, eps, mu, frac) #Store all parameter values and save
save(parameters, file="Large_RegPhy_Parameters.Rdata")

#################################
#3.) Simulate Regional Trees and save them all
trees <- c()		#empty vector to hold all simulated trees
for (i in 1:sims){
	tree <- sim.bd.taxa(n[i], 1, lambda[i], mu[i], frac=frac[i], complete=FALSE)
	trees <- c(trees, tree)
}
save(trees, file="Large_RegionalTrees.Rdata")

