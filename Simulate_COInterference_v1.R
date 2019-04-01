# TODO: Add comment
# 
# Strength of CO interference (nu) can be estimated by fitting distribution of COs by gamma model (R package 'xoi'). 
# I want to know how the number of COs per chromosome affect the estimate of nu while maintaining prior value of nu (or m when simulating COs by package 'qtl').
#
# infile 
# gendist... total genetic distance (ie, recombination rate with fixed number of markers)
# rep ... number of replicates
# int ... strength of crossover interference
#
# Author: takikawakami
###############################################################################



#Clear the workspace
rm(list=ls())

# library(stringr)
# options(scipen=999)
# library(simpleboot)
library(devtools)
library(xoi)


comarg    <- commandArgs()   
gendist <- as.integer(comarg[6])              # total genetic distance to simulate. Average genetic distance across 16 honeybee chromosome is about 300 cM.
rep <- as.integer(comarg[7])    			#  the number of replicates
int <- as.numeric(comarg[8])    			#  the number of replicates

df.out <- data.frame()

for (i in 1:rep) {
	
# 1. simulate genetic map. genetic distance = 100 cM, the number of markers = 101
	map1 <- sim.map(gendist, n.mar=3001, anchor=TRUE, include.x=FALSE, eq=TRUE)
	
# 2. simulate backcross samples with the number of individuals = 10, interference m = 6
	x <- sim.cross(map1, n.ind=1000, m=int, type="bc")
	
# 3. Convert the format of data on crossover locations to that needed for the function fitGamma.
	
	xoloc <- convertxoloc(find.breaks(x))
	
	
	### run fitGamma to estimate nu
# point estimate of nu
	fG.all.mle <- fitGamma(xoloc, lo=1, hi=10, se=FALSE)
	
# store results
	df.out <- rbind(df.out, fG.all.mle)
	
}


write.table(df.out, file=paste("xoi.interferece.sim.", gendist, "cM.", int, "int.", rep, "rep.out", sep=""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")






