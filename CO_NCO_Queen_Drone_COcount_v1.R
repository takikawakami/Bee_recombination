# TODO: Add comment
# 
#
# After detecting CO and NCO by script CO_NCO_Queen_Drone_v1.R
# 
# We want to first get a big picture of CO events by filtering out unlikely NCO events
# Filter sets:
# 	1) removing regions <10kb and <10 SNPs
# 	2) removing regions <100kb and <50 SNPs
#	3) removing regions <1Mb and <500SNPs
# 
# Author: takikawakami
# 
###############################################################################


#Clear the workspace
rm(list=ls())

library(stringr)
options(scipen=999)

comarg    <- commandArgs()   
geno.for.file <- comarg[6]              # geno.for.file <- "Group11.diploid.vcf.FB_B12.recode.4R.hap.out.hapBlock.ind.1.filter.10SNP.10kb.out" 
#geno.rev.file <- comarg[7]              # geno.rev.file <- "Group9.diploid.vcf.FB_B12.recode.4R.rev.hap.out.hapBlock.ind.10.filter.10SNP.10kb.out" 

# setwd("/Users/takikawakami/Documents/Bee_Project/sandbag")

for.df <- read.table(geno.for.file, header=FALSE)
#rev.df <- read.table(geno.rev.file, header=FALSE)

colnames(for.df) <- c("chrom", "POS1", "POS2", "hap", "Nsnp", "size")
#colnames(rev.df) <- c("chrom", "POS1", "POS2", "hap", "Nsnp", "size")

if (nrow(for.df)>1) {
	# store lines here
	subset.for.df <- c()
	
	for (i in 1:nrow(for.df)) {
		if (for.df$POS1[i]!=-9) {
			if (i<nrow(for.df)) {
				if (for.df$POS1[i+1]!=-9) {		# if the line (haplotype block) has sufficient number of snps and size, keep it.
					subset.for.df <- rbind(subset.for.df, for.df[i,])
				} else if (for.df$POS1[i+1]==-9) {
					if (i>1) {
						if (for.df$POS1[i-1]!=-9) {
							subset.for.df <- rbind(subset.for.df, for.df[i,])
						}
					} 
				} else {
					print(i)
					writeLines("# Unexpected haplotype ")
					break
				}
			} else if (i==nrow(for.df)) {
				if (for.df$POS1[i-1]!=-9) {
					subset.for.df <- rbind(subset.for.df, for.df[i,])
				}
			}
		} else if (for.df$POS1[i]==-9) {		# if the line is a gap, insert gap
			subset.for.df <- rbind(subset.for.df, for.df[i,])
		}
	}
	
	
	### further refine subset.for.df by removing redundant gaps
	subset.for.df2 <- c()
	
# remove -9 if they are at the first consequtive lines
	start <- 1
	
	for (i in 1:nrow(subset.for.df)) {
		if (subset.for.df$POS1[i]==-9) {
			start <- start +1
		} else {
			break
		}
	} 
	
	#### if COs are observed ####
	if (start < nrow(subset.for.df)) {
		# scan the data, starting from the first non-gap line
		for (i in start:nrow(subset.for.df)) {
			if (subset.for.df$POS1[i]!=-9) {
				subset.for.df2 <- rbind(subset.for.df2, subset.for.df[i,])
			} else if (i<nrow(subset.for.df)) {
				if (subset.for.df$POS1[i+1]!=-9) {
					subset.for.df2 <- rbind(subset.for.df2, subset.for.df[i,])
				}
			}
		}
		
		# CO table
		chrom.name <- as.character(subset.for.df2[1,1])
		chrom <- rep(chrom.name, nrow(subset.for.df2)-1)
		POS1 <- subset.for.df2$POS2[1:(nrow(subset.for.df2)-1)]
		POS2 <- subset.for.df2$POS1[2:(nrow(subset.for.df2))]
		
		CO.for.df <- data.frame(chrom=chrom, POS1=POS1, POS2=POS2)
		CO.for.df <- subset(CO.for.df, POS1!=-9 & POS2!=-9)
		CO.for.df$dist <- CO.for.df$POS2 - CO.for.df$POS1	
	} else if (start > nrow(subset.for.df)) {	# there is no observable CO
		# CO table
		chrom.name <- as.character(for.df[1,1])
		chrom <- chrom.name
		POS1 <- -9
		POS2 <- -9
		
		CO.for.df <- data.frame(chrom=chrom, POS1=POS1, POS2=POS2)
		CO.for.df$dist <- CO.for.df$POS2 - CO.for.df$POS1
		
		# haplotype block table... empty
		subset.for.df2 <- data.frame(chrom=chrom, POS1=POS1, POS2=POS2, hap=-9, Nsnp=-9, size=-9)
		
		
	} else {
		writeLines("# Still error!!!! ")
	}
	
	
	
	
# write table
	write.table(subset.for.df2, file=paste(geno.for.file,".hapBlock.out", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
	write.table(CO.for.df, file=paste(geno.for.file,".CO.out", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
	
	
} else {
	# no observable CO by having only one entry
	# create an empty df with -9
	chrom.name <- as.character(for.df[1,1])
	chrom <- chrom.name
	POS1 <- -9
	POS2 <- -9	
	
	CO.for.df <- data.frame(chrom=chrom, POS1=POS1, POS2=POS2)
	CO.for.df$dist <- CO.for.df$POS2 - CO.for.df$POS1
	
	# haplotype block table... copy and paste the input file since this is the only one entry
	subset.for.df2 <- for.df
	
	write.table(subset.for.df2, file=paste(geno.for.file,".hapBlock.out", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
	write.table(CO.for.df, file=paste(geno.for.file,".CO.out", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
}

