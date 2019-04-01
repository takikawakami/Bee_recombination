# TODO: Add comment
# 
#
# After detecting CO and NCO by script CO_NCO_Queen_Drone_v1.R
#
# First we want to count COs, by applying the following filter
# Filter sets:
# 	1) removing regions <10kb and <10 SNPs
# 
# Next, we want to detect gene conversion events, one associated with CO events and the other NCOs
#
# version 4. Changes in NCO part.
# First, I want to keep track of NCOs supported by 1 SNP. This could well be a sequencing error, but could also be a real gene conversion. 
# To test if GC bases are preferentially transmitted in these small NCO regions, keep all putative NCOs.
#
# Second, scaffold coordinate will be kept in the final output for later analysis (the coordinate in the bam files are scaffold-based). 
#
# Author: takikawakami
# 
###############################################################################


#Clear the workspace
rm(list=ls())

library(stringr)
options(scipen=999)

comarg    <- commandArgs()   
geno.file <- comarg[6]              # geno.file <- "Group16.diploid.vcf.FB_B12.recode.4R.hap.out.hapBlock.ind.1"  
Nsnp_th <- as.integer(comarg[7])				# threshold for the number of snps, default 10
size_th <- as.integer(comarg[8])				# threshold for the distance between the beginning and end of haplotype (haplotype length)... default 10000

# setwd("/Users/takikawakami/Documents/Bee_Project/sandbag")

df <- read.table(geno.file, header=FALSE)
colnames(df) <- c("chrom", "POS1", "POS2", "hap", "Nsnp","scaf", "sPOS1", "sPOS2", "size")

subset_ct <- 1

# store lines here
subset.df <- c()

# filter input file by removing haplotype blocks composed of only N SNP.
Nsnp_th2 <- 1

for (i in 1:nrow(df)) {
	if (df$POS1[i]!=-9) {
		if (df$Nsnp[i]>Nsnp_th2) {		
			subset.df <- rbind(subset.df, df[i,])
		}
	} else if (df$POS1[i]==-9) {		# if the line is a gap, insert gap
		subset.df <- rbind(subset.df, df[i,])
	}
}



# collapse consequtive lines with the same haplotypes... function... NOTE: scaf coordinate is NOT considered here
collapse.hap <- function(input.df) {
	line_ct <-0		# stored line counter
	chrom <- as.character(input.df[1,1])
	POS1 <- c()
	POS2 <- c()
	hap.id <- c()
	Nsnp.sum <- c() 
	Nsnp.temp <- 0
	
	for (i in 1:nrow(input.df)) {
		if (line_ct==0 & input.df$POS1[i]==-9) {		# if the input.df starts from a gap, skip the line
			next
		} else if (line_ct==0) {						# very first line with data
			POS1 <- input.df$POS1[i]
			prev.hap.id <- as.character(input.df$hap[i])
			Nsnp.sum.temp <- input.df$Nsnp[i]
			line_ct <- 1
		} else if (input.df$POS1[i-1]==-9 & input.df$POS1[i]!=-9) {			# if right after the gap, reset data as well.
			POS1 <- c(POS1, input.df$POS1[i])
			prev.hap.id <- as.character(input.df$hap[i])
			Nsnp.sum.temp <- input.df$Nsnp[i]
			line_ct <- 1
		} else if (input.df$POS1[i]!=-9) {				# regular lines
			if (prev.hap.id==as.character(input.df$hap[i])) {		# if the current line has the same haplotype with the previous one, add it
				Nsnp.sum.temp <- Nsnp.sum.temp + input.df$Nsnp[i]
			} else if (prev.hap.id!=as.character(input.df$hap[i])) {
				# store data
				POS2 <- c(POS2, input.df$POS2[i-1])
				Nsnp.sum <- c(Nsnp.sum, Nsnp.sum.temp)
				hap.id <- c(hap.id, prev.hap.id)
				
				# reset data for the following loop
				POS1 <- c(POS1, input.df$POS1[i])
				prev.hap.id <- as.character(input.df$hap[i])
				Nsnp.sum.temp <- input.df$Nsnp[i]
				
				# counter
				line_ct <- line_ct + 1
				
			} else {
				print(i)
				writeLines("# Unexpected haplotype ")
			}
		} else if (input.df$POS1[i]==-9) {				# if the line is a gap, store data and reset data.
			if (input.df$POS1[i-1]==-9) {				# if the previous line is also a gap, do nothing
				next
			} else {
				# store data
				POS1 <- c(POS1, -9)						# insert gap
				POS2 <- c(POS2, input.df$POS2[i-1], -9)
				Nsnp.sum <- c(Nsnp.sum, Nsnp.sum.temp, -9)
				hap.id <- c(hap.id, prev.hap.id, -9)
				
				# counter
				line_ct <- line_ct + 1
			}
		} else {
			print(i)
			writeLines("# WTF!!! first hap line not found... Line 92")
		}
	}
	
	
	# after finishing loop, complete the POS2 and Nsnp.sum
	chrom.names <- rep(chrom, length(POS1))
	
	if (POS1[length(POS1)]!=-9) {
		POS2 <- c(POS2, input.df$POS2[i])
		hap.id <- c(hap.id, prev.hap.id)
		Nsnp.sum <- c(Nsnp.sum, Nsnp.sum.temp)
	}
	
	# summary df
	df.summary <- data.frame(chr=chrom.names, POS1=POS1, POS2=POS2, hap.id=hap.id, Nsnp=Nsnp.sum)
	df.summary$bp <- (df.summary$POS2 - df.summary$POS1)
	
	return(df.summary)
	
}

df.CO.NCO <- collapse.hap(subset.df)

# separate COs and NCOs
df.summary.NCO <- subset(df.CO.NCO, bp<size_th & bp>0)

df.summary.CO.infile <- subset(df.CO.NCO, bp>=size_th | bp==0)
df.summary.CO <- collapse.hap(df.summary.CO.infile)		# run the collapse function again after removing NCOs

# putative NCOs supported by only one SNP are added.
df.1snp.NCO <- subset(df, scaf!="gap" & size==0)
df.1snp.NCO <- df.1snp.NCO[,c(1:5)]

if (nrow(df.1snp.NCO)>0) {
	df.1snp.NCO$bp <- 1
	colnames(df.1snp.NCO) <- colnames(df.summary.NCO)
	df.summary.NCO <- rbind(df.summary.NCO, df.1snp.NCO )
	df.summary.NCO <- df.summary.NCO[ order(df.summary.NCO[,2]), ]
}


# determine whether given gene conversions are CO associated or not
# Also retrieve scaffold coordinate for NCOs

if (nrow(df.summary.NCO)>0) {
	# store test results for CO-associated NCO or not
	withCO.logical <- vector(length=nrow(df.summary.NCO))
	
	# store scaf coordinate
	scaf.out <- vector(length=nrow(df.summary.NCO))
	sPOS1 <- vector(length=nrow(df.summary.NCO))
	sPOS2 <- vector(length=nrow(df.summary.NCO))
	
	df.summary.CO.POS <- c(subset(df.summary.CO, POS1>0)$POS1, subset(df.summary.CO, POS2>0)$POS2)
	
	for (i in 1:nrow(df.summary.NCO)) {
		NCO.POS1 <- df.summary.NCO$POS1[i]
		NCO.POS2 <- df.summary.NCO$POS2[i]
		
		# if a NCO is within 10-kb (size_th) from COs, they are called 'CO-associated gene conversions' and marked as '1'. Else '0'.
		if (any(abs(df.summary.CO.POS - NCO.POS1) <size_th) | any(abs(df.summary.CO.POS - NCO.POS2) <size_th)) {
			withCO.logical[i] <- 1		# this is a CO-associated gene conversion
		} else {
			withCO.logical[i] <- 0		# this is a CO-associated gene conversion
		}
		
		# retrieve scaf coordinate
		df.sub <- subset(df, POS1<=NCO.POS1 & NCO.POS2<=POS2)
		
		if (nrow(df.sub)==1) {
			scaf.out[i] <-  as.character(df.sub$scaf)
			sPOS1[i] <- df.sub$sPOS1
			sPOS2[i] <- df.sub$sPOS2
			
		} else {
			scaf.list <- unique(as.character(df.sub$scaf))
			
			if (length(scaf.list)==1) {
				scaf.out[i] <-  unique(as.character(df.sub$scaf))
				sPOS1[i] <- df.sub$sPOS1[1]
				sPOS2[i] <- df.sub$sPOS2[nrow(df.sub)]
			} else {
				# This is the case where NCOs are detected at scaffold gaps. This should be filtered out in later analyses.
				scaf.out[i] <-  -9
				sPOS1[i] <- -9
				sPOS2[i] <- -9
			}
		}
	}
	
	# output
	df.summary.NCO$scaf <- scaf.out
	df.summary.NCO$sPOS1 <- sPOS1
	df.summary.NCO$sPOS2 <- sPOS2
	df.summary.NCO$COassoc <- withCO.logical	
	
} else {
	# if there is no NCO detected, make an empty output
	df.summary.NCO <- as.data.frame(matrix(c(as.character(df[1,1]), -9, -9, -9, -9, -9,-9,-9,-9,-9), nrow=1))
	colnames(df.summary.NCO) <- c("chr", "POS1", "POS2", "hap.id", "Nsnp", "bp", "scaf", "sPOS1", "sPOS2", "COassoc")
}


# write table
write.table(df.summary.NCO, file=paste(geno.file,".filter.",Nsnp_th,"SNP.",size_th/1000,"kb.NCO.out", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
write.table(df.summary.CO, file=paste(geno.file,".filter.",Nsnp_th,"SNP.",size_th/1000,"kb.CO.out", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
















## 
## # collapse consequtive lines with the same haplotypes... function
## collapse.hap <- function(input.df) {
##     line_ct <-0		# stored line counter
##     chrom <- as.character(input.df[1,1])
##     POS1 <- c()
##     POS2 <- c()
##     hap.id <- c()
##     Nsnp.sum <- c() 
##     Nsnp.temp <- 0
##     scaf <- c()
##     sPOS1 <- c()
##     sPOS2 <- c()	
## 
##     for (i in 1:nrow(input.df)) {
##         if (line_ct==0 & input.df$POS1[i]==-9) {		# if the input.df starts from a gap, skip the line
##             next
##         } else if (line_ct==0) {						# very first line with data
##             POS1 <- input.df$POS1[i]
##             sPOS1 <- input.df$sPOS1[i]
##             prev.hap.id <- as.character(input.df$hap[i])
##             Nsnp.sum.temp <- input.df$Nsnp[i]
##             line_ct <- 1
##         } else if (input.df$POS1[i-1]==-9 & input.df$POS1[i]!=-9) {			# if right after the gap, reset data as well.
##             POS1 <- c(POS1, input.df$POS1[i])
##             sPOS1 <- c(sPOS1, input.df$sPOS1[i])
##             prev.hap.id <- as.character(input.df$hap[i])
##             Nsnp.sum.temp <- input.df$Nsnp[i]
##             line_ct <- 1
##         } else if (input.df$POS1[i]!=-9) {				# regular lines
##             if (prev.hap.id==as.character(input.df$hap[i])) {		# if the current line has the same haplotype with the previous one, add it
##                 Nsnp.sum.temp <- Nsnp.sum.temp + input.df$Nsnp[i]
##             } else if (prev.hap.id!=as.character(input.df$hap[i])) {
##                 # store data
##                 POS2 <- c(POS2, input.df$POS2[i-1])
##                 sPOS2 <- c(sPOS2, input.df$sPOS2[i-1])
##                 Nsnp.sum <- c(Nsnp.sum, Nsnp.sum.temp)
##                 hap.id <- c(hap.id, prev.hap.id)
## 
##                 # reset data for the following loop
##                 POS1 <- c(POS1, input.df$POS1[i])
##                 sPOS1 <- c(sPOS1, input.df$sPOS1[i])
##                 prev.hap.id <- as.character(input.df$hap[i])
##                 Nsnp.sum.temp <- input.df$Nsnp[i]
## 
##                 # counter
##                 line_ct <- line_ct + 1
## 
##             } else {
##                 print(i)
##                 writeLines("# Unexpected haplotype ")
##             }
##         } else if (input.df$POS1[i]==-9) {				# if the line is a gap, store data and reset data.
##             if (input.df$POS1[i-1]==-9) {				# if the previous line is also a gap, do nothing
##                 next
##             } else {
##                 # store data
##                 POS1 <- c(POS1, -9)						# insert gap
##                 sPOS1 <- c(sPOS1, -9)						# insert gap
##                 POS2 <- c(POS2, input.df$POS2[i-1], -9)
##                 sPOS2 <- c(sPOS2, input.df$sPOS2[i-1], -9)
##                 Nsnp.sum <- c(Nsnp.sum, Nsnp.sum.temp, -9)
##                 hap.id <- c(hap.id, prev.hap.id, -9)
## 
##                 # counter
##                 line_ct <- line_ct + 1
##             }
##         } else {
##             print(i)
##             writeLines("# WTF!!! first hap line not found... Line 92")
##         }
##     }
## 
## 
##     # after finishing loop, complete the POS2 and Nsnp.sum
##     chrom.names <- rep(chrom, length(POS1))
## 
##     if (POS1[length(POS1)]!=-9) {
##         POS2 <- c(POS2, input.df$POS2[i])
##         sPOS2 <- c(sPOS2, input.df$sPOS2[i])
##         hap.id <- c(hap.id, prev.hap.id)
##         Nsnp.sum <- c(Nsnp.sum, Nsnp.sum.temp)
##     }
## 
##     # summary df
##     df.summary <- data.frame(chr=chrom.names, POS1=POS1, POS2=POS2, hap.id=hap.id, Nsnp=Nsnp.sum, scaf=scaf, sPOS1=sPOS1, sPOS2=sPOS2)
##     df.summary$bp <- (df.summary$POS2 - df.summary$POS1)
## 
##     return(df.summary)
## 
## }


