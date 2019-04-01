# TODO: Add comment
# 
# THis is to calculate rho/kb in 100 bp bins
#
# version 2: PB vcf files contains multiple "Group"s. Say Group1_27, Group1_0, etc... These are accommodated.  
#
# Author: takikawakami
# 
###############################################################################


#Clear the workspace
rm(list=ls())

# library(stringr)
options(scipen=999)

comarg    <- commandArgs()   
geno.file <- comarg[6]              # geno.file <- "Group14.diploid.vcf.FB_B12.recode.4R" 

# setwd("/Users/takikawakami/Documents/Bee_Project/sandbag")

df <- read.table(geno.file, header=FALSE)

df.geno <- df[,c(5:ncol(df))]
group <- as.character(df$V1)
pos <- as.integer(df$V2)
ref.allele <- as.character(df$V3)
alt.allele <- as.character(df$V4)

# the number of individuals
Nind <- ncol(df.geno)

#snp1.info <- df[1,c(1:4)]
snp1 <- as.integer(df[1,c(5:ncol(df))])		# previous line

# initialize
hap1 <- 0
hap2 <- 1

# counter
gap_ct <- 0	# the number of haplotype gaps

# temporary end of haplotypes
hap1_end <- 0
hap2_end <- 1

# imputed genotypes
df.impute <- as.matrix(df[1,])

# removed genotypes
df.remove <- matrix(nrow = 0, ncol = Nind+4)

# function to compare two genotypes

compare.func <- function(line1, line2){
	
	mismatch_ct <- 0
	
	for (z in 1:Nind) {	
		if (line1[z]!=999 & line2[z]!=999) {
			if (line1[z]!=line2[z]) {
				mismatch_ct <- mismatch_ct +1
			}
		}
	}
	
	if (mismatch_ct>0) {
		return(1)	# if 1, there is one or more mismatching genotype
	} else {
		return(0)	# if 0, all individuals have identical genotype between line1 and line2
	}
}


for (i in 2:nrow(df.geno)) {
	#print(i)
	# set of snps at current row
	snpx <- as.numeric(df.geno[i,])
	snpx_next <- as.numeric(df.geno[i+1,])
	snpx_next2 <- as.numeric(df.geno[i+2,])
	
	# counter
	zerozero_ct <- 0
	oneone_ct <- 0
	zeroone_ct <- 0
	onezero_ct <- 0
	
	for (j in 1:ncol(df.geno)) {
		if (snp1[j] == 0 & snpx[j] == 0) {
			zerozero_ct <- zerozero_ct + 1
		} else if (snp1[j] == 1 & snpx[j] == 1) {
			oneone_ct <- oneone_ct + 1
		} else if (snp1[j] == 0 & snpx[j] == 1) {
			zeroone_ct <- zeroone_ct + 1
		} else if (snp1[j] == 1 & snpx[j] == 0) {
			onezero_ct <- onezero_ct + 1
		}
	}
	
	# make a data frame and 
	snp1.allele <- c(0,1,0,1)
	snpx.allele <- c(0,1,1,0)
	ct <- c(zerozero_ct, oneone_ct, zeroone_ct, onezero_ct)
	ct.df <- data.frame(snp1.allele=snp1.allele, snpx.allele=snpx.allele, ct=ct)
	ct.df <- ct.df[order(ct.df$ct),]
	
	# find the two most common combinations
	common1 <- as.numeric(ct.df[nrow(ct.df),])
	common_comp <- vector(length=3)
	
	# create complementary haplotype of common1 hap
	if (common1[1]==0) {
		common_comp[1] <- 1
	} else if (common1[1]==1) {
		common_comp[1] <- 0
	}
	
	if (common1[2]==0) {
		common_comp[2] <- 1
	} else if (common1[2]==1) {
		common_comp[2] <- 0
	}
	
	common_comp <- as.numeric(subset(ct.df, snp1.allele==common_comp[1] & snpx.allele==common_comp[2]))
	
	# count the number of expected and unexpected haplotypes by assuming that common1 is the correct one.
	exp_sum <- common1[3] + common_comp[3] 
	unexp_sum <- sum(ct.df$ct) - exp_sum
	
	
	########### if conditional upon the number of unexpected haplotypes ###########
	
	# make complement of the current line
	snp1_comp <- replace(snp1, snp1==0, 777)
	snp1_comp <- replace(snp1_comp, snp1_comp==1, 0)
	snp1_comp <- replace(snp1_comp, snp1_comp==777, 1)
	
	snpx_comp <- replace(snpx, snpx==0, 777)
	snpx_comp <- replace(snpx_comp, snpx_comp==1, 0)
	snpx_comp <- replace(snpx_comp, snpx_comp==777, 1)
	
	# if unexpected haplotypes >2, break haplotype (by inserting -9)
	if (unexp_sum>2) {
		if (nrow(df.geno)==i) {
			# If this is the last line with unexpected haplotype>2, remove the last line
			df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
			break
		} else if (nrow(df.geno)==i+1) {
			# If this is the 2nd from the last line with unexpected haplotype>2, remove the last line
			df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
			next
		} else {
			# if the following line "snp_next" is identical to the current line...
			if (compare.func(snpx, snpx_next)==0 || compare.func(snpx_comp, snpx_next)==0) {
				if (compare.func(snp1, snpx_next2)==0 || compare.func(snp1_comp, snpx_next2)==0) {
					# if the 3rd line from the previous line (snp1) is identical to snp1, remove the current line
					
					df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
					next
				} else if (compare.func(snpx, snpx_next2)==0 || compare.func(snpx_comp, snpx_next2)==0) {
					# if current line and the following 2 lines are identical, insert haplotype gap 
					
					hap1 <- c(hap1,-9, 0)
					hap2 <- c(hap2,-9, 1)
					df.impute <- rbind(df.impute, c("gap", -9, -9, -9, rep(-9,Nind)), c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
					
					# counter
					gap_ct <- gap_ct + 1
					
					# renew the previous SNP line
					snp1 <- snpx
					
					# renew hap end snps... restart haplotype reconstruction
					hap1_end <- 0
					hap2_end <- 1
					
					next
				} else {
					# if the 3rd line from the previous line (snp1) is identical to snp1, remove the current line
					
					df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
					next
				}
			} else {			# if not, remove the current line, and go back to the loop head
				df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
				
				next
			}
		}
			
	} else if (unexp_sum==2) {
		if (nrow(df.geno)==i) {
			# If this is the last line with unexpected haplotype>2, remove the last line
			df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
			break
		} else if (nrow(df.geno)==i+1) {
			# If this is the 2nd from the last line with unexpected haplotype>2, remove the last line
			df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
			next
		} else {
			# if the following line is different from the current line, remove the current line
			
			if (compare.func(snpx, snpx_next)==0 || compare.func(snpx_comp, snpx_next)==0) {
				# if the following line is identical to the current line, move on.
			} else {
				# remove current line
				df.remove <- rbind(df.remove, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
				
				next
			}
		}
	}
	
	
	
	
	# most common snp pair	
	if (common1[1] == 0) {		# if previous hap1 or hap2
		if (hap1_end == 0) {
			hap1 <- c(hap1,common1[2])
			hap2 <- c(hap2,common_comp[2])
		} else if (hap2_end == 0) {
			hap2 <- c(hap2,common1[2])
			hap1 <- c(hap1,common_comp[2])
		}
	} else if (common1[1] == 1) {
		if (hap1_end == 1) {
			hap1 <- c(hap1,common1[2])
			hap2 <- c(hap2,common_comp[2])
		} else if (hap2_end == 1) {
			hap2 <- c(hap2,common1[2])
			hap1 <- c(hap1,common_comp[2])
		}
	} else {
		writeLines("# WTF!!! haplotype not found!!!!")
	}
	
	### impute missing genotypes
	if (999 %in% snpx) {			# if snpx contains missing 999
		
		# identify which individual contain missing value
		missing.index <- which(snpx %in% 999)
		
		# get allele from snp1 list
		snp1.allele <- snp1[missing.index]
		
		if (snp1.allele == common1[1]) {
			snpx[missing.index] <- common1[2]
		} else if (snp1.allele == common_comp[1]) {
			snpx[missing.index] <- common_comp[2]
		}
		
		# add line
		df.impute <- rbind(df.impute, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
	} else {								# if there is no missing data
		df.impute <- rbind(df.impute, c(group[i], pos[i], ref.allele[i], alt.allele[i], snpx))
	}
	
	# renew the previous SNP line
	snp1 <- snpx
			
	# renew hap end snps
	hap1_end <- hap1[length(hap1)]
	hap2_end <- hap2[length(hap2)]
}



### if the first line contains missing values
df.impute <- cbind(df.impute, hap1, hap2)

if (999 %in% df.impute) {
	missing.mat <- which(df.impute == 999, arr.ind = T)
	rownames(missing.mat) <- seq(1,nrow(missing.mat))
	
	missing.pos <- as.data.frame(missing.mat, row.names=NULL, col.names=NULL)
	missing.pos <- missing.pos[order(-missing.pos$row),]	#descending order by adding minus sign
	
	# scan 
	for (l in 1:nrow(missing.pos)) {
		missing.row <- missing.pos$row[l]
		missing.ind <- missing.pos$col[l]
		
		if (missing.row!=nrow(df.impute)) {		# if the line with missing data is the very last one, skip it
			missing.allele_next <- df.impute[missing.row+1,missing.ind]
			hap1.allele_next <- df.impute[missing.row+1, Nind+5]
			hap2.allele_next <- df.impute[missing.row+1, Nind+6]
			
			if (missing.allele_next!=-9 & missing.allele_next!=999) {		# if the one line below is a gap (with continuous 999 missing genotypes), do nothing.
				if (missing.allele_next == hap1.allele_next) {
					df.impute[missing.row, missing.ind] <- df.impute[missing.row, Nind+5]
				} else if (missing.allele_next == hap2.allele_next) {
					df.impute[missing.row, missing.ind] <- df.impute[missing.row, Nind+6]
				}
			}	
		}
	}
}


### combine df.impute and queen haplotypes

write.table(df.impute, file=paste(geno.file,".hap.out",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

if (nrow(df.remove)>0) {
	write.table(df.remove, file=paste(geno.file,".remove.out",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
}


### outout gap positions

gap.summary <- matrix(nrow = 0, ncol = 3)

if ("gap" %in% df.impute) {
	gap.pos.mat <- which(df.impute == "gap", arr.ind = T)
	rownames(gap.pos.mat) <- seq(1,nrow(gap.pos.mat))
	
	gap.pos.df <- as.data.frame(gap.pos.mat, row.names=NULL, col.names=NULL)
	
	for (m in 1:nrow(gap.pos.df)) {
		gap.centre <- gap.pos.df[m,1]
		
		gap.chrom.temp <- df.impute[c(gap.centre-1,gap.centre+1),1]
		if (gap.chrom.temp[1]==gap.chrom.temp[2]) {
			gap.chrom <- gap.chrom.temp[1]
		} else {
			gap.chrom <- paste(gap.chrom.temp[1], gap.chrom.temp[2], sep=":")
		}
		
		gap.pos <- df.impute[c(gap.centre-1,gap.centre+1),2]
		gap.summary <- rbind(gap.summary, c(gap.chrom,gap.pos))
	}
}

write.table(gap.summary, file=paste(geno.file,".gap.out",sep=""), row.names=F, col.names=FALSE, sep="\t", quote=F)

### output summary result

df.impute <- as.data.frame(df.impute, row.names=FALSE, col.names=NULL)

df.summary <- data.frame(snpIN=nrow(df), snpOUT=nrow(subset(df.impute, V1!="gap")), removed=nrow(df.remove), gap=gap_ct)
write.table(df.summary, file=paste(geno.file,".summary.out",sep=""), row.names=F, col.names=TRUE, sep="\t", quote=F)

