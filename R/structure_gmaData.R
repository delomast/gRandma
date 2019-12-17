#' implementing a structure to hold data
#' 
#' baseline: the dataframe of the baseline to use for inferring relationships, with genotypes
#'   recoded to integers
#' mixture: the dataframe of the mixture to use for inferring relationships, with genotypes
#'   recoded to integers
#' baselineParams: A list (on entry per baseline population) of lists of vectors 
#'   contining alphas for Dirichlet distributions for allele frequencies
#' unsamParams:  A list (on entry per baseline population) of lists of vectors 
#'   contining alphas for Dirichlet distributions for allele frequencies in the unsampled breeding 
#'   population where the corresponding baseline would be found
#' lociErr: a list of matrices (one for each locus). Each matrix gives the probability
#'   if the true genotype is X (rows), the observed genotype is Y. So, rows should always sum to 1.
#'   entries are named for the locus, rows and columns are named for the genotype. Row and column 
#'   names are just added to make human reading/editing easier. the order is what's used by the 
#'   internal functions. 
#' @keywords internal
#' @noRd

construct_grandma <- function(x){

	# do some checking here of required values
	gmaDataNames <- c("baseline", "mixture", "baselineParams", "unsamParams", "lociErr")
	if(sum(names(x) != gmaDataNames) > 0) stop("Wrong names to make a gmaData object")
	
	class(x) <- "gmaData"
	return(x)
}

#' print method for gmaData
gmaData.print <- function(x){
	cat("\nclass of gmaData with\n")
	
	#some more here
	
	
}

#' This creates a gmaData structure from user supplied dataframes
#' 
#' @param baseline the baseline individuals to use to infer relationships. The first column should
#'   be the population each individual is from. The second column is the individual's identifier.
#'   Following columns are genotypes, with columns 3 adn 4 being Allele1 and 2 for locus 1, columns
#'   5 and 6 beign locus 2, ... This is a "two column per call" type of organization. Genotypes 
#'   should be given as concatenated basecalls, with each base represented by a single character. Some examples of
#'   microhaplotype alleles are "ACA", "AAD", "CTCTGGA". Some SNP alleles (which are just microhaplotypes
#'   with a length of one base) are "A", "G", "D". In these examples, deletions are represented by "D". Missing genotypes
#'   should be NA.
#' @param mixture the mixture individuals to infer relationships for. The first column  is the individual's identifier.
#'   Following columns are genotypes, in the same manner as for \code{baseline}. The order and column names of the loci 
#'   must be the same in both the baseline and mixture dataframes.
#' @param unsampliedPops
#' @param perSNPerror list of length equal to number of loci. Each entry is vector of error rates for each SNP in a locus.
#' @param dropoutProb dataframe with column 1 being locus name, column 2 being allele, and column 3 being dropout probability.
#'   the locus name must be the column name of allele 1 for the corresponding locus in the baseline and mixture dataframes
#' @importFrom Rcpp evalCpp
#' @useDynLib gRandma, .registration=TRUE
#' @export

createGmaInput <- function(baseline, mixture, unsampledPops = NULL, perSNPerror = NULL, dropoutProb = NULL){
	
	#type check
	if (!is.data.frame(baseline)){
		warning("Coercing baseline to a dataframe")
		baseline <- as.data.frame(baseline)
	}
	if (!is.data.frame(mixture)){
		warning("Coercing mixture to a dataframe")
		mixture <- as.data.frame(mixture)
	}
	
	# check that if one is NA, all are NA
	
	# first, get names of all markers
	markers <- colnames(baseline)[seq(3, (ncol(baseline) - 1), 2)]
	markersMix <- colnames(mixture)[seq(2, (ncol(mixture) - 1), 2)]
	if(markers != markersMix) stop("error, the loci are either not the same or not in the same order between the baseline and mixture.")
	if(names(perSNPerror) != markers){
		warning("reorering perSNPerror to match the order of the markers in baseline and mixture")
		perSNPerror <- perSNPerror[markers]
	}
	
	
	
	# now, for each locus
	for(m in seq(3, (ncol(baseline) - 1), 2)){
		
		# add default options for perSNPerror and dropoutProb
		
		# get all alleles
		alleles <- sort(unique(c(baseline[,m], baseline[,m+1], mixture[,m-1], mixture[,m])))
		alleles <- alleles[!is.na(alleles)]
		lenLocus <- unique(nchar(alleles))
		if(length(lenLocus) != 1) stop("not all alleles are the same length for locus ", colnames(baseline)[m])
		if(lenLocus != length(perSNPerror[[((m-1)/2)]])) stop("not all SNPs have error rates for locus ", colnames(baseline)[m])
		
		numAlleles <- length(alleles)
		# assign ints to each allele
		alleleKey <- data.frame(alleles = alleles, intAlleles = 1:numAlleles) # dataframe b/c combining char and num vectors
		
		#recode genotypes
		baseline[,m] <- recodeAlleles(baseline[,m], alleleKey)
		baseline[,m+1] <- recodeAlleles(baseline[,m+1], alleleKey)
		mixture[,m-1] <- recodeAlleles(mixture[,m-1], alleleKey)
		mixture[,m] <- recodeAlleles(mixture[,m], alleleKey)
		#flip alleles to all be consistent
		baseline[,m:(m+1)] <- flipHets(baseline[,m:(m+1)])
		mixture[,(m-1):m] <- flipHets(mixture[,(m-1):m])
		
		#recode dropoutProb
		tempDropoutProb <- dropoutProb[dropoutProb[,1] == colnames(baseline)[m],]
		tempDropoutProb[,2] <- recodeAlleles(tempDropoutProb[,2], alleleKey)
		
		# calc per allele error rates
		alleleErr <- matrix(nrow = numAlleles, ncol = numAlleles)
		for(i in 1:(numAlleles - 1)){
			for(j in (i+1):numAlleles){
				tempErr <- 1
				for(c in 1:lenLocus){
					if(substr(alleles[i],c,c) != substr(alleles[j],c,c)){
						tempErr <- tempErr * perSNPerror[[((m-1)/2)]][c]
					} else {
						tempErr <- tempErr * (1 - perSNPerror[[((m-1)/2)]][c])
					}
				}
				
				alleleErr[i,j] <- tempErr
				alleleErr[j,i] <- tempErr
			}
		}
		diag(alleleErr) <- prod(1 - perSNPerror[[((m-1)/2)]])

		# now per genotype error rates
		# general models:
		# for a true heterozygote, AB: 
		#  P(CD|AB) = P(dropout of A) * P(CD|AA) + P(dropout of B) * P(CD|BB) +  P(CD|AB, no dropout)
		# where,
		#  P(CD|AB, no dropout) = (1 - P(dropout of A) - P(dropout of B)) * (P(A obs as C) * P(B obs as D) + P(A obs as D) * P(B obs as C).
		# for a true homozygote, AA:
		#  P(CD|AA) = P(A obs as C) * P(A obs as D) * 2
		
		# assign ints to each genotype
		genotypeKey <- matrix(0,0,2)
		for(i in 1:numAlleles){ #build all combinations
			genotypeKey <- rbind(genotypeKey, cbind(i,1:numAlleles))
		}
		#remove repetitive entries
		genotypeKey <- genotypeKey[genotypeKey[,1] <= genotypeKey[,2],]

		#calculate genotype to genotype error probabilities
		# rows are true genotype, columns are observed genotype, entries are probability of observing given true
		genoErr <- matrix(nrow = nrow(genotypeKey), ncol = nrow(genotypeKey))
		# first calculate all for true homozygotes (these are used to then calculate for true heterozygotes)
		for(i in which(genotypeKey[,1] == genotypeKey[,2])){ #for each true genotype
			for(j in 1:nrow(genotypeKey)){ #for each possible observed genotype
				trueA <- genotypeKey[i,1]
				obsC <- genotypeKey[j,1]
				obsD <- genotypeKey[j,2]
			
				
				if(obsC != obsD){
					# 2 * P(A obs as C) * P(A obs as D)
					genoErr[i,j] <- 2 * alleleErr[trueA,obsC] * alleleErr[trueA,obsD]
				} else {
					# P(A obs as C)^2
					genoErr[i,j] <- alleleErr[trueA,obsC]^2
				}
			}
		}
		# then calculate for true heterozygotes
		
		############
		## something messed up here
		############
		
		for(i in which(genotypeKey[,1] != genotypeKey[,2])){ #for each true genotype
			for(j in 1:nrow(genotypeKey)){ #for each possible observed genotype
				trueA <- genotypeKey[i,1]
				trueB <- genotypeKey[i,2]
				obsC <- genotypeKey[j,1]
				obsD <- genotypeKey[j,2]
				
				if(obsC == obsD && (obsC == trueA || obsC == trueB)){
					# 	P(AA|AB, no dropout) = (1 - P(dropout of A) - P(dropout of B)) * (P(A obs as A) * P(B obs as C)
					PnoDropout <- (1 - tempDropoutProb[tempDropoutProb[,2] == trueA,3] - tempDropoutProb[tempDropoutProb[,2] == trueB,3]) *
						(alleleErr[trueA,obsC] * alleleErr[trueB,obsD])
				} else {
					# P(CD|AB, no dropout) = (1 - P(dropout of A) - P(dropout of B)) * (P(A obs as C) * P(B obs as D) + P(A obs as D) * P(B obs as C)
					PnoDropout <- (1 - tempDropoutProb[tempDropoutProb[,2] == trueA,3] - tempDropoutProb[tempDropoutProb[,2] == trueB,3]) *
						(alleleErr[trueA,obsC] * alleleErr[trueB,obsD] + alleleErr[trueA,obsD] * alleleErr[trueB,obsC])
				}

			

				trueAA <- which(genotypeKey[,1] == trueA & genotypeKey[,2] == trueA)
				trueBB <- which(genotypeKey[,1] == trueB & genotypeKey[,2] == trueB)
				# P(CD|AB) = P(dropout of A) * P(CD|BB) + P(dropout of B) * P(CD|AA) +  P(CD|AB, no dropout)
				genoErr[i,j] <- (tempDropoutProb[tempDropoutProb[,2] == trueA,3] * genoErr[trueBB,j]) +
					(tempDropoutProb[tempDropoutProb[,2] == trueB,3] * genoErr[trueAA,j]) + 
					PnoDropout
		}
		}
		
		#testing
		return(genoErr)
		
		
		#recode genotypes
		
		### lots of nested for loops.....
		###  may have to move some/all of this to c++
		
				
		# recode genotypes
		newBase <- baseline[,i:(i+1)]
		newMix <- mixture[,(i-1):i]
		

		
	}
	
}
