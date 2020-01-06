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
	gmaDataNames <- c("baseline", "mixture", "unsampledPops", "genotypeErrorRates", "genotypeKeys", 
							"alleleKeys", "baselineParams", "unsampledPopsParams")
	if(sum(names(x) != gmaDataNames) > 0) stop("Wrong names to make a gmaData object")
	
	class(x) <- "gmaData"
	return(x)
}

#' print method for gmaData
print.gmaData <- function(x){
	cat("\ngmaData object with\n\t")
	cat(length(unique(x$baseline[,1])), " baseline population(s)\n\t")
	if(is.null(x$unsampledPops)){
		cat("No unsampledPops\n\t")
	} else{
		cat(length(unique(x$unsampledPops[,1])), " unsampled population(s)\n\t")
	}
	cat(nrow(x$mixture), " mixture individuals\n\t")
	cat(length(x$genotypeErrorRates), " loci\n")
	cat("The number of loci by number of unique alleles is:\n")
	print(table(sapply(x$genotypeKeys, function(y){
		return(length(unique(y[,2])))
	})))
	cat("\n\n")
	
}

#' This creates a gmaData structure from user supplied dataframes
#' 
#' @param baseline a dataframe of the baseline individuals to use to infer relationships. The first column should
#'   be the population each individual is from. The second column is the individual's identifier.
#'   Following columns are genotypes, with columns 3 adn 4 being Allele1 and 2 for locus 1, columns
#'   5 and 6 beign locus 2, ... This is a "two column per call" type of organization. Genotypes 
#'   should be given as concatenated basecalls, with each base represented by a single character. Some examples of
#'   microhaplotype alleles are "ACA", "AAD", "CTCTGGA". Some SNP alleles (which are just microhaplotypes
#'   with a length of one base) are "A", "G", "D". In these examples, deletions are represented by "D". Missing genotypes
#'   should be NA.
#' @param mixture a dataframe of the mixture individuals to infer relationships for. The first column  is the individual's identifier.
#'   Following columns are genotypes, in the same manner as for \code{baseline}. The order and column names of the loci 
#'   must be the same in both the baseline and mixture dataframes.
#' @param unsampledPops a dataframe of the individuals sampled from the "unsampled populations" used to estimate allele frequencies 
#'   in these populations. See details for more explanation. Column 1 has the baseline population that an individual corresponds to, 
#'   column 2 is the individual's identifier
#'   (not currently used, but must be present). The following columns are genotypes, in the same manner as for \code{baseline}. 
#'   The order and column names of the loci must be the same as in the baseline and mixture dataframes.
#' @param perSNPerror dataframe with each row representing a SNP. Column 1 is the locus name, column 2 is the order 
#'   of the SNP in the locus (ie 1, 2, 3, ...), column 3 is the number of alleles at that SNP (not the locus as a whole),
#'   and column 4 is the error rate (probability of observing any allele other than the correct allele at that SNP - not the 
#'   locus as a whole).
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
	if(!is.null(unsampledPops) && is.na(unsampledPops)) unsampledPops <- NULL
	useUnsamp <- !is.null(unsampledPops)
	if (useUnsamp & is.data.frame(unsampledPops)){
		warning("Coercing mixture to a dataframe")
		unsampledPops <- as.data.frame(unsampledPops)
	}
	
	# check that if one is NA, all are NA?
	
	# first, get names of all markers
	markers <- colnames(baseline)[seq(3, (ncol(baseline) - 1), 2)]
	markersMix <- colnames(mixture)[seq(2, (ncol(mixture) - 1), 2)]
	if(any(markers != markersMix)) stop("error, the loci are either not the same or not in the same order between the baseline and mixture.")
	if(any(!(markers %in% perSNPerror[,1]))) stop("not all markers have entries in perSNPerror")
	if(any(perSNPerror[,3] < 2)) stop("all SNPs in perSNPerror must have 2 or more alleles")
	if(useUnsamp){
		markersUnsamp <- colnames(unsampledPops)[seq(3, (ncol(unsampledPops) - 1), 2)]
		if(markers != markersUnsamp) stop("error, the loci are either not the same or not in the same order between the baseline and unsampledPops.")
	}
	
	# set up storage
	genotypeErrorList <- list() # list of matrices giving P(obs geno | true geno) for each locus
	genotypeKeyList <- list() # key of genotype code and alleles for each locus
	alleleKeyList <- list()
	baselineParams <- list() # list, entry for each population, of lists, for each locus, a named vector giving Dirichlet posterior of
		# allele frequencies given a 1/N prior, N = number of alleles
	unsampledPopsParams <- list() # same as baselineParams, but for unsampledPops
	
	# initiate storage for different pops
	basePops <- unique(baseline[,1])
	for(p in basePops) baselineParams[[as.character(p)]] <- list()
	if(useUnsamp){
		unsamPops <- unique(unsampledPops[,1])
		for(p in unsamPops) unsampledPopsParams[[as.character(p)]] <- list()
	}
	
	# now, for each locus
	for(m in seq(3, (ncol(baseline) - 1), 2)){
		### lots of nested for loops.....
		###  may have to move some/all of this to c++
		mName <- colnames(baseline)[m]
		# add default options for perSNPerror and dropoutProb ??
		
		tempPerSNPerror <- perSNPerror[perSNPerror[,1] == mName,]
		tempPerSNPerror <- tempPerSNPerror[order(tempPerSNPerror[,2]),]
		
		# get all alleles
		alleles <- sort(unique(c(baseline[,m], baseline[,m+1], mixture[,m-1], mixture[,m])))
		if(useUnsamp) alleles <- sort(unique(c(alleles, unsampledPops[,m], unsampledPops[,m+1])))
		alleles <- alleles[!is.na(alleles)]
		
		# check that locus is variable and not all missing
		if(length(alleles) == 0) stop("no genotypes for locus ", mName)
		if(length(alleles) == 1) stop("no variation for locus ", mName)
		
		lenLocus <- unique(nchar(alleles))
		if(length(lenLocus) != 1) stop("not all alleles are the same length for locus ", mName)
		if(lenLocus < nrow(tempPerSNPerror)) stop("not all SNPs have error rates for locus ", mName)
		if(lenLocus > nrow(tempPerSNPerror)) stop("more error rates for locus ", mName, " than there are SNPs")
		
		numAlleles <- length(alleles)
		# assign ints to each allele
		alleleKey <- data.frame(alleles = alleles, intAlleles = 1:numAlleles, stringsAsFactors = FALSE) # dataframe b/c combining char and num vectors
		
		# recode alleles
		baseline[,m] <- recodeAlleles(baseline[,m], alleleKey)
		baseline[,m+1] <- recodeAlleles(baseline[,m+1], alleleKey)
		mixture[,m-1] <- recodeAlleles(mixture[,m-1], alleleKey)
		mixture[,m] <- recodeAlleles(mixture[,m], alleleKey)
		# flip alleles to all be consistent
		baseline[,m:(m+1)] <- flipHets(baseline[,m:(m+1)])
		mixture[,(m-1):m] <- flipHets(mixture[,(m-1):m])
		
		# recode dropoutProb
		tempDropoutProb <- dropoutProb[dropoutProb[,1] == colnames(baseline)[m],]
		tempDropoutProb[,2] <- recodeAlleles(tempDropoutProb[,2], alleleKey)
		
		# check that no two dropout rates sum to greater than 1
		# this would violate assumption that dropouts are disjoint
		if(checkSums(tempDropoutProb[,3])) stop("two dropout probabilities for locus ", mName, " sum to greater than 1.",
															 " this violates the assumption that dropouts are disjoint.")
		
		# calc per allele error rates
		alleleErr <- matrix(nrow = numAlleles, ncol = numAlleles)
		for(i in 1:(numAlleles - 1)){
			for(j in (i+1):numAlleles){
				tempErr <- 1
				for(c in 1:lenLocus){
					if(substr(alleles[i],c,c) != substr(alleles[j],c,c)){
						tempErr <- tempErr * tempPerSNPerror[c,4] / (tempPerSNPerror[c,3] - 1)
					} else {
						tempErr <- tempErr * (1 - tempPerSNPerror[c,4])
					}
				}
				
				alleleErr[i,j] <- tempErr
				alleleErr[j,i] <- tempErr
			}
		}
		diag(alleleErr) <- prod(1 - tempPerSNPerror[,4])

		# now per genotype error rates

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
		for(i in which(genotypeKey[,1] != genotypeKey[,2])){ #for each true genotype
			for(j in 1:nrow(genotypeKey)){ #for each possible observed genotype
				trueA <- genotypeKey[i,1]
				trueB <- genotypeKey[i,2]
				obsC <- genotypeKey[j,1]
				obsD <- genotypeKey[j,2]
				
				#prob of no dropout
				PnoDropout <- 1 - tempDropoutProb[tempDropoutProb[,2] == trueA,3] - tempDropoutProb[tempDropoutProb[,2] == trueB,3]
				
				if(obsC == obsD && (obsC == trueA || obsC == trueB)){
					# 	P(AA|AB, no dropout) = (1 - P(dropout of A) - P(dropout of B)) * (P(A obs as A) * P(B obs as C)
					PnoDropout <- PnoDropout * alleleErr[trueA,obsC] * alleleErr[trueB,obsD]
				} else if (obsC == obsD && obsC != trueA && obsC != trueB) {
					# P(AA|BC, no dropout) =  (1 - P(dropout of A) - P(dropout of B)) * P(B obs as A) * P(C obs as A)
					PnoDropout <- PnoDropout * alleleErr[trueA,obsC] * alleleErr[trueB,obsD]
				} else {
					# P(CD|AB, no dropout) = (1 - P(dropout of A) - P(dropout of B)) * (P(A obs as C) * P(B obs as D) + P(A obs as D) * P(B obs as C)
					PnoDropout <- PnoDropout * (alleleErr[trueA,obsC] * alleleErr[trueB,obsD] + alleleErr[trueA,obsD] * alleleErr[trueB,obsC])
				}

				trueAA <- which(genotypeKey[,1] == trueA & genotypeKey[,2] == trueA)
				trueBB <- which(genotypeKey[,1] == trueB & genotypeKey[,2] == trueB)
				# P(CD|AB) = P(dropout of A) * P(CD|BB) + P(dropout of B) * P(CD|AA) +  P(CD|AB, no dropout)
				genoErr[i,j] <- (tempDropoutProb[tempDropoutProb[,2] == trueA,3] * genoErr[trueBB,j]) +
					(tempDropoutProb[tempDropoutProb[,2] == trueB,3] * genoErr[trueAA,j]) + 
					PnoDropout
		}
		}
		
		# check sums and issue warning about alleles if pertinent
		if(!isTRUE(all.equal(rowSums(genoErr), rep(1, nrow(genoErr))))) warning("Not all observed genotype probabilities sum to 1 at locus", mName, 
				". This can happen when the number of alleles observed at a SNP in the data do not match the number of SNPs given ",
				"in perSNPerror. If it is simply because some known alleles were not observed in this population, this warning can probably be ignored.")
		
		# paramaters of Dirichlet posterior for allele frequency
		for(p in basePops){
			baselineParams[[as.character(p)]][[mName]] <- countAlleles(baseline[baseline[,1] == p, m:(m+1)], alleleKey[,2]) + (1/nrow(alleleKey))
			names(baselineParams[[as.character(p)]][[mName]]) <- alleleKey[,2] - 1 # -1 for 0 index in c++
		}
		if(useUnsamp){
			for(p in unsamPops){
				unsampledPopsParams[[as.character(p)]][[mName]] <- countAlleles(unsampledPops[unsampledPops[,1] == p, m:(m+1)], alleleKey[,2]) + (1/nrow(alleleKey))
				names(unsampledPopsParams[[as.character(p)]][[mName]]) <- alleleKey[,2] - 1 # -1 for 0 index in c++
			}
		}
		
		#recode genotypes - 0 index for c++
		baseline[,m] <- recodeGenotypes(baseline[,m:(m+1)], genotypeKey) - 1
		mixture[,m-1] <- recodeGenotypes(mixture[,(m-1):m], genotypeKey) - 1
		if(useUnsamp){
			unsampledPops[,m] <- recodeAlleles(unsampledPops[,m], alleleKey)
			unsampledPops[,m+1] <- recodeAlleles(unsampledPops[,m+1], alleleKey)
			unsampledPops[,m:(m+1)] <- flipHets(unsampledPops[,m:(m+1)])
			unsampledPops[,m] <- recodeGenotypes(unsampledPops[,m:(m+1)], genotypeKey) - 1
		}
		
		# create as dataframe
		genotypeKey <- data.frame(genotypeCode = (1:nrow(genotypeKey)) - 1,
										  allele1 = genotypeKey[,1] - 1,
										  allele2 = genotypeKey[,2] - 1, stringsAsFactors = FALSE)
		
		# add row and column names to genoErr for human readability
		rownames(genoErr) <- paste0(genotypeKey$allele1, "/", genotypeKey$allele2)
		colnames(genoErr) <- rownames(genoErr)
		
		# 0 index for alleleKey
		alleleKey[,2] <- alleleKey[,2] - 1
		
		# store results
		genotypeErrorList[[mName]] <- genoErr
		genotypeKeyList[[mName]] <- genotypeKey
		alleleKeyList[[mName]] <- alleleKey

	} # end for each locus
	
	# drop alleles and keep genotype codes
	baseline <- baseline[,c(1, 2, seq(3, (ncol(baseline) - 1), 2))]
	mixture <- mixture[,c(1, seq(2, (ncol(mixture) - 1), 2))]
	if(useUnsamp) {
		unsampledPops <- unsampledPops[,c(1, 2, seq(3, (ncol(baseline) - 1), 2))]
	} else {
		unsampledPopsParams <- NULL
	}
	
	gmaData <- list(
		baseline = baseline,
		mixture = mixture,
		unsampledPops = unsampledPops, # returning unsampled Pops, but don't have a use for it right now
		genotypeErrorRates = genotypeErrorList,
		genotypeKeys = genotypeKeyList,
		alleleKeys = alleleKeyList,
		baselineParams = baselineParams,
		unsampledPopsParams = unsampledPopsParams
	)
	
	return(construct_grandma(gmaData))
}
