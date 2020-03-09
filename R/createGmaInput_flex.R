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
#' @param perAlleleError either a constant value representing the per allele error rate (probability of 
#'   observing any allele other than the correct allele)
#'   to use across all loci, or a dataframe with each row representing a locus. Column 1 is the locus name,
#'   and column 2 is the error rate (probability of observing any allele other than the correct allele).
#'   The locus name must be the column name of allele 1 for the corresponding locus in the baseline and mixture dataframes
#' @param dropoutProb either a constant value representing teh probability that any allele in any locus drop out, or
#'   a dataframe with column 1 being locus name, column 2 being allele, and column 3 being dropout probability.
#'   The locus name must be the column name of allele 1 for the corresponding locus in the baseline and mixture dataframes
#' @param alleleDistFunc a function that takes the distance between two non-identical alleles (either number of 
#'   basepair distances for snps/microhaps,
#'   or the numeric distance for microsats) as input, and outputs the weight to give the probabilty of misgenotyping one allele
#'   as the other. The default is weight = 1/x where x is the distance. The probabilty of misgenotyping is the normalized weight
#'   multiplied by the perAlleleError for that locus.
#' @importFrom Rcpp evalCpp
#' @useDynLib gRandma, .registration=TRUE
#' @export

# this version is being written to be more flexible with error models

createGmaInput <- function(baseline, mixture = NULL, unsampledPops = NULL, perAlleleError = NULL, dropoutProb = NULL,
										  markerType = c("microhaps", "snps", "microsats"), 
										  alleleDistFunc = NULL){
	
	mType <- match.arg(markerType)
	if(mType == "snps") mType <- "microhaps" # snps are just microhaps with one position, treating the same for now
	
	if(is.null(mixture)) {
		mixture <- as.data.frame(matrix(0,0,ncol(baseline) - 1))
		colnames(mixture)[seq(2, (ncol(mixture) - 1), 2)] <- colnames(baseline)[seq(3, (ncol(baseline) - 1), 2)]
	}
	
	if(is.null(alleleDistFunc)) alleleDistFunc <- defaultAlleleDistFunc
	
	# type check
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
		warning("Coercing unsampledPops to a dataframe")
		unsampledPops <- as.data.frame(unsampledPops)
	}
	if(any(apply(baseline[,1:2],2, function(x) any(is.na(x))))) stop("NA values are not allowed in baseline columns 1 and 2")
	if(useUnsamp && any(apply(unsampledPops[,1:2],2, function(x) any(is.na(x))))) stop("NA values are not allowed in unsampledPops columns 1 and 2")
	if(any(is.na(mixture[,1]))) stop("NA values are not allowed in mixture column 1")
	
	# first, get names of all markers
	markers <- colnames(baseline)[seq(3, (ncol(baseline) - 1), 2)]
	markersMix <- colnames(mixture)[seq(2, (ncol(mixture) - 1), 2)]
	if(any(markers != markersMix)) stop("error, the loci are either not the same or not in the same order between the baseline and mixture.")
	if(useUnsamp){
		markersUnsamp <- colnames(unsampledPops)[seq(3, (ncol(unsampledPops) - 1), 2)]
		if(markers != markersUnsamp) stop("error, the loci are either not the same or not in the same order between the baseline and unsampledPops.")
	}
	
	# if perAlleleError is a constant
	if(is.numeric(perAlleleError))	{
		if(length(perAlleleError) != 1) stop("perAlleleError must either be a single number or a data frame")
		perAlleleError <- data.frame(locus = markers, error = perAlleleError, stringsAsFactors = FALSE)
	}
	if(any(perAlleleError[,2] < 0 | perAlleleError[,2] > 1)) stop("perAlleleError must be between 0 and 1, inclusive")
	
	if(any(!(markers %in% perAlleleError[,1]))) stop("not all markers have entries in perSNPerror")

	# input error check
	if(is.numeric(dropoutProb) && length(dropoutProb) != 1) stop("dropoutProb must either be a single number or a data frame")
	if(is.numeric(dropoutProb)){
		if(dropoutProb < 0) stop("dropoutProb must be 0 or greater")
	} else {
		if(any(dropoutProb[,3] < 0)) stop("all probabilities in dropoutProb must be 0 or greater")
	}

	# set up storage
	genotypeErrorList <- list() # list of matrices giving P(obs geno | true geno) for each locus
	missingParams <- list() # list of matrices (one each locus) giving parameters of a beta distribution for missing genotypes
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

		mName <- colnames(baseline)[m]
		
		tempperAlleleError <- perAlleleError[perAlleleError[,1] == mName,2]
		if(length(tempperAlleleError) != 1) stop(mName, " should only have one entry in perAlleleError")

		# get all alleles
		alleles <- sort(unique(c(baseline[,m], baseline[,m+1], mixture[,m-1], mixture[,m])))
		if(useUnsamp) alleles <- sort(unique(c(alleles, unsampledPops[,m], unsampledPops[,m+1])))
		alleles <- alleles[!is.na(alleles)]
		
		# check that locus is variable and not all missing
		if(length(alleles) == 0) stop("no genotypes for locus ", mName)
		if(length(alleles) == 1) stop("no variation for locus ", mName)
		
		if(is.numeric(dropoutProb)){
			# if dropoutProb is a constant
			tempDropoutProb <- data.frame(locus = mName, allele = alleles, error = dropoutProb, stringsAsFactors = FALSE)
		} else {
			tempDropoutProb <- dropoutProb[dropoutProb[,1] == mName,]
		}
		if(any(!(alleles %in% tempDropoutProb[,2]))) stop("Not all alleles for ", mName, " are in dropoutProb.")

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
		tempDropoutProb[,2] <- recodeAlleles(tempDropoutProb[,2], alleleKey)
		
		# check that no two dropout rates sum to greater than 1
		# this would violate assumption that dropouts are disjoint
		if(checkSums(tempDropoutProb[,3])) stop("two dropout probabilities for locus ", mName, " sum to greater than 1.",
															 " this violates the assumption that dropouts are disjoint.")
		
		# calculate distances between alleles
		if (mType == "microhaps"){
			alleleDist <- calcDist_hap(mName, alleleKey)
		} else if (mType == "microsats") {
			alleleDist <- calcDist_sat(mName, alleleKey)
		} else {
			stop("Unrecognized mType.")
		}
		
		# calculate proportional error rates for each allele
		# not vectorized to make user supplied functions simpler to write
		# need to skip 0's
		for(i in 1:(numAlleles - 1)){
			for(j in (i + 1):numAlleles){
				alleleDist[i,j] <- alleleDistFunc(alleleDist[i,j])
				alleleDist[j,i] <- alleleDist[i,j]
			}
		}
		
		# calc per allele error rates
		# function of the weights based on distance between the alleles
		alleleErr <- matrix(nrow = numAlleles, ncol = numAlleles)
		# normalize and multiply by per locus error rate
		# NOTE we are assigning results to alleleErr and leaving alleleDist unchanged
		for(i in 1:numAlleles){
			alleleErr[i,] <- (alleleDist[i,] / sum(alleleDist[i,])) * tempperAlleleError
			alleleErr[i,i] <- 1 - tempperAlleleError # prob of not making an error
		}

		# now per genotype error rates
		
		# assign ints to each genotype
		genotypeKey <- matrix(0,0,2)
		for(i in 1:numAlleles){ #build all combinations
			genotypeKey <- rbind(genotypeKey, cbind(i,1:numAlleles))
		}
		#remove repetitive entries
		genotypeKey <- genotypeKey[genotypeKey[,1] <= genotypeKey[,2],]

		# calculate genotype to genotype error probabilities
		# rows are true genotype, columns are observed genotype, entries are probability of observing given true
		genoErr <- matrix(nrow = nrow(genotypeKey), ncol = nrow(genotypeKey))
		# first calculate all for true homozygotes (these are used to then calculate for true heterozygotes with dropout)
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
		
		# check sums and issue warning if pertinent
		# this shouldn't happen... just here for my testing
		if(!isTRUE(all.equal(rowSums(genoErr), rep(1, nrow(genoErr))))) warning("Potential internal error: Not", 
				" all observed genotype probabilities sum to 1 at locus", mName, ".")
		
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
		
		# calculate parameters of beta for whether a genotype is missing or not
		# posterior with Beta(.5,.5) prior
		nMiss <- sum(is.na(baseline[,m]), is.na(mixture[,m-1]))
		nTotal <- nrow(baseline) + nrow(mixture)
		if(useUnsamp){
			nMiss <- nMiss + sum(is.na(unsampledPops[,m]))
			nTotal <- nTotal + nrow(unsampledPops)
		}
		missingParams[[mName]] <- c(nMiss, nTotal - nMiss) + .5
		
		# recode genotypes - 0 index for c++
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
		unsampledPopsParams = unsampledPopsParams,
		missingParams = missingParams
	)
	
	return(construct_grandma(gmaData))
}
