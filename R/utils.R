#' pmf of beta binomial function
#' 
#' @param x number of successes
#' @param n number of trials
#' @param alpha alpha parameter of beta distribution for prob of success
#' @param beta beta parameter of beta distribution for prob of success
#' @param log TRUE to retun the log of the probability

pmfBetaBinomial <- function(x, n, alpha, beta, log = FALSE){
	if(log) return(lchoose(n,x) + lbeta(x + alpha, n - x + beta) - lbeta(alpha, beta))
	return(exp(lchoose(n,x) + lbeta(x + alpha, n - x + beta) - lbeta(alpha, beta)))
}


#' gives conditional probability of an offspring for one biallelic locus assuming
#' genotypes are true
#' 
#' 
#' @param ma genotype of one parent
#' @param pa fenotype of the other parent
#' @param kid genotype of the offspring
#' 

condProbTrio <- function(ma, pa, kid){
	temp <- NA
	if(ma == 1 && pa == 1){
		if(kid == 1){
			temp <- .5
		} else{
			temp <- .25
		}
	} else if (ma == pa){
		if(kid == ma){
			temp <- 1
		} else {
			temp <- 0
		}
	} else if (ma == 1 || pa == 1) {
			if(kid == 1){
				temp <- .5
			} else if (kid == ma || kid == pa) {
				temp <- .5
			} else {
				temp <- 0
			}
	} else if (pa %in% c(0,2) && pa %in% c(0,2)){
		if (kid == 1){
			temp <- 1
		} else {
			temp <- 0
		}
	}
	return(temp)
}

#' recode alleles based on key
#' @param genotypes character vector of genotypes with missing as NA
#' @param key dataframe with first column being alleles (as characters) and second column being
#'   integers representing them 
#' @keywords internal
#' @noRd
recodeAlleles <- function(genotypes, key){
	
	newGenotypes <- genotypes
	for(i in 1:nrow(key)){
		newGenotypes[!is.na(genotypes) & genotypes == key[i,1]] <- key[i,2]
	}
	return(as.numeric(newGenotypes))
}

#' flip alleles represented as integers to have all hets referred to the same (small int first)
#' @param genotypes dataframe with allele1 in column 1 and allele 2 in column 2, missing as NA
#' @keywords internal
#' @noRd
flipHets <- function(genotypes){
	
	newGenotypes <- genotypes
	boolSelect <- (!is.na(genotypes[,1])) & (genotypes[,1] > genotypes[,2])
	newGenotypes[boolSelect,1] <- genotypes[boolSelect,2]
	newGenotypes[boolSelect,2] <- genotypes[boolSelect,1]
	
	return(newGenotypes)
}
