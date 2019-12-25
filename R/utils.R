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

#' check that no combination of two different values from one vectors sum to greater than 1
#' used to make sure dropout rates don't violate assumptions
#' returns TRUE is a combination is found that does sum to greater than 1
#' @keywords internal
#' @noRd
checkSums <- function(x){
	if(length(x) < 2) stop("internal error in checkSums")
	for(i in 1:(length(x) - 1)){
		for(j in (i+1):length(x)){
			if(x[i] + x[j] > 1) return(TRUE)
		}
	}
	return(FALSE)
}

#' recode genotypes based on allele calls, does differentiate between "AB" and "BA", so make sure
#' everything is "flipped" appropriately
#' @param genotypes dataframe with allele1 in column 1 and allele 2 in column 2, missing as NA
#' @param key dataframe or matrix with row number being genotype code, allele 1 being column1 and
#'   allele 2 being column 2
#' @keywords internal
#' @noRd
recodeGenotypes <- function(genotypes, key){
	newGenotypes <- genotypes[,1]
	boolSelect <- !is.na(genotypes[,1])
	for(i in 1:nrow(key)){
		newGenotypes[boolSelect & genotypes[,1] == key[i,1] & genotypes[,2] == key[i,2]] <- i
	}
	
	return(newGenotypes)
}

#' count number of each allele to derive Dirichlet posterior
#' @param genotypes dataframe with allele1 in column 1 and allele 2 in column 2, missing as NA
#' @param alleles list of alleles to count
#' @keywords internal
#' @noRd
countAlleles <- function(genotypes, alleles){
	counts <- as.vector(table(c(genotypes[,1], genotypes[,2]), useNA = "no")[as.character(alleles)])
	counts[is.na(counts)] <- 0
	return(counts)
}
