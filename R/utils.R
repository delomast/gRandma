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

