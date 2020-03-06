# R helper functions


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

#' convert NA entries to -9 for c++ to recognize easily
#' #param v vector of entries for which to convert
#' @keywords internal
#' @noRd
convertMissing <- function(v){
	v[is.na(v)] <- -9
	return(v)
}

defaultAlleleDistFunc <- function(x){
	return(1/x)
}

