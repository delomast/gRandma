# Calculating distances between alleles


#' for microhaplotypes (and snps)
#' @param mName marker name, just used for error messages
#' @param alleleKey data.frame(alleles = alleles, intAlleles = 1:numAlleles, stringsAsFactors = FALSE)
#' @keywords internal
#' @noRd
#' @export

calcDist_hap <- function(mName, alleleKey){

	lenLocus <- unique(nchar(alleleKey[,1]))
	if(length(lenLocus) != 1) stop("not all alleles are the same length for locus ", mName)
	
	distM <- matrix(NA, nrow = nrow(alleleKey), ncol = nrow(alleleKey))
	for(i in 1:nrow(alleleKey)){
		for(j in i:nrow(alleleKey)){
			if(i == j) {
				distM[i,j] <- 0
			} else {
				# distance is number of base pairs
				distM[i,j] <- sum(strsplit(alleleKey[i,1], "")[[1]] !=  strsplit(alleleKey[j,1], "")[[1]])
				distM[j,i] <- distM[i,j]
			}
		}
	}
	return(distM)
}

#' for microsatellites
#' @param mName marker name, just used for error messages
#' @param alleleKey data.frame(alleles = alleles, intAlleles = 1:numAlleles, stringsAsFactors = FALSE)
#' @keywords internal
#' @noRd
#' @export
#' 
calcDist_sat <- function(mName, alleleKey){
	
	alleleKey[,1] <- as.numeric(alleleKey[,1])
	if(any(is.na(alleleKey[,1]))) stop("not all alleles were able to be coerced to numeric for locus ", mName)
	
	distM <- matrix(NA, nrow = nrow(alleleKey), ncol = nrow(alleleKey))
	for(i in 1:nrow(alleleKey)){
		for(j in i:nrow(alleleKey)){
			if(i == j) {
				distM[i,j] <- 0
			} else {
				# distance is difference in lengths
				distM[i,j] <- abs(alleleKey[i,1] - alleleKey[j,1])
				distM[j,i] <- distM[i,j]
			}
		}
	}
	return(distM)
}
