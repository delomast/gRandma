#' Calculate log-likelihood ratios to infer relationships
#' 
#' some more here
#' 
#' and some more here
#' 
#' @param gmaData
#' @param relationship
#' @param crossRecords dataframe or matrix: column 1 is population, column 2 is individual
#'  1 identifier, column 3 is individual 2 identifier
#' 
#' 
#' 
#' @export
inferGrandma <- function(gmaData, relationship = c("ssGP", "test1", "test2"), crossRecords = NULL){
	rel <- match.arg(relationship)
	
	# need to calc llr for all individuals in the mixture compared to all individuals in the baseline
	# this is really just a wrapper for subfunctions for each relationship
	
	# turn pop names and indiv names into ints
	#  speed and avoid headache of dealing with strings
	
	# turn alleles into ints
	# recode genotype key and make sure params are in order
	#### maybe just keep this in gmaData instead of not returning it...

	if(rel == "ssGP"){
		if(is.null(crossRecords)){
			crossRecords <- matrix(0, 0, 0)
		} else{
			crossRecords <- as.matrix(crossRecords)
		}
		
		
		results <- ssGP()
		
	}
	
	# decode everything from ints back into strings
	
}
